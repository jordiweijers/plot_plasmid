from pygenomeviz import GenomeViz
import mysql.connector
import pandas as pd
from typing import List, Dict, Tuple
import os
from time import time
import numpy as np
import argparse
import itertools

from plot_plasmid.utils import setup_logging, run_command, run_in_parallel
from plot_plasmid.database import fetch_proteins, fetch_amr_for_proteins, fetch_mcl_for_proteins, fetch_replicons_for_contig, fetch_pseudogenes, fetch_mcl_for_pseudogenes
from plot_plasmid.blast import build_makeblastdb_command, build_blastp_command
from plot_plasmid.parse import parse_blast_output
from plot_plasmid.plot import plot_contigs
from plot_plasmid.order import canonical_pair, swap_query_subject, score_pair, order_by_clustering

CONTIG_IDS = [3889, 9353, 10599]
OUTPUT_FILE = f"/zfshome/sunam274/compare_plasmids/results/plots/{'_'.join(map(str, CONTIG_IDS))}.svg"
ARO_INDEX_FILE = "/zfshome/sunam274/compare_plasmids/card-data/aro_index.tsv"
REPLICON_MCL_FILE = "/zfshome/sunam274/compare_plasmids/results/replicon_mcl.csv"
CONJUGATION_MCL_FILE = "/zfshome/sunam274/compare_plasmids/results/conjugation_mcl.csv"

DB_CONFIG = {
    "host": "localhost",
    "user": "gmg",
    "database": "kes2021"
}

MAKEBLASTDB_PARAMS = {
    "dbtype": "prot",
}

BLASTP_PARAMS = {
    "outfmt": 6,
    "evalue": 1e-4,
    "word_size": 3,
    "max_target_seqs": 1000000,
    "num_threads": 1
}

MIN_IDENTITY = 60.0  

LOG_DIR = "/zfshome/sunam274/compare_plasmids/results/plotting_data/logs"
os.makedirs(LOG_DIR, exist_ok=True)
FASTA_DIR = "/zfshome/sunam274/compare_plasmids/results/plotting_data/faa"
os.makedirs(FASTA_DIR, exist_ok=True)
BLAST_DB_DIR = "/zfshome/sunam274/compare_plasmids/results/plotting_data/blast_db"
os.makedirs(BLAST_DB_DIR, exist_ok=True)
BLAST_OUTPUT_DIR = "/zfshome/sunam274/compare_plasmids/results/plotting_data/blast_output"
os.makedirs(BLAST_OUTPUT_DIR, exist_ok=True)

def load_aro_index(file_path: str) -> pd.DataFrame:
    """
    Load the ARO index from a TSV file and return it as a pandas DataFrame.
    Args:
        file_path (str): The path to the ARO index TSV file.
    Returns:
        pd.DataFrame: A DataFrame containing the ARO index.
    """
    df = pd.read_csv(file_path, sep='\t')
    df['ARO Accession'] = df['ARO Accession'].str.replace("ARO:", "").astype(int)
    return df

def load_mcl_clusters(file_path: str) -> pd.DataFrame:
    """
    Load MCL cluster data (e.g. replicon or conjugation clusters) from a CSV file and return it as a
    pandas DataFrame.
    Args:
        file_path (str): The path to the MCL clusters CSV file.
    Returns:
        pd.DataFrame: A DataFrame containing the MCL cluster data.
    """
    df = pd.read_csv(file_path)
    df['mcl_id'] = df['mcl_id'].astype(int)
    return df

def write_fasta(sequences: Dict[str, str], output_file: str):
    """
    Write a dictionary of sequences and corresponding headers to a FASTA file.
    Args:
        sequences (Dict[str, str]): A dictionary where keys are headers and values are sequences.
        output_file (str): The path to the output FASTA file.
    """
    with open(output_file, 'w') as f:
        for header, seq in sequences.items():
            f.write(f">{header}\n{seq}\n")

def _run_single_blastp(query: int, subject: int, faa: str, db_prefix: str, output_file: str, log_path: str) -> Tuple[int, int]:
    """
    Run a single BLASTP command for one (query, subject) pair. Meant to be dispatched via
    `run_in_parallel`, so it must be a top-level, picklable function.
    Args:
        query (int): The query contig ID.
        subject (int): The subject contig ID.
        faa (str): The path to the query protein FASTA file.
        db_prefix (str): The BLAST database prefix for the subject.
        output_file (str): The path to write the BLASTP output to.
        log_path (str): The path to write the BLASTP log to.
    Returns:
        Tuple[int, int]: The (query, subject) pair that was BLASTed, once complete.
    """
    cmd = build_blastp_command(faa, db_prefix, output_file, BLASTP_PARAMS)
    run_command(cmd, log_path)
    return (query, subject)

def run_blastp_pairs(pairs: List[Tuple[int, int]], faa_paths: Dict[int, str], db_prefixes: Dict[int, str], logger, cpus: int = 1) -> Dict[Tuple[int, int], pd.DataFrame]:
    """
    Run (or reuse cached) BLASTP for each (query, subject) contig pair and parse the results.
    Args:
        pairs (List[Tuple[int, int]]): A list of (query_contig_id, subject_contig_id) pairs to BLAST.
        faa_paths (Dict[int, str]): A dictionary mapping contig_id to its protein FASTA path.
        db_prefixes (Dict[int, str]): A dictionary mapping contig_id to its BLAST database prefix.
        logger: A logger instance for logging progress.
        cpus (int): Number of CPUs to use to run BLASTP pairs in parallel.
    Returns:
        Dict[Tuple[int, int], pd.DataFrame]: A dictionary mapping each (query, subject) pair to its parsed BLASTP results.
    """
    output_files = {}
    pending_args = []
    for query, subject in pairs:
        output_file = os.path.join(BLAST_OUTPUT_DIR, f"{query}_vs_{subject}.tsv")
        output_files[(query, subject)] = output_file
        if os.path.exists(output_file):
            logger.info(f"Using existing BLASTP output for {query} vs {subject} at {output_file}")
        else:
            logger.info(f"Queuing BLASTP for {query} vs {subject}, outputting to {output_file}")
            log_path = os.path.join(LOG_DIR, f"blastp_{query}_vs_{subject}_{int(time())}.log")
            pending_args.append((query, subject, faa_paths[query], db_prefixes[subject], output_file, log_path))

    if pending_args:
        run_in_parallel(_run_single_blastp, pending_args, cpus)

    return {pair: parse_blast_output(output_file) for pair, output_file in output_files.items()}

def parse_arguments():
    parser = argparse.ArgumentParser(description="Plot plasmid contigs with PyGenomeViz")
    input_group = parser.add_mutually_exclusive_group()
    input_group.add_argument(
        "-c", "--contig_ids", type=int, nargs="+", help="List of contig IDs to plot (space-separated)"
    )
    input_group.add_argument(
        "-f", "--contig_file", type=str, help="Path to a file with one contig ID per line, optionally followed by a tab and a category for label coloring"
    )
    parser.add_argument(
        "-o", "--output_file", type=str, default=None, help="Output file path for the plot (e.g. contigs.svg)"
    )
    parser.add_argument(
        "-s", "--order-by-similarity", action="store_true",
        help="Reorder contigs by all-vs-all BLASTP similarity clustering instead of the input order"
    )
    parser.add_argument(
        "-p", "--cpus", type=int, default=1, help="Number of CPUs to use for running BLASTP pairs in parallel"
    )
    return parser.parse_args()

def run_plot(
        contig_ids: List[int],
        output_file: str,
        order_by_similarity: bool = False,
        contig_categories: Dict[int, str] = None,
        cpus: int = 1,
    ):
    """
    Fetch DB data for the given contigs, run/reuse BLASTP as needed, and plot them.

    This is `main()`'s body minus argument parsing, so it can be called directly (e.g. from
    another script) with explicit arguments instead of going through argparse/`sys.argv`.
    Args:
        contig_ids (List[int]): A list of contig IDs to plot.
        output_file (str): The path to the output file.
        order_by_similarity (bool): Whether to reorder contigs by all-vs-all BLASTP similarity clustering.
        contig_categories (Dict[int, str]): A dictionary mapping contig_id to a category label for track-label coloring.
        cpus (int): Number of CPUs to use for running BLASTP pairs in parallel.
    Returns:
        Tuple[GenomeViz, Dict[int, FeatureTrack]]: as returned by `plot_contigs`.
    """
    CONTIG_IDS = contig_ids
    contig_categories = contig_categories or {}
    OUTPUT_FILE = output_file
    logger = setup_logging(os.path.join(LOG_DIR, f"plot_plasmid_{int(time())}.log"))
    conn = mysql.connector.connect(**DB_CONFIG)
    aro_index_df = load_aro_index(ARO_INDEX_FILE)
    replicon_df = load_mcl_clusters(REPLICON_MCL_FILE)
    conjugation_df = load_mcl_clusters(CONJUGATION_MCL_FILE)
    all_contigs_df = pd.DataFrame()
    plasann_replicons = {}
    faa_paths = {}
    db_prefixes = {}
    for contig_id in CONTIG_IDS:

        # ---------- Get contig data ----------
        try:
            proteins_df = fetch_proteins(contig_id, conn)
        except Exception as e:
            raise ValueError(f"Failed to fetch proteins for contig {contig_id}: {e}.")
        proteins_df = proteins_df.copy()
        proteins_df['type'] = 'protein'
        try:
            pseudogenes_df = fetch_pseudogenes(contig_id, conn)
        except Exception as e:
            logger.error(f"Error fetching pseudogenes for contig {contig_id}: {e}")
            pseudogenes_df = pd.DataFrame()
        if not pseudogenes_df.empty:
            pseudogenes_df = pseudogenes_df.copy()
            pseudogenes_df['type'] = 'pseudo'
            pseudogenes_df['pseudo_id'] = pseudogenes_df['pseudo_id'].astype(str)
            pseudogenes_df['id'] = 'pseudo_' + pseudogenes_df['pseudo_id']
            contig_df = pd.concat([proteins_df, pseudogenes_df], ignore_index=True)
        else:
            contig_df = proteins_df.copy()
        try:
            plasann_df = fetch_replicons_for_contig(contig_id, conn)
        except Exception as e:
            logger.error(f"Error fetching PlasAnn replicon data for contig {contig_id}: {e}")
            plasann_df = pd.DataFrame()
        plasann_replicons[contig_id] = plasann_df

        
        # ---------- Get protein annotations ----------
        try:
            amr_df = fetch_amr_for_proteins(proteins_df['id'].tolist(), conn)
        except Exception as e:
            logger.error(f"Error fetching AMR data for contig {contig_id}: {e}")
            amr_df = pd.DataFrame()
        if not amr_df.empty:
            amr_df['ARO'] = amr_df['ARO'].astype(int)
            contig_df = contig_df.merge(amr_df, how='left', on='id')
            contig_df = contig_df.merge(aro_index_df, how='left', left_on='ARO', right_on='ARO Accession')
        try:
            mcl_df = fetch_mcl_for_proteins(proteins_df['id'].tolist(), conn)
        except Exception as e:
            logger.error(f"Error fetching MCL data for contig {contig_id}: {e}")
            mcl_df = pd.DataFrame()
        if not mcl_df.empty:
            mcl_df['mcl_id'] = mcl_df['mcl_id'].astype(int)
        try:
            pmcl_df = fetch_mcl_for_pseudogenes(pseudogenes_df['id'].tolist(), conn)
        except Exception as e:
            logger.error(f"Error fetching MCL data for pseudogenes in contig {contig_id}: {e}")
            pmcl_df = pd.DataFrame()
        if not pmcl_df.empty and not pseudogenes_df.empty:
            pmcl_df['mcl_id'] = pmcl_df['mcl_id'].astype(int)
            pmcl_df['id'] = 'pseudo_' + pmcl_df['pseudo_id'].astype(str)
        combined_mcl_df = pd.concat([mcl_df, pmcl_df], ignore_index=True)
        if not combined_mcl_df.empty:
            replicon_matches = combined_mcl_df.merge(replicon_df, how='inner', on='mcl_id')
            replicon_matches = replicon_matches.rename(columns={'cluster_name': 'replicon'})
            contig_df = contig_df.merge(replicon_matches[['id', 'mcl_id', 'replicon']], how='left', on='id')

            conjugation_matches = combined_mcl_df.merge(conjugation_df, how='inner', on='mcl_id')
            conjugation_matches = conjugation_matches.rename(columns={'cluster_name': 'conjugation'})
            contig_df = contig_df.merge(conjugation_matches[['id', 'conjugation']], how='left', on='id')
        all_contigs_df = pd.concat([all_contigs_df, contig_df], ignore_index=True)

        # ---------- Write FASTA ----------
        faa_path = os.path.join(FASTA_DIR, f"{contig_id}.faa")
        if os.path.exists(faa_path):
            logger.info(f"Using existing FASTA for contig {contig_id} at {faa_path}")
        else:
            logger.info(f"Writing FASTA for contig {contig_id} to {faa_path}")
            sequences = {f"{row['id']}": row['sequence'] for _, row in contig_df.iterrows() if row['type'] == 'protein'}
            write_fasta(sequences, faa_path)
        faa_paths[contig_id] = faa_path

        # ---------- Create BLAST database ----------
        db_prefix = os.path.join(BLAST_DB_DIR, f"{contig_id}")
        check_path = os.path.join(db_prefix + ".pin")
        if os.path.exists(check_path):
            logger.info(f"Using existing BLAST database for contig {contig_id} at {db_prefix}")
        else:
            logger.info(f"Creating BLAST database for contig {contig_id} at {db_prefix}")
            cmd, _ = build_makeblastdb_command(faa_path, BLAST_DB_DIR, MAKEBLASTDB_PARAMS)
            run_command(cmd, os.path.join(LOG_DIR, f"makeblastdb_{contig_id}_{int(time())}.log"))
        db_prefixes[contig_id] = db_prefix
    conn.close()

    # ---------- Run BLASTP ----------
    if order_by_similarity:
        all_pairs = list(itertools.combinations(sorted(set(CONTIG_IDS)), 2))
        blast_results = run_blastp_pairs(all_pairs, faa_paths, db_prefixes, logger, cpus)

        pair_scores = {}
        for pair, df in blast_results.items():
            try:
                pair_scores[pair] = score_pair(df, MIN_IDENTITY)
            except ValueError as e:
                logger.warning(f"No BLASTP hits above identity threshold for pair {pair}: {e}")
                pair_scores[pair] = 0.0

        CONTIG_IDS = order_by_clustering(CONTIG_IDS, pair_scores)
        logger.info(f"Reordered contigs by similarity clustering: {CONTIG_IDS}")

        adjacent_pairs = [(CONTIG_IDS[i], CONTIG_IDS[i+1]) for i in range(len(CONTIG_IDS)-1)]
        adjacent_blast_results = {}
        for query, subject in adjacent_pairs:
            canon = canonical_pair(query, subject)
            df = blast_results[canon]
            adjacent_blast_results[(query, subject)] = df if canon == (query, subject) else swap_query_subject(df)
        blast_results = adjacent_blast_results
    else:
        adjacent_pairs = [(CONTIG_IDS[i], CONTIG_IDS[i+1]) for i in range(len(CONTIG_IDS)-1)]
        blast_results = run_blastp_pairs(adjacent_pairs, faa_paths, db_prefixes, logger, cpus)

    # ---------- Plot contigs ----------
    return plot_contigs(all_contigs_df, CONTIG_IDS, blast_results, plasann_replicons, MIN_IDENTITY, OUTPUT_FILE, contig_categories)

def main():
    args = parse_arguments()
    contig_categories: Dict[int, str] = {}
    if args.contig_file:
        CONTIG_IDS = []
        with open(args.contig_file) as fh:
            for line in fh:
                line = line.rstrip("\n")
                if not line.strip():
                    continue
                parts = line.split("\t")
                contig_id = int(parts[0].strip())
                CONTIG_IDS.append(contig_id)
                if len(parts) > 1 and parts[1].strip():
                    contig_categories[contig_id] = parts[1].strip()
    elif args.contig_ids:
        CONTIG_IDS = args.contig_ids
    # else: fall back to module-level CONTIG_IDS default
    if not args.output_file:
        if args.contig_file:
            base = os.path.splitext(os.path.basename(args.contig_file))[0]
        else:
            base = '_'.join(map(str, CONTIG_IDS))
        OUTPUT_FILE = f"/zfshome/sunam274/compare_plasmids/results/plots/{base}.svg"
    else:
        OUTPUT_FILE = args.output_file
    run_plot(CONTIG_IDS, OUTPUT_FILE, args.order_by_similarity, contig_categories, args.cpus)

if __name__ == "__main__":
    main()