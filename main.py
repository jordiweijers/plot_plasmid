from pygenomeviz import GenomeViz
import mysql.connector
import pandas as pd
from typing import List, Dict
import os
from time import time
import numpy as np
import argparse

from plot_plasmid.utils import setup_logging, run_command
from plot_plasmid.database import fetch_proteins, fetch_amr_for_proteins, fetch_replicons_from_plasann
from plot_plasmid.blast import build_makeblastdb_command, build_blastp_command
from plot_plasmid.parse import parse_blast_output
from plot_plasmid.plot import plot_contigs

CONTIG_IDS = [1732, 2823, 2591]
OUTPUT_FILE = f"/zfshome/sunam274/compare_plasmids/results/contig_{'_'.join(map(str, CONTIG_IDS))}.svg"
ARO_INDEX_FILE = "/zfshome/sunam274/compare_plasmids/card-data/aro_index.tsv"

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

MIN_IDENTITY = 50.0  

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

def assign_replicons_to_proteins(proteins_df: pd.DataFrame, replicons_df: pd.DataFrame) -> pd.DataFrame:
    """
    Assign replicon from replicons_df to proteins_df.
    A protein is labeld as replicon if:
        - Its ID matches rpelicons_df.kes2021_id
        - Or it is the closest protein to a replicon on the same contig, that does not have a direct match with replicons_df.kes2021_id. This is to account for cases where the replicon prediction may not have been perfectly accurate in terms of start and end positions.
    Args:
        proteins_df (pd.DataFrame): A DataFrame containing the proteins to assign replicons to. Must have columns 'id', 'contig_ID', 'start', and 'end'.
        replicons_df (pd.DataFrame): A DataFrame containing the replicons to assign. Must have columns 'kes2021_id', 'contig_ID', 'start', and 'end'.
    Returns:
        pd.DataFrame: A DataFrame containing the proteins with an additional column 'is_replicon' indicating whether each protein is a replicon or not.
    """
    proteins_df = proteins_df.copy()
    proteins_df['replicon'] = pd.NA
    if replicons_df.empty:
        return proteins_df
    # ---------- First assign replicons based on direct ID match ----------
    direct_replicons = replicons_df[replicons_df['kes2021_id'].notna()]
    for _, rep in direct_replicons.iterrows():
        mask = proteins_df['id'] == rep['kes2021_id']
        proteins_df.loc[mask, 'replicon'] = rep['gene_name']
    # ---------- Coordinate assignment for reamining replicons ----------
    coord_replicons = replicons_df[replicons_df['kes2021_id'].isna()]
    for _, replicon in coord_replicons.iterrows():
        contig_proteins = proteins_df[proteins_df['contig_ID'] == replicon['contig_ID']]
        if contig_proteins.empty:
            continue
        contig_proteins['distance'] = np.minimum(
            abs(contig_proteins['start'] - replicon['start']),
            abs(contig_proteins['end'] - replicon['end'])
        )
        closest_idx = contig_proteins['distance'].idxmin()
        proteins_df.loc[closest_idx, 'replicon'] = replicon['gene_name']
    return proteins_df

def parse_arguments():
    """
    Parse command-line arguments.
    Returns:
        argparse.Namespace: The parsed arguments.
    """
    parser = argparse.ArgumentParser(description="Plot plasmid contigs with PyGenomeViz")
    parser.add_argument(
        "-c", "--contig_ids", type=int, nargs="+", default=CONTIG_IDS, help="List of contig IDs to plot (space-separated)"
    )
    parser.add_argument(
        "-o", "--output_file", type=str, default=OUTPUT_FILE, help="Output file path for the plot (e.g. contigs.svg)"
    )
    return parser.parse_args()

def main():
    args = parse_arguments()
    CONTIG_IDS = args.contig_ids
    OUTPUT_FILE = args.output_file
    logger = setup_logging(os.path.join(LOG_DIR, f"plot_plasmid_{int(time())}.log"))
    conn = mysql.connector.connect(**DB_CONFIG)
    aro_index_df = load_aro_index(ARO_INDEX_FILE)
    all_contigs_df = pd.DataFrame()
    faa_paths = {}
    db_prefixes = {}
    for contig_id in CONTIG_IDS:

        # ---------- Get contig data ----------
        proteins_df = fetch_proteins(contig_id, conn)
        amr_df = fetch_amr_for_proteins(proteins_df['id'].tolist(), conn)
        if not amr_df.empty:
            amr_df['ARO'] = amr_df['ARO'].astype(int)
            contig_df = proteins_df.merge(amr_df, how='left', on='id')
            contig_df = contig_df.merge(aro_index_df, how='left', left_on='ARO', right_on='ARO Accession')
        else:
            contig_df = proteins_df.copy()
        replicons_df = fetch_replicons_from_plasann(contig_id, conn)
        contig_df = assign_replicons_to_proteins(contig_df, replicons_df)
        all_contigs_df = pd.concat([all_contigs_df, contig_df], ignore_index=True)

        # ---------- Write FASTA ----------
        faa_path = os.path.join(FASTA_DIR, f"{contig_id}.faa")
        if os.path.exists(faa_path):
            logger.info(f"Using existing FASTA for contig {contig_id} at {faa_path}")
        else:
            logger.info(f"Writing FASTA for contig {contig_id} to {faa_path}")
            sequences = {f"{row['id']}": row['sequence'] for _, row in contig_df.iterrows()}
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
    adjacent_pairs = [(CONTIG_IDS[i], CONTIG_IDS[i+1]) for i in range(len(CONTIG_IDS)-1)]
    blast_results = {}
    for query, subject in adjacent_pairs:
        faa = faa_paths[query]
        db_prefix = db_prefixes[subject]
        output_file = os.path.join(BLAST_OUTPUT_DIR, f"{query}_vs_{subject}.tsv")
        if os.path.exists(output_file):
            logger.info(f"Using existing BLASTP output for {query} vs {subject} at {output_file}")
        else:
            logger.info(f"Running BLASTP for {query} vs {subject}, outputting to {output_file}")
            cmd = build_blastp_command(faa, db_prefix, output_file, BLASTP_PARAMS)
            run_command(cmd, os.path.join(LOG_DIR, f"blastp_{query}_vs_{subject}_{int(time())}.log"))
        blast_results[(query, subject)] = parse_blast_output(output_file)

    # ---------- Plot contigs ----------
    plot_contigs(all_contigs_df, CONTIG_IDS, blast_results, MIN_IDENTITY,OUTPUT_FILE)

if __name__ == "__main__":
    main()