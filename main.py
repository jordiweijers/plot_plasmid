from pygenomeviz import GenomeViz
import mysql.connector
import pandas as pd
from typing import List, Dict
import os
from time import time
import numpy as np
import argparse

from plot_plasmid.utils import setup_logging, run_command
from plot_plasmid.database import fetch_proteins, fetch_amr_for_proteins, fetch_mcl_for_proteins, fetch_replicons_for_contig
from plot_plasmid.blast import build_makeblastdb_command, build_blastp_command
from plot_plasmid.parse import parse_blast_output
from plot_plasmid.plot import plot_contigs

CONTIG_IDS = [1732, 2823, 2591]
OUTPUT_FILE = f"/zfshome/sunam274/compare_plasmids/results/plots/{'_'.join(map(str, CONTIG_IDS))}.svg"
ARO_INDEX_FILE = "/zfshome/sunam274/compare_plasmids/card-data/aro_index.tsv"
REPLICON_MCL_FILE = "/zfshome/sunam274/compare_plasmids/results/replicon_mcl.csv"

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

def load_replicon_mcl(file_path: str) -> pd.DataFrame:
    """
    Load the replicon MCL data from a CSV file and return it as a pandas DataFrame.
    Args:
        file_path (str): The path to the replicon MCL CSV file.
    Returns:
        pd.DataFrame: A DataFrame containing the replicon MCL data.
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
        "-o", "--output_file", type=str, default=None, help="Output file path for the plot (e.g. contigs.svg)"
    )
    return parser.parse_args()

def main():
    args = parse_arguments()
    CONTIG_IDS = args.contig_ids
    if not args.output_file:
        OUTPUT_FILE = f"/zfshome/sunam274/compare_plasmids/results/plots/{'_'.join(map(str, CONTIG_IDS))}.svg"
    else:
        OUTPUT_FILE = args.output_file
    logger = setup_logging(os.path.join(LOG_DIR, f"plot_plasmid_{int(time())}.log"))
    conn = mysql.connector.connect(**DB_CONFIG)
    aro_index_df = load_aro_index(ARO_INDEX_FILE)
    replicon_df = load_replicon_mcl(REPLICON_MCL_FILE)
    all_contigs_df = pd.DataFrame()
    plasann_replicons = {}
    faa_paths = {}
    db_prefixes = {}
    for contig_id in CONTIG_IDS:

        # ---------- Get contig data ----------
        proteins_df = fetch_proteins(contig_id, conn)
        try:
            amr_df = fetch_amr_for_proteins(proteins_df['id'].tolist(), conn)
        except Exception as e:
            logger.error(f"Error fetching AMR data for contig {contig_id}: {e}")
            amr_df = pd.DataFrame()
        if not amr_df.empty:
            amr_df['ARO'] = amr_df['ARO'].astype(int)
            contig_df = proteins_df.merge(amr_df, how='left', on='id')
            contig_df = contig_df.merge(aro_index_df, how='left', left_on='ARO', right_on='ARO Accession')
        else:
            contig_df = proteins_df.copy()
        try:
            mcl_df = fetch_mcl_for_proteins(proteins_df['id'].tolist(), conn)
        except Exception as e:
            logger.error(f"Error fetching MCL data for contig {contig_id}: {e}")
            mcl_df = pd.DataFrame()
        if not mcl_df.empty:
            mcl_df['mcl_id'] = mcl_df['mcl_id'].astype(int)
            mcl_rep_df = replicon_df.merge(mcl_df, how='inner', on='mcl_id')
            mcl_rep_df = mcl_rep_df.rename(columns={'cluster_name': 'replicon'})
            contig_df = contig_df.merge(mcl_rep_df[['id', 'replicon', 'mcl_id']], how='left', on='id')
        else:
            contig_df = contig_df.copy()
        try:
            plasann_df = fetch_replicons_for_contig(contig_id, conn)
        except Exception as e:
            logger.error(f"Error fetching PlasAnn replicon data for contig {contig_id}: {e}")
            plasann_df = pd.DataFrame()
        plasann_replicons[contig_id] = plasann_df
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
    plot_contigs(all_contigs_df, CONTIG_IDS, blast_results, plasann_replicons, MIN_IDENTITY, OUTPUT_FILE)

if __name__ == "__main__":
    main()