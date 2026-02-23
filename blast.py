import os
from typing import List, Tuple, Dict, Any

def build_makeblastdb_command(fasta_file: str, db_dir: str, params: Dict[str, Any]) -> Tuple[List[str], str]:
    """
    Build the makeblastdb command for a given FASTA file.
    Args:
        fasta_file (str): The path towards the FASTA file.
        db_dir (str): The path towards a directory where to store the BLAST database.
        params (dict): A dictionary containing the parameters to use.
    Returns:
        command (list[str]): The command to run.
        db_prefix (str): The prefix for the BLAST database.
    """
    fasta_name = os.path.splitext(os.path.basename(fasta_file))[0]
    os.makedirs(db_dir, exist_ok=True)
    db_prefix = os.path.join(db_dir, fasta_name)
    cmd = [
        "makeblastdb",
        "-in", fasta_file,
        "-dbtype", str(params["dbtype"]),
        "-out", db_prefix
    ]
    return cmd, db_prefix

def build_blastp_command(query_fasta: str, sub_db_prefix: str, output_file: str, params: Dict[str, Any]) -> List[str]:
    """
    Build the blastp command for a given query fasta and subject blastdb.
    Args:
        query_fasta (str): The path towards the query FASTA file.
        sub_db_prefix (str): The path towards the subject BLAST database prefix.
        output_file (str): The path where to store the output file.
        
    Returns:
        command (list[str]): The command to run.
    """
    cmd = [
        "blastp",
        "-query", query_fasta,
        "-db", sub_db_prefix,
        "-out", output_file,
        "-outfmt", str(params["outfmt"]),
        "-evalue", str(params["evalue"]),
        "-word_size", str(params["word_size"]),
        "-max_target_seqs", str(params["max_target_seqs"]),
        "-num_threads", str(params["num_threads"])
    ]
    return cmd