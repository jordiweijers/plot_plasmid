import pandas as pd
import numpy as np
from typing import List

def fetch_proteins(contig_id: int, conn) -> pd.DataFrame:
    """
    Fetch all proteins for a given contig_ID from the MySQL database and return them as a pandas DataFrame.
    Args:
        contig_id (int): The ID of the contig to fetch proteins for.
        conn: A MySQL database connection object.
    Returns:
        pd.DataFrame: A DataFrame containing the proteins for the given contig_ID.
    """
    cursor = conn.cursor(dictionary=True)
    query = """
    SELECT *
    FROM proteins
    WHERE contig_ID = %s
    order by start
    """
    cursor.execute(query,(int(contig_id),))
    rows = cursor.fetchall()
    df = pd.DataFrame(rows)
    return df

def fetch_amr_for_proteins(protein_ids: List[str], conn) -> pd.DataFrame:
    """
    Fetch AMR information for a list of protein IDs from the MySQL database and return it as a pandas DataFrame.
    Args:
        protein_ids (List[str]): A list of protein IDs to fetch AMR information for.
        conn: A MySQL database connection object.
    Returns:
        pd.DataFrame: A DataFrame containing the AMR information for the given protein IDs.
    """
    if not protein_ids:
        raise ValueError("protein_ids list cannot be empty.")
    cursor = conn.cursor(dictionary=True)
    placeholders = ', '.join(['%s'] * len(protein_ids))
    query = f"""
    SELECT id, ARO 
    FROM amr
    WHERE id IN ({placeholders})
    """
    cursor.execute(query, protein_ids)
    rows = cursor.fetchall()
    return pd.DataFrame(rows)

def fetch_replicons_from_plasann(contig_id: str, conn) -> pd.DataFrame:
    """
    Fetch replicon information for a given contig_ID from the MySQL database and return it as a pandas DataFrame.
    Args:
        contig_id (str): The ID of the contig to fetch replicon information for.
        conn: A MySQL database connection object.
    Returns:
        pd.DataFrame: A DataFrame containing the replicon information for the given contig_ID.
    """
    cursor = conn.cursor(dictionary=True)
    query = """
    SELECT id, kes2021_id, contig_ID, start, end, product, gene_name
    FROM plasann_kes2021.proteins
    WHERE contig_ID = %s
        AND product = "Predicted replicon"
    """
    cursor.execute(query,(int(contig_id),))
    rows = cursor.fetchall()
    if not rows:
        return pd.DataFrame(columns=['id', 'kes2021_id', 'contig_ID', 'start', 'end', 'product', 'gene_name'])
    df = pd.DataFrame(rows)
    df['start'] = df['start'].astype(int)
    df['end'] = df['end'].astype(int)
    return df
