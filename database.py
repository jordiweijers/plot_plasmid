import pandas as pd

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