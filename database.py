import mysql.connector
import pandas as pd
from typing import Dict, List
from dataclasses import dataclass

from plot_plasmid.models import Contig, Feature
from plot_plasmid.config import DatabaseSchema, ContigSource, FeatureSource


class ContigLoader:
    """
    Loads contigs and their features from a MySQL database, using the table and column names in its DatabaseSchema.
    """
    def __init__(self, database: str, schema: DatabaseSchema, host: str, user: str):
        self.database = database
        self.schema = schema
        self.conn = mysql.connector.connect(host=host, user=user, database=database)

    def __enter__(self) -> "ContigLoader":
        return self

    def __exit__(self, *exc) -> None:
        self.close()

    def close(self) -> None:
        self.conn.close()

    def _query(self, query: str, params: tuple) -> List[dict]:
        cursor = self.conn.cursor(dictionary=True)
        cursor.execute(query, params)
        rows = cursor.fetchall()
        cursor.close()
        return rows

    def fetch_contig(self, contig_id: str) -> Contig:
        source = self.schema.contigs
        rows = self._query(
            f"SELECT {source.length_column} AS length FROM {source.table} WHERE {source.id_column} = %s",
            (contig_id,),
        )
        if not rows:
            raise ValueError(f"Contig '{contig_id}' not found in database '{self.database}'.")
        return Contig(
            id=str(contig_id),
            length=int(rows[0]["length"]),
            features=self.fetch_features(contig_id),
        )
    
    def fetch_features(self, contig_id: str) -> List[Feature]:
        features = []
        for source in self.schema.features:
            rows = self._query(
                f"""
                SELECT f.{source.id_column} AS id, f.start, f.end, f.strand,
                       ({source.is_pseudo}) AS is_pseudo,
                       m.{source.mcl_family_column} AS protein_family
                FROM {source.table} f
                LEFT JOIN {source.mcl_table} m ON m.{source.mcl_id_column} = f.{source.id_column}
                WHERE f.{source.contig_column} = %s
                """,
                (contig_id,),
            )
            aro_by_id = self.fetch_aro(source, contig_id)
            for row in rows:
                feature_id = str(row["id"])
                features.append(Feature(
                    id=feature_id,
                    type="pseudo" if row["is_pseudo"] else "protein",
                    start=int(row["start"]),
                    end=int(row["end"]),
                    strand=row["strand"],
                    protein_family=str(row["protein_family"]) if row["protein_family"] is not None else None,
                    aro=aro_by_id.get(feature_id),
                ))
        features.sort(key=lambda feature: feature.start)
        return features

    def fetch_aro(self, source: FeatureSource, contig_id: str) -> Dict[str, str]:
        if source.amr_table is None:
            return {}
        rows = self._query(
            f"""
            SELECT a.{source.amr_id_column} AS id, a.{source.amr_aro_column} AS aro
            FROM {source.amr_table} a
            JOIN {source.table} f ON f.{source.id_column} = a.{source.amr_id_column}
            WHERE f.{source.contig_column} = %s
            """,
            (contig_id,),
        )
        return {str(row["id"]): str(row["aro"]) for row in rows}


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

def fetch_pseudogenes(contig_id: int, conn) -> pd.DataFrame:
    """
    Fetch all pseudogenes for a given contig_ID from the MySQL database and return them as a pandas DataFrame.
    Args:
        contig_id (int): The ID of the contig to fetch pseudogenes for.
        conn: A MySQL database connection object.
    Returns:
        pd.DataFrame: A DataFrame containing the pseudogenes for the given contig_ID.
    """
    cursor = conn.cursor(dictionary=True)
    query = """
    SELECT *
    FROM pseudogenes
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
    if not  rows:
        raise ValueError("No AMR information found for the provided protein IDs.")
    return pd.DataFrame(rows)

def fetch_mcl_for_proteins(protein_ids: List[str], conn) -> pd.DataFrame:
    """
    Fetch MCL cluster information for a list of protein IDs from the MySQL database and return it as a pandas DataFrame.
    Args:
        protein_ids (List[str]): A list of protein IDs to fetch MCL cluster information for.
        conn: A MySQL database connection object.
    Returns:
        pd.DataFrame: A DataFrame containing the MCL cluster information for the given protein IDs.
    """
    if not protein_ids:
        raise ValueError("protein_ids list cannot be empty.")
    cursor = conn.cursor(dictionary=True)
    placeholders = ', '.join(['%s'] * len(protein_ids))
    query = f"""
    SELECT id, clust AS mcl_id
    FROM mcl30
    WHERE id IN ({placeholders})
    """
    cursor.execute(query, protein_ids)
    rows = cursor.fetchall()
    if not rows:
        raise ValueError("No MCL cluster information found for the provided protein IDs.")
    return pd.DataFrame(rows)

def fetch_mcl_for_pseudogenes(pseudo_ids: List[str], conn) -> pd.DataFrame:
    """
    Fetch MCL cluster information for a list of pseudogene IDs from the MySQL database and return it as a pandas DataFrame.
    Args:
        pseudo_ids (List[str]): A list of pseudogene IDs to fetch MCL cluster information for.
        conn: A MySQL database connection object.
    Returns:
        pd.DataFrame: A DataFrame containing the MCL cluster information for the given pseudogene IDs.
    """
    if not pseudo_ids:
        raise ValueError("pseudo_ids list cannot be empty.")
    cursor = conn.cursor(dictionary=True)
    placeholders = ', '.join(['%s'] * len(pseudo_ids))
    query = f"""
    SELECT pseudo_id, clust AS mcl_id
    FROM pmcl
    WHERE pseudo_id IN ({placeholders})
    """
    cursor.execute(query, pseudo_ids)
    rows = cursor.fetchall()
    if not rows:
        raise ValueError("No MCL cluster information found for the provided pseudogene IDs.")
    return pd.DataFrame(rows)

def fetch_replicons_for_contig(contig_id: int, conn) -> pd.DataFrame:
    """
    Fetch replicon information for a given contig_ID from the MySQL database and return it as a pandas DataFrame.
    Args:
        contig_id (int): The ID of the contig to fetch replicon information for.
        conn: A MySQL database connection object.
    Returns:
        pd.DataFrame: A DataFrame containing the replicon information for the given contig_ID.
    """
    cursor = conn.cursor(dictionary=True)
    query = """
    SELECT id, contig_ID, category, gene_name, start, end
    FROM plasann_kes2021.proteins
    WHERE contig_ID = %s
        AND category = "Replicon"
    """
    cursor.execute(query,(int(contig_id),))
    rows = cursor.fetchall()
    if not rows:
        raise ValueError("No replicon information found for the provided contig_ID.")
    return pd.DataFrame(rows)
