import mysql.connector
from typing import List

from plot_plasmid.models import Contig, Feature
from plot_plasmid.config import DatabaseSchema


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
            for row in rows:
                features.append(Feature(
                    id=str(row["id"]),
                    type="pseudo" if row["is_pseudo"] else "protein",
                    start=int(row["start"]),
                    end=int(row["end"]),
                    strand=row["strand"],
                    protein_family=str(row["protein_family"]) if row["protein_family"] is not None else None,
                ))
        features.sort(key=lambda feature: feature.start)
        return features
