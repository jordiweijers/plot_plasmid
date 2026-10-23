from typing import Literal, List

from plot_plasmid.database import Database
from plot_plasmid.models import Contig, Feature, FamilyID


class ContigLoader:
    def __init__(self, database: Database):
        self.database = database

    def fetch_contig(self, contig_id: str, topology: Literal["linear", "circular"] = "circular") -> Contig:
        """
        Fetch a contig and its features from the database.
        Args:
            contig_id (str): The ID of the contig to fetch.
            topology (Literal["linear", "circular"]): The topology of the contig. Defaults to "circular".
        Returns:
            Contig: The contig object with its features populated.
        Raises:
            ValueError: If the contig is not found in the database.
        """
        source = self.database.schema.contigs
        rows = self.database.query(
            f"SELECT {source.length_column} AS length FROM {source.table} WHERE {source.id_column} = %s",
            (contig_id,)
        )
        if not rows:
            raise ValueError(f"Contig '{contig_id}' not found in database '{self.database.name}'.")
        return Contig(
            source_id=str(contig_id),
            start=1,
            end=int(rows[0]["length"]),
            topology=topology,
            features=self.fetch_features(contig_id),
        )

    def fetch_features(self, contig_id: str) -> List[Feature]:
        """
        Fetch features for a given contig from the database.
        Args:
            contig_id (str): The ID of the contig to fetch features for.
        Returns:
            features (List[Feature]): A list of features associated with the contig.
        """
        features = []
        for source in self.database.schema.features:
            rows = self.database.query(
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
                    source_id=str(row["id"]),
                    type="pseudo" if row["is_pseudo"] else "protein",
                    start=int(row["start"]),
                    end=int(row["end"]),
                    strand=row["strand"],
                    protein_family=FamilyID(source_id=str(row["protein_family"])) if row["protein_family"] is not None else None,
                ))
        return features
