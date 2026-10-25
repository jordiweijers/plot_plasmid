from typing import Literal, List

from contigviz.database import Database
from contigviz.models import Contig, ContigID, Feature, FamilyID, ProteinID, PseudoID


class ContigLoader:
    def __init__(self, database: Database):
        self.database = database

    def fetch_contig(self, contig_id: ContigID, topology: Literal["linear", "circular"] = "circular") -> Contig:
        """
        Fetch a contig and its features from the database.
        Args:
            contig_id (ContigID): The ID of the contig to fetch.
            topology (Literal["linear", "circular"]): The topology of the contig (default = "circular").
        Returns:
            Contig: The contig object with its features populated.
        Raises:
            ValueError: If the contig is not found in the database.
        """
        source = self.database.schema.contigs
        rows = self.database.query(
            f"SELECT {source.length_column} AS length FROM {source.table} WHERE {source.id_column} = %s",
            (contig_id.source_id,)
        )
        if not rows:
            raise ValueError(f"Contig '{contig_id}' not found in database '{self.database.name}'.")
        length = int(rows[0]["length"])
        if length < 1:
            raise ValueError(f"Contig '{contig_id}' has length {length} in database '{self.database.name}'.")
        features = self.fetch_features(contig_id)
        if topology == "circular":
            # Some databases, e.g. ESKAPEE25, store a feature that wraps around the origin as running past the end of the contig. Counting on from 1 after the end turns it into start > end; positions inside the contig stay.
            for feature in features:
                feature.start = (feature.start - 1) % length + 1
                feature.end = (feature.end - 1) % length + 1
        return Contig(
            id=contig_id,
            start=1,
            end=length,
            topology=topology,
            features=features,
        )

    def fetch_features(self, contig_id: ContigID) -> List[Feature]:
        """
        Fetch features for a given contig from the database.
        Args:
            contig_id (ContigID): The ID of the contig to fetch features for.
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
                (contig_id.source_id,),
            )
            for row in rows:
                feature_class = PseudoID if row["is_pseudo"] else ProteinID
                features.append(Feature(
                    id=feature_class(source_id=str(row["id"])),
                    start=int(row["start"]),
                    end=int(row["end"]),
                    strand=row["strand"],
                    protein_family=FamilyID(source_id=str(row["protein_family"])) if row["protein_family"] is not None else None,
                ))
        return features
