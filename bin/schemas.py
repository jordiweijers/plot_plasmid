"""
The schemas of the databases of our group, ready to use: from plot_plasmid.schemas import kes2021.

To add a database, copy one of the schemas below and change the names; SCHEMAS picks it up by itself. To connect to a
schema as another user or on another host, copy it with dataclasses.replace(kes2021, host=..., user=...).
"""
from dataclasses import dataclass
from typing import List


@dataclass(frozen=True)
class FeatureSource:
    """
    Describes how one table of features is structured in the database. The table and column names, and is_pseudo, are
    inserted into the SQL queries as-is, so is_pseudo is evaluated by MySQL, e.g. "TRUE" for a table of pseudogenes.
    Attributes:
        table (str): The name of the table containing the features.
        id_column (str): The name of the column containing feature IDs.
        contig_column (str): The name of the column containing contig IDs.
        is_pseudo (str): The SQL boolean expression to determine if a feature is a pseudogene.
        mcl_table (str): The name of the MCL table.
        mcl_id_column (str): The name of the column containing feature IDs in the MCL table.
        mcl_family_column (str): The name of the column containing MCL family IDs.
    """
    table: str
    id_column: str
    contig_column: str
    is_pseudo: str
    mcl_table: str
    mcl_id_column: str
    mcl_family_column: str


@dataclass(frozen=True)
class ContigSource:
    """
    Describes how one table of contigs is structured in the database.
    Attributes:
        table (str): The name of the table containing the contigs.
        id_column (str): The name of the column containing contig IDs.
        length_column (str): The name of the column containing contig lengths.
    """
    table: str
    id_column: str
    length_column: str


@dataclass(frozen=True)
class DatabaseSchema:
    """
    Describes where a database is and its schema, including contigs and features tables.
    Attributes:
        database (str): The name of the database.
        host (str): The MySQL host the database is on.
        user (str): The MySQL user to connect as.
        contigs (ContigSource): The source of the contigs table.
        features (List[FeatureSource]): The sources of the features tables.
    """
    database: str
    host: str
    user: str
    contigs: ContigSource
    features: List[FeatureSource]


HOST = "localhost"
USER = "gmg"


kes2021 = DatabaseSchema(
    database="kes2021",
    host=HOST,
    user=USER,
    contigs=ContigSource(table="contigs", id_column="contig_ID", length_column="Size"),
    features=[
        FeatureSource(
            table="proteins", id_column="id", contig_column="contig_ID", is_pseudo="FALSE",
            mcl_table="mcl30", mcl_id_column="id", mcl_family_column="clust",
        ),
        FeatureSource(
            table="pseudogenes", id_column="pseudo_id", contig_column="contig_ID", is_pseudo="TRUE",
            mcl_table="pmcl", mcl_id_column="pseudo_id", mcl_family_column="clust",
        ),
    ],
)

ESKAPEE25 = DatabaseSchema(
    database="ESKAPEE25",
    host=HOST,
    user=USER,
    contigs=ContigSource(table="contigs", id_column="contig_id", length_column="size"),
    features=[
        FeatureSource(
            table="proteins", id_column="id", contig_column="contig_id", is_pseudo="FALSE",
            mcl_table="mcl", mcl_id_column="id", mcl_family_column="clust",
        ),
        FeatureSource(
            table="pseudogenes", id_column="pseudo_id", contig_column="contig_id", is_pseudo="TRUE",
            mcl_table="pmcl", mcl_id_column="pseudo_id", mcl_family_column="clust",
        ),
    ],
)

acinetobacter_2023 = DatabaseSchema(
    database="acinetobacter_2023",
    host=HOST,
    user=USER,
    contigs=ContigSource(table="contigs", id_column="contig_ID", length_column="Size"),
    features=[
        FeatureSource(
            table="proteins", id_column="id", contig_column="contig_ID", is_pseudo="FALSE",
            mcl_table="mcl30", mcl_id_column="id", mcl_family_column="clust",
        ),
    ],
)

pantoea2025_v2 = DatabaseSchema(
    database="pantoea2025_v2",
    host=HOST,
    user=USER,
    contigs=ContigSource(table="contigs", id_column="contig_id", length_column="contig_length"),
    features=[
        FeatureSource(
            table="features_cds", id_column="cds_id", contig_column="contig_id",
            is_pseudo="FIND_IN_SET('pseudo', REPLACE(attributes, ';', ',')) > 0",
            mcl_table="mcls", mcl_id_column="cds_id", mcl_family_column="mcl_id",
        ),
    ],
)


# Every schema defined above, by database name, so a new schema only needs to be written once.
SCHEMAS = {value.database: value for value in list(globals().values()) if isinstance(value, DatabaseSchema)}
