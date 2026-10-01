from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Dict, Literal, Optional, Tuple

import yaml
from matplotlib.colors import is_color_like

CATEGORY_TYPES = ("protein", "pseudo", "family", "contig")


@dataclass(frozen=True)
class FeatureSource:
    """
    Describes how one table of features is structured in the database.
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
    Describes the schema of the database, including contigs and features tables.
    Attributes:
        contigs (ContigSource): The source of the contigs table.
        features (Tuple[FeatureSource, ...]): The sources of the features tables.
    """
    contigs: ContigSource
    features: Tuple[FeatureSource, ...]


@dataclass(frozen=True)
class CategoryRule:
    """
    Assigns a category to the proteins, pseudogenes, protein families or contigs whose IDs are listed in a file or returned by an SQL query.
    Attributes:
        name (str): The category label shown in the plot legend.
        color (str): The color of the features or contig labels in this category.
        type (str): What the IDs refer to: "protein", "pseudo", "family" or "contig".
        database (str): The database this category applies to.
        query (Optional[Path]): The path to an SQL file that returns one column of IDs, if any.
        path (Optional[Path]): The path to a file with one ID per line, if any. With a query, the query's IDs are saved here and reused.
    """
    name: str
    color: str
    type: Literal["protein", "pseudo", "family", "contig"]
    database: str
    query: Optional[Path] = None
    path: Optional[Path] = None


@dataclass(frozen=True)
class Config:
    """
    The complete plot_plasmid configuration, as read from a YAML file by load_config.
    Attributes:
        host (str): The MySQL host.
        user (str): The MySQL user.
        databases (Dict[str, DatabaseSchema]): The schema of each database, by database name.
        categories (Tuple[CategoryRule, ...]): The category rules, in priority order (first match wins).
    """
    host: str
    user: str
    databases: Dict[str, DatabaseSchema]
    categories: Tuple[CategoryRule, ...]


def _parse_schema(database: str, entry: dict) -> DatabaseSchema:
    """
    Build a DatabaseSchema from one entry under 'databases' in the config file.
    Args:
        database (str): The name of the database, used in error messages.
        entry (dict): The parsed YAML entry of the database.
    Returns:
        DatabaseSchema: The schema of the database.
    Raises:
        ValueError: If a key is missing or unknown, or is_pseudo is not a string.
    """
    try:
        contigs = ContigSource(**entry["contigs"])
        features = tuple(FeatureSource(**source) for source in entry["features"])
    except (KeyError, TypeError) as e:
        raise ValueError(f"Invalid schema for database '{database}': {e}") from e
    for source in features:
        if not isinstance(source.is_pseudo, str):
            raise ValueError(
                f"is_pseudo of table '{source.table}' in database '{database}' must be a string with an SQL expression, "
                f"got {source.is_pseudo!r}. Put quotes around it in the config file, e.g. \"FALSE\"."
            )
    return DatabaseSchema(contigs=contigs, features=features)


def _parse_category(entry: dict, config_dir: Path, databases: Dict[str, DatabaseSchema]) -> CategoryRule:
    """
    Build a CategoryRule from one entry under 'categories' in the config file.
    Args:
        entry (dict): The parsed YAML entry of the category.
        config_dir (Path): The directory of the config file, which the query and path are relative to.
        databases (Dict[str, DatabaseSchema]): The databases in the config file, to check the category's database against.
    Returns:
        CategoryRule: The category rule, with its query and path resolved relative to the config file.
    Raises:
        ValueError: If a key is missing or unknown, the type, database or color is invalid, neither query nor path is given,
            a file does not exist, or the query uses @contig_id.
    """
    try:
        name, color, type_, database = entry["name"], entry["color"], entry["type"], entry["database"]
    except KeyError as e:
        raise ValueError(f"Category {entry.get('name', entry)!r} is missing the key {e}") from e
    unknown = set(entry) - {"name", "color", "type", "database", "query", "path"}
    if unknown:
        raise ValueError(f"Category '{name}' has unknown keys: {', '.join(sorted(unknown))}")
    if type_ not in CATEGORY_TYPES:
        raise ValueError(f"Category '{name}' has type '{type_}', which must be one of: {', '.join(CATEGORY_TYPES)}.")
    if database not in databases:
        raise ValueError(f"Category '{name}' has database '{database}', which is not under 'databases' in the config file.")
    if not is_color_like(color):
        raise ValueError(f"Category '{name}' has an invalid color {color!r}. Put quotes around colors in the config file, e.g. '#5e2bff'.")

    query = config_dir / entry["query"] if entry.get("query") else None
    path = config_dir / entry["path"] if entry.get("path") else None
    if query is None and path is None:
        raise ValueError(f"Category '{name}' needs a query, a path, or both.")
    if query is not None:
        if not query.is_file():
            raise ValueError(f"The query of category '{name}' does not exist: {query}")
        if "@contig_id" in query.read_text():
            raise ValueError(
                f"The query of category '{name}' uses @contig_id, but queries run on the whole database, not per contig. "
                f"Remove @contig_id from {query}."
            )
    elif not path.is_file():
        raise ValueError(f"The path of category '{name}' does not exist: {path}")
    return CategoryRule(name=name, color=color, type=type_, database=database, query=query, path=path)


def load_config(path: Path) -> Config:
    """
    Load and validate the plot_plasmid configuration from a YAML file.
    Args:
        path (Path): The path to the YAML config file.
    Returns:
        Config: The parsed configuration.
    Raises:
        ValueError: If the config file is missing a section or contains an invalid schema or category.
    """
    with open(path) as fh:
        raw = yaml.safe_load(fh)
    try:
        connection, databases, categories = raw["connection"], raw["databases"], raw["categories"]
    except KeyError as e:
        raise ValueError(f"Config file {path} is missing the section {e}") from e
    schemas = {database: _parse_schema(database, entry) for database, entry in databases.items()}
    return Config(
        host=connection["host"],
        user=connection["user"],
        databases=schemas,
        categories=tuple(_parse_category(entry, path.parent, schemas) for entry in categories),
    )
