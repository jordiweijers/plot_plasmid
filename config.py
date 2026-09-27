from __future__ import annotations

import re
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, Optional, Tuple, Union

import yaml
from matplotlib.colors import is_color_like


@dataclass(frozen=True)
class FeatureSource:
    """
    Describes how one table of features is structured in the database.
    Attributes:
        table (str): The name of the table containing the features.
        id_column (str): The name of the column containing feature IDs.
        contig_column (str): The name of the column containing contig IDs.
        name_column (Optional[str]): The name of the column containing feature names.
        symbol_column (Optional[str]): The name of the column containing feature symbols, or None if the table has none.
        is_pseudo (str): The SQL boolean expression to determine if a feature is a pseudogene.
        mcl_table (str): The name of the MCL table.
        mcl_id_column (str): The name of the column containing MCL IDs.
        mcl_family_column (str): The name of the column containing MCL family IDs.
        amr_table (Optional[str]): The name of the AMR table, if any.
        amr_id_column (Optional[str]): The name of the column containing AMR IDs, if any.
        amr_aro_column (Optional[str]): The name of the column containing AMR ARO IDs, if any.
    """
    table: str
    id_column: str
    contig_column: str
    is_pseudo: str
    mcl_table: str
    mcl_id_column: str
    mcl_family_column: str
    name_column: Optional[str] = None
    symbol_column: Optional[str] = None
    amr_table: Optional[str] = None
    amr_id_column: Optional[str] = None
    amr_aro_column: Optional[str] = None


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
        replicon_mcl_file (Optional[Path]): The path to the replicon MCL file, if any.
        conjugation_mcl_file (Optional[Path]): The path to the conjugation MCL file, if any.
    """
    contigs: ContigSource
    features: Tuple[FeatureSource, ...]
    replicon_mcl_file: Optional[Path] = None
    conjugation_mcl_file: Optional[Path] = None


@dataclass(frozen=True)
class CategoryRule:
    """
    Assigns a category to protein families whose main name or symbol matches a pattern.
    Attributes:
        name (str): The category label shown in the plot legend.
        column (str): Which annotation to match against, either "name" or "symbol".
        pattern (re.Pattern): The compiled Python regular expression, searched anywhere in the value.
        color (str): The color of the features in this category.
    """
    name: str
    column: str
    pattern: re.Pattern
    color: str


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


def _parse_category(entry: dict) -> CategoryRule:
    """
    Build a CategoryRule from one entry under 'categories' in the config file.
    Args:
        entry (dict): The parsed YAML entry of the category.
    Returns:
        CategoryRule: The category rule, with its pattern compiled.
    Raises:
        ValueError: If a key is missing or unknown, the column is not "name" or "symbol", the pattern is not a valid regular expression, or the color is not a valid color.
    """
    try:
        name, column, pattern, color = entry["name"], entry["column"], entry["pattern"], entry["color"]
    except KeyError as e:
        raise ValueError(f"Category {entry.get('name', entry)!r} is missing the key {e}") from e
    unknown = set(entry) - {"name", "column", "pattern", "color"}
    if unknown:
        raise ValueError(f"Category '{name}' has unknown keys: {', '.join(sorted(unknown))}")
    if column not in ("name", "symbol"):
        raise ValueError(f"Category '{name}' has column '{column}', which must be 'name' or 'symbol'.")
    if any(ord(char) < 32 for char in pattern):
        raise ValueError(
            f"The pattern of category '{name}' contains a control character, probably from a backslash inside double quotes. "
            f"Use single quotes around patterns in the config file, e.g. '\\bIS3\\b'."
        )
    try:
        compiled = re.compile(pattern)
    except re.error as e:
        raise ValueError(f"The pattern of category '{name}' is not a valid regular expression: {e}") from e
    if not is_color_like(color):
        raise ValueError(f"Category '{name}' has an invalid color {color!r}. Put quotes around colors in the config file, e.g. '#5e2bff'.")
    return CategoryRule(name=name, column=column, pattern=compiled, color=color)


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
        connection, databases = raw["connection"], raw["databases"]
    except KeyError as e:
        raise ValueError(f"Config file {path} is missing the section {e}") from e
    return Config(
        host=connection["host"],
        user=connection["user"],
        databases={database: _parse_schema(database, entry) for database, entry in databases.items()},
        categories=tuple(_parse_category(entry) for entry in raw.get("categories") or []),
    )
