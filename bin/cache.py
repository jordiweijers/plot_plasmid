import hashlib
import logging
import os
from datetime import datetime, time
from pathlib import Path
from typing import Iterable, Optional, Set, Tuple, Dict

import mysql.connector

from plot_plasmid.database import Database

DATE_FORMAT = "%Y-%m-%d %H:%M:%S"

logger = logging.getLogger(__name__)


def load_ids(database: Database, query: str, cache: Optional[Path] = None, refresh: bool = False) -> Set[str]:
    """
    Get the IDs a query returns, from the cache if they were saved before and are still up to date, otherwise by running the query.
    Args:
        database (Database): The database to run the query on.
        query (str): The SQL query, which must return one column of IDs.
        cache (Optional[Path]): The Path to save the IDs in, or None to always run the query and save nothing.
        refresh (bool): Whether to run the query even if its saved IDs are up to date.
    Returns:
        Set[str]: The IDs, as strings.
    """
    if cache is None:
        start = time.monotonic()
        ids = _query_ids(database, query)
        finish = time.monotonic() - start
        if finish > 60:
            logger.warning(f"The query took {finish:.2f} seconds. Give a cache path to save its IDs, so it only needs to run once.")
        return ids

    fresh, reason = _is_fresh(cache, database, query, refresh)
    if fresh:
        ids = read_ids(cache)
        logger.info(f"Read {len(ids)} IDs from {cache}, because {reason}")
        return ids

    logger.info(f"Running the query for {cache}, because {reason}. This can take a long time.")
    created = database.now()
    ids = _query_ids(database, query)
    _write_ids(cache, ids, database, query, created)
    logger.info(f"Saved {len(ids)} IDs at {cache}")
    return ids


def read_ids(path: Path) -> Set[str]:
    """
    Read a list of IDs from a file.
    Args:
        path (Path): The path to a file with one ID per line. Blank lines and lines starting with # are skipped.
    Returns:
        Set[str]: The IDs in the file.
    """
    with open(path) as fh:
        return {line.strip() for line in fh if line.strip() and not line.startswith("#")}


def _query_ids(database: Database, query: str) -> Set[str]:
    """
    Run a query that returns one column of IDs.
    Args:
        database (Database): The database to run the query on.
        query (str): The SQL query.
    Returns:
        Set[str]: The IDs, as strings, so they compare equal to the IDs in the pipeline. NULL values are skipped.
    Raises:
        ValueError: If the query uses @contig_id or does not return exactly one column.
        RuntimeError: If MySQL fails to run the query.
    """
    if "@contig_id" in query:
        raise ValueError("The query uses @contig_id, which is not allowed. Queries run on the whole database, not per contig")
    columns = database.columns(query)
    if len(columns) != 1:
        raise ValueError(f"The query must return exactly one column, but returns {len(columns)}: {', '.join(columns)}.")
    try:
        rows = database.query(query)
    except mysql.connector.Error as e:
        raise RuntimeError(f"The query failed on database '{database.name}': {e}") from e
    return {str(value) for row in rows for value in row.values() if value is not None}


def _is_fresh(cache: Path, database: Database, query: str, refresh: bool) -> Tuple[bool, str]:
    """
    Check whether the saved IDs of a query can be used instead of running it.
    Args:
        cache (Path): The file with the saved IDs.
        database (Database): The database the query runs on.
        query (str): The SQL query.
        refresh (bool): Whether to run the query even if its saved IDs are up to date.
    Returns:
        Tuple[bool, str]: Whether the saved IDs can be used, and the reason, for logging.
    """
    if refresh:
        return False, "a refresh was asked"
    if not cache.is_file():
        return False, "it is not saved yet"
    header = _read_header(cache)
    if not {"database", "created", "query_sha256", "ids_sha256"} <= header.keys():
        return False, "its header is incomplete"
    _check_cache(cache, header["ids_sha256"])
    if header["database"] != database.name:
        return False, f"it was saved for database '{header['database']}', not '{database.name}'"
    if header["query_sha256"] != _sha256(query):
        return False, "the query has changed since the cache was created"
    created = datetime.strptime(header["created"], DATE_FORMAT)
    changed = database.tables_changed(query)
    if changed is not None and changed > created:
        return False, f"a table it uses changed on {changed:{DATE_FORMAT}}, after it was saved on {created:{DATE_FORMAT}}"
    return True, f"it is up to date (saved on {created:{DATE_FORMAT}})"


def _sha256(text: str) -> str:
    return hashlib.sha256(text.encode()).hexdigest()


def _ids_sha256(ids: Iterable[str]) -> str:
    """The hash of a set of IDs, as saved in and checked against the cache header."""
    return _sha256('\n'.join(sorted(ids)))


def _read_header(cache: Path) -> Dict[str, str]:
    """
    Read the header of a cache file into a dictionary.
    Args:
        cache (Path): The file with the saved IDs.
    Returns:
        Dict[str, str]: The header key-value pairs.
    """
    header = {}
    with open(cache) as fh:
        for line in fh:
            if line.startswith("#   "):
                continue
            if not line.startswith("#"):
                break
            if ":" in line:
                key, _, value = line[1:].partition(":")
                header[key.strip()] = value.strip()
    return header


def _check_cache(cache: Path, ids_sha256: str) -> None:
    """
    Check that the IDs in a cache file match the expected hash.
    Args:
        cache (Path): The file with the saved IDs.
        ids_sha256 (str): The expected hash of the IDs.
    Raises:
        ValueError: If the IDs do not match the expected hash.
    """
    if _ids_sha256(read_ids(cache)) != ids_sha256:
        raise ValueError(f"The cache file in {cache} was changed unexpectedly.\nTo use your own IDs, make the rule with CategoryRule.from_file() or CategoryRule(category, ids).\nTo overwrite the file use refresh=True.")


def _write_ids(cache: Path, ids: Iterable[str], database: Database, query: str, created: datetime) -> None:
    """
    Save the IDs of a query, with a header recording the database, the creation date and the query. The file is written
    under a temporary name and then renamed, so an interrupted write never leaves a partial list that looks complete.
    Args:
        cache (Path): The file to save the IDs in.
        ids (Iterable[str]): The IDs to save.
        database (Database): The database the query ran on.
        query (str): The SQL query that returned the IDs.
        created (datetime): When the query started, so table changes made while it ran make the IDs out of date.
    """
    cache.parent.mkdir(parents=True, exist_ok=True)
    tmp_path = cache.with_name(cache.name + ".tmp")
    with open(tmp_path, "w") as fh:
        fh.write(f"# database: {database.name}\n")
        fh.write(f"# created: {created:{DATE_FORMAT}}\n")
        fh.write(f"# query_sha256: {_sha256(query)}\n")
        fh.write(f"# ids_sha256: {_ids_sha256(ids)}\n")
        fh.write("# query:\n")
        for line in query.strip().splitlines():
            fh.write(f"#   {line}\n")
        for id in sorted(ids):
            fh.write(f"{id}\n")
    os.replace(tmp_path, cache)
