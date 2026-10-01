import hashlib
import os
from datetime import datetime
from pathlib import Path
from typing import Iterable, Optional, Set, Tuple

DATE_FORMAT = "%Y-%m-%d %H:%M:%S"


def query_hash(query: Path) -> str:
    """
    Compute the hash of an SQL file, so a saved list can be checked against the query that made it.
    Args:
        query (Path): The path to the SQL file.
    Returns:
        str: The SHA-256 hash of the file contents.
    """
    return hashlib.sha256(query.read_bytes()).hexdigest()



def read_header(path: Path) -> Optional[Tuple[str, datetime]]:
    """
    Read the query hash and creation date from the header of a saved list, as written by write_ids.
    Args:
        path (Path): The path to the saved list.
    Returns:
        Optional[Tuple[str, datetime]]: The query hash and the creation date, or None if the file does not exist or its header lacks either.
    """
    if not path.is_file():
        return None
    header = {}
    with open(path) as fh:
        for line in fh:
            if not line.startswith("#"):
                break
            key, _, value = line[1:].partition(":")
            header[key.strip()] = value.strip()
    try:
        return header["query_sha256"], datetime.strptime(header["created"], DATE_FORMAT)
    except (KeyError, ValueError):
        return None


def read_ids(path: Path) -> Set[str]:
    """
    Read a list of IDs from a file.
    Args:
        path (Path): The path to a file with one ID per line. Blank lines and lines starting with # are skipped.
    Returns:
        Set[str]: The IDs in the file.
    """
    with open(path) as f:
        return {line.strip() for line in f if line.strip() and not line.startswith("#")}


def write_ids(path: Path, ids: Iterable[str], query: Path, created: datetime) -> None:
    """
    Save a list of IDs with a header recording the query that made it. The file is written under a temporary name and then
    renamed, so an interrupted write never leaves a partial list that looks complete.
    Args:
        path (Path): The path to save the list at. Missing directories are created.
        ids (Iterable[str]): The IDs to save.
        query (Path): The path to the SQL file that returned the IDs.
        created (datetime): When the query started, so table changes made while it ran make the list out of date.
    """
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp_path = path.with_name(path.name + ".tmp")
    with open(tmp_path, "w") as f:
        f.write(f"# query: {query}\n")
        f.write(f"# query_sha256: {query_hash(query)}\n")
        f.write(f"# created: {created.strftime(DATE_FORMAT)}\n")
        for id in sorted(ids):
            f.write(f"{id}\n")
    os.replace(tmp_path, path)

def is_fresh(
        path: Optional[Path],
        query: Optional[Path],
        tables_changed: Optional[datetime],
        refresh_cache: bool = False,
    ) -> Tuple[bool, str]:
    """
    Check whether a category's list can be read from its path instead of running its query.
    Args:
        path (Optional[Path]): The path to the list, or None if the category has none.
        query (Optional[Path]): The path to the SQL file that makes the list, or None if the category has none.
        tables_changed (Optional[datetime]): The latest change to the tables the query uses, or None if unknown.
        refresh_cache (bool): Whether to rerun the query even if its saved list is up to date.
    Returns:
        Tuple[bool, str]: Whether the list can be read from path, and the reason, for logging.
    """
    if query is None:
        return True, "it has no query"
    if path is None:
        return False, "it has no path to save its IDs"
    if refresh_cache:
        return False, "--refresh-cache was given"
    if not path.is_file():
        return False, "it does not exist yet"
    header = read_header(path)
    if header is None:
        return False, "it has no header with the query hash and creation date"
    saved_hash, created = header
    if saved_hash != query_hash(query):
        return False, f"{query.name} changed since the list was created"
    if tables_changed is not None and tables_changed > created:
        return False, f"a table in {query.name} changed on {tables_changed:{DATE_FORMAT}}, after the list was created on {created:{DATE_FORMAT}}"
    return True, f"it is up to date (created {created:{DATE_FORMAT}})"
