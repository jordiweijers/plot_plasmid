import logging
import re
import threading
from datetime import datetime
from pathlib import Path
from typing import Iterable, List, Optional
import mysql.connector

from plot_plasmid.cache import is_fresh, read_ids, write_ids
from plot_plasmid.config import CategoryRule, DatabaseSchema
from plot_plasmid.models import Category, Contig, Feature

logger = logging.getLogger("main_logger")


class ContigLoader:
    """
    Loads categories, contigs and their features from a MySQL database, using the table and column names in its DatabaseSchema.
    """
    def __init__(self, database: str, schema: DatabaseSchema, host: str, user: str):
        self.database = database
        self.schema = schema
        self.host = host
        self.user = user
        self.conn = mysql.connector.connect(host=host, user=user, database=database)

    def __enter__(self) -> "ContigLoader":
        return self

    def __exit__(self, *exc) -> None:
        self.close()

    def close(self) -> None:
        self.conn.close()

    def _kill_connection(self) -> None:
        """
        Stop the running query on the MySQL server by killing this loader's connection from a second connection.
        Stopping Python alone leaves the query running on the server.
        """
        killer = mysql.connector.connect(host=self.host, user=self.user)
        killer.cursor().execute(f"KILL {self.conn.connection_id}")
        killer.close()
        logger.info(f"Stopped the running query on MySQL (connection {self.conn.connection_id})")

    def _query(self, query: str, params: tuple = ()) -> List[dict]:
        """
        Run a query and return its rows. The query runs in a separate thread, so Ctrl+C reaches Python while MySQL is
        still busy; the query is then stopped on the server too.
        Args:
            query (str): The SQL query, with %s for each parameter.
            params (tuple): The parameters of the query.
        Returns:
            List[dict]: The rows, as dictionaries from column name to value.
        Raises:
            KeyboardInterrupt: If Ctrl+C is pressed while the query runs.
        """
        result = {}
        def run():
            try:
                cursor = self.conn.cursor(dictionary=True)
                cursor.execute(query, params)
                result["rows"] = cursor.fetchall()
                cursor.close()
            except Exception as e:
                result["error"] = e
        thread = threading.Thread(target=run, daemon=True)
        thread.start()
        try:
            while thread.is_alive():
                thread.join(0.1)
        except KeyboardInterrupt:
            self._kill_connection()
            raise
        if "error" in result:
            raise result["error"]
        return result["rows"]

    def _tables_changed(self, query: Optional[Path]) -> Optional[datetime]:
        """
        Find the latest change to the tables a query uses. A table counts as used when its name appears as a word in the
        query text, so a word that only looks like a table name makes the check more careful, never less.
        Args:
            query (Optional[Path]): The path to the SQL file, or None if there is no query.
        Returns:
            Optional[datetime]: The latest create_time or update_time of those tables, or None if there is no query or MySQL knows neither.
        """
        if query is None:
            return None
        text = query.read_text()
        tables = self._query(
            "SELECT table_name AS name, create_time, update_time FROM information_schema.tables WHERE table_schema = %s",
            (self.database,),
        )
        changes = [
            time
            for table in tables if re.search(rf"\b{re.escape(table['name'])}\b", text, re.IGNORECASE)
            for time in (table["create_time"], table["update_time"]) if time is not None
        ]
        return max(changes, default=None)

    def load_categories(self, rules: Iterable[CategoryRule], refresh_cache: bool = False) -> List[Category]:
        """
        Collect the IDs of the categories for this database.
        Args:
            rules (Iterable[CategoryRule]): The category rules from the config file. Rules for other databases are skipped.
            refresh_cache (bool): Whether to rerun the queries even if their saved lists are up to date.
        Returns:
            List[Category]: The categories with their IDs, in config order.
        """
        return [self._load_category(rule, refresh_cache) for rule in rules if rule.database == self.database]

    def _load_category(self, rule: CategoryRule, refresh_cache: bool = False) -> Category:
        """
        Collect the IDs of a category from its path, its query, or its query saved at its path.
        The saved list is used if it is fresh: see cache.is_fresh.
        Args:
            rule (CategoryRule): The category rule from the config file.
            refresh_cache (bool): Whether to rerun the query even if its saved list is up to date.
        Returns:
            Category: The category with its IDs.
        Raises:
            RuntimeError: If MySQL fails to run the query.
            ValueError: If the query does not return exactly one column.
        """
        fresh, reason = is_fresh(rule.path, rule.query, self._tables_changed(rule.query), refresh_cache)

        if fresh:
            ids = read_ids(rule.path)
            logger.info(f"Read {len(ids)} IDs for category '{rule.name}' from {rule.path}, because {reason}")
        else:
            logger.info(f"Running {rule.query.name} for category '{rule.name}', because {reason}. This can take a long time.")
            # The server's clock, so it compares correctly with the tables' create_time and update_time.
            time = self._query("SELECT NOW() AS now")[0]["now"]
            try:
                rows = self._query(rule.query.read_text())
            except mysql.connector.Error as e:
                raise RuntimeError(f"The query of category '{rule.name}' ({rule.query}) failed: {e}") from e
            if rows and len(rows[0]) != 1:
                raise ValueError(f"The query of category '{rule.name}' ({rule.query}) must return one column of IDs, but returns {len(rows[0])}.")
            ids = {str(value) for row in rows for value in row.values() if value is not None}
            if rule.path is not None:
                write_ids(rule.path, ids, rule.query, time)
                logger.info(f"Saved {len(ids)} IDs for category '{rule.name}' at {rule.path}")
        return Category(name=rule.name, color=rule.color, type=rule.type, ids=ids)

    def fetch_contig(self, contig_id: str, categories: List[Category]) -> Contig:
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
            features=self.fetch_features(contig_id, categories),
            category=next((category.name for category in categories if category.type == "contig" and str(contig_id) in category.ids), None),
        )

    def fetch_features(self, contig_id: str, categories: List[Category]) -> List[Feature]:
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
        for feature in features:
            feature.category = next((category.name for category in categories if _feature_matches(category, feature)), None)
        features.sort(key=lambda feature: feature.start)
        return features


def _feature_matches(category: Category, feature: Feature) -> bool:
    """
    Check whether a feature belongs to a category.
    Args:
        category (Category): The category.
        feature (Feature): The feature.
    Returns:
        bool: True if the category is a family category and contains the feature's protein family, or has the feature's
            type (protein or pseudo) and contains its ID. Contig categories never match a feature.
    """
    if category.type == "family":
        return feature.protein_family in category.ids
    return category.type == feature.type and feature.id in category.ids
