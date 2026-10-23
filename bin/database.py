import logging
import re
import threading
from dataclasses import dataclass
from datetime import datetime
from typing import List, Optional, Tuple
import mysql.connector
from mysql.connector.constants import ClientFlag

logger = logging.getLogger(__name__)


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
        database (str): The name of the database.
        contigs (ContigSource): The source of the contigs table.
        features (List[FeatureSource]): The sources of the features tables.
    """
    database: str
    contigs: ContigSource
    features: List[FeatureSource]


class Database:
    """
    Represents a connection to a MySQL database and provides methods to read from it. It only runs single SELECT
    queries, so it cannot change the database, even with a user that may.
    """
    def __init__(self, schema: DatabaseSchema, host: str, user: str):
        self.schema = schema
        self.name = schema.database
        self.host = host
        self.user = user
        # Without this flag the server runs every statement in a query, so "SELECT ...; DROP ..." would drop.
        self.conn = mysql.connector.connect(host=host, user=user, database=self.name, client_flags=[-ClientFlag.MULTI_STATEMENTS])

    def __enter__(self) -> "Database":
        return self

    def __exit__(self, *exc) -> None:
        self.close()

    def close(self) -> None:
        self.conn.close()

    def _kill_connection(self) -> None:
        """
        Stop the running query on the MySQL server by killing this database's connection from a second connection.
        Stopping Python alone leaves the query running on the server.
        """
        killer = mysql.connector.connect(host=self.host, user=self.user)
        killer.cursor().execute(f"KILL {self.conn.connection_id}")
        killer.close()
        logger.info(f"Stopped the running query on MySQL (connection {self.conn.connection_id})")

    def columns(self,  query: str, params: Tuple = ()) -> List[str]:
        """
        Get the column names of the result of a query. And additionally check that the query is a single SELECT statement and correct syntax
        Args:
            query (str): The SQL query, with %s for each parameter.
            params (tuple): The parameters of the query.
        Returns:
            List[str]: The column names of the result.
        Raises:
            ValueError: If the query is not a single SELECT statement or has incorrect syntax.
        """
        checked = f"SELECT * FROM (\n{query.strip().rstrip(';')}\n) AS checked LIMIT 0"
        with self.conn.cursor() as cursor:
            try:
                cursor.execute(checked, params)
                columns = [desc[0] for desc in cursor.description]
                cursor.fetchall()
            except mysql.connector.Error as e:
                raise ValueError(f"Query is not a single SELECT statement or has incorrect syntax that database '{self.name}' cannot execute.") from e
        return columns

    def query(self, query: str, params: Tuple = ()) -> List[dict]:
        """
        Run a query and return its rows. The query runs in a separate thread, so Ctrl+C reaches Python while MySQL is still busy; the query is then stopped on the server too.
        Args:
            query (str): The SQL query, with %s for each parameter.
            params (tuple): The parameters of the query.
        Returns:
            List[dict]: The rows, as dictionaries from column name to value.
        Raises:
            KeyboardInterrupt: If Ctrl+C is pressed while the query runs.
        """
        result = {}
        self.columns(query, params)
        def run():
            try:
                with self.conn.cursor(dictionary=True) as cursor:
                    cursor.execute(query, params)
                    result["rows"] = cursor.fetchall()
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

    def tables_changed(self, query: str) -> Optional[datetime]:
        """
        Find the latest change to the tables a query uses. A table counts as used when its name appears as a word in the query text, so a word that only looks like a table name makes the check more careful, never less.
        Args:
            query (str): The SQL query as a string.
        Returns:
            Optional[datetime]: The latest create_time or update_time of those tables, or None if MySQL knows neither.
                MySQL often leaves update_time empty.
        """
        if not query:
            raise ValueError("Query must be a non-empty string.")
        tables = self.query(
            "SELECT table_name AS name, create_time, update_time FROM information_schema.tables WHERE table_schema = %s",
            (self.name,),
        )
        changes = [
            time
            for table in tables if re.search(rf"\b{re.escape(table['name'])}\b", query, re.IGNORECASE)
            for time in (table["create_time"], table["update_time"]) if time is not None
        ]
        return max(changes, default=None)

    def now(self) -> datetime:
        """
        Get the current date and time from the MySQL server, so it compares correctly with the tables' create_time and
        update_time.
        Returns:
            datetime: The current date and time according to the MySQL server.
        """
        result = self.query("SELECT NOW() AS now")
        return result[0]["now"]
