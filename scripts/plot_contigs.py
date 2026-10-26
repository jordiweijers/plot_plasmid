import argparse
import logging
from collections import Counter
from pathlib import Path
from typing import List, Optional

from contigviz import (
    Category, CategoryRule, ContigID, ContigLoader, Database, FamilyID, ProteinID,
    assign_categories, order_by_clustering, plot_contigs,
)
from contigviz.schemas import SCHEMAS

logger = logging.getLogger(__name__)

QUERIES = Path(__file__).resolve().parent.parent / "queries"
CATEGORIES = [
    (Category("AMR", "#33ffff"), ProteinID, "amr_proteins.sql"),
    (Category("Replication", "#5e2bff"), FamilyID, "replication_families.sql"),
    (Category("Conjugation", "#8b0000"), FamilyID, "conjugation_families.sql"),
]


def setup_logging(log_file: Path) -> None:
    """
    Log the messages of this script and contigviz to the console (INFO and up) and to a file (everything). Other packages, such as matplotlib, only log their warnings.
    Args:
        log_file (Path): The file to append the log to.
    """
    console_handler = logging.StreamHandler()
    console_handler.setLevel(logging.INFO)
    console_handler.setFormatter(logging.Formatter("%(message)s"))
    file_handler = logging.FileHandler(log_file, mode="a")
    file_handler.setFormatter(logging.Formatter("[%(asctime)s] %(levelname)-8s %(name)s: %(message)s", datefmt="%Y-%m-%d %H:%M:%S"))
    # The handlers go on the root logger, which the messages of every logger reach, also those of contigviz.
    logging.basicConfig(level=logging.WARNING, handlers=[console_handler, file_handler])
    for name in (__name__, "contigviz"):
        logging.getLogger(name).setLevel(logging.DEBUG)


def build_rules(database: Database, cache_dir: Optional[Path] = None, refresh: bool = False) -> List[CategoryRule]:
    """
    Make a rule for each category in CATEGORIES that has a query file for this database.
    Args:
        database (Database): The database to run the queries on.
        cache_dir (Optional[Path]): The directory to save the IDs of each query in, under the name of the database, or None to always run the queries.
        refresh (bool): Whether to run the queries even if their saved IDs are up to date.
    Returns:
        List[CategoryRule]: The rules, in the priority order of CATEGORIES.
    """
    rules = []
    for category, id_class, query_file in CATEGORIES:
        query_path = QUERIES / database.name / query_file
        if not query_path.exists():
            logger.info(f"Category '{category.label}': no {query_file} for database '{database.name}', so it is skipped")
            continue
        cache = cache_dir / database.name / f"{query_path.stem}.txt" if cache_dir is not None else None
        rules.append(CategoryRule.from_query(database, category, id_class, query_path.read_text(), cache, refresh))
    return rules

def read_contig_file(path: Path) -> List[str]:
    """
    Read contig IDs from a file.
    Args:
        path (Path): A file with one contig ID per line. Blank lines are skipped, and anything after a tab is ignored.
    Returns:
        List[str]: The contig IDs, in file order.
    """
    contig_ids = []
    with open(path) as fh:
        for line in fh:
            contig_id = line.split("\t")[0].strip()
            if contig_id:
                contig_ids.append(contig_id)
    return contig_ids

def parse_arguments():
    parser = argparse.ArgumentParser(description="Plot contigs from one database with pygenomeviz, colored by the AMR, replication and conjugation categories of that database")
    input_group = parser.add_mutually_exclusive_group(required=True)
    input_group.add_argument(
        "-c", "--contig_ids", type=str, nargs="+", help="List of contig IDs to plot (space-separated)"
    )
    input_group.add_argument(
        "-f", "--contig_file", type=Path, help="Path to a file with one contig ID per line"
    )
    parser.add_argument(
        "-d", "--database", type=str, required=True, choices=sorted(SCHEMAS), help="Name of the database to load the contigs from"
    )
    parser.add_argument(
        "-o", "--output_file", type=Path, required=True,
        help="Output file path for the plot; the extension sets the format (e.g. contigs.svg or contigs.png). A log file <name>.log is written next to it"
    )
    parser.add_argument(
        "-s", "--order-by-similarity", action="store_true",
        help="Reorder contigs by similarity of their protein families instead of the input order"
    )
    cache_group = parser.add_mutually_exclusive_group()
    cache_group.add_argument(
        "--cache-dir", type=Path, default=Path("cache"),
        help="Directory to save the IDs of the category queries in, so they only run once (default: cache/ in the current directory)"
    )
    cache_group.add_argument(
        "--no-cache", action="store_true", help="Do not save the IDs of the category queries, but run them on every plot"
    )
    parser.add_argument(
        "--refresh-cache", action="store_true",
        help="Rerun the category queries and overwrite their saved IDs, even if they are up to date"
    )
    return parser.parse_args()


def main():
    args = parse_arguments()
    args.output_file.parent.mkdir(parents=True, exist_ok=True)
    setup_logging(args.output_file.with_suffix(".log"))
    contig_ids = read_contig_file(args.contig_file) if args.contig_file else args.contig_ids
    cache_dir = None if args.no_cache else args.cache_dir

    # ---------- Load the contigs and assign their categories ----------
    with Database(SCHEMAS[args.database]) as database:
        rules = build_rules(database, cache_dir, args.refresh_cache)
        loader = ContigLoader(database)
        contigs = [loader.fetch_contig(ContigID(contig_id)) for contig_id in contig_ids]
    assign_categories(contigs, rules)
    for contig in contigs:
        category_counts = Counter(feature.category.label for feature in contig.features if feature.category is not None)
        logger.info(
            f"Loaded contig {contig.id} from {args.database}: {contig.length} bp, {len(contig.features)} features"
            + "".join(f", {count} {label}" for label, count in category_counts.items())
        )

    # ---------- Order contigs ----------
    if args.order_by_similarity:
        contigs = order_by_clustering(contigs)
        logger.info(f"Reordered contigs by protein family similarity: {', '.join(str(contig.id) for contig in contigs)}")

    # ---------- Plot contigs ----------
    plot_contigs(contigs, args.output_file)
    logger.info(f"Saved plot to {args.output_file}")


if __name__ == "__main__":
    main()
