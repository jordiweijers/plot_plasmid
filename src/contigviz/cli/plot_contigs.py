import argparse
import logging
from collections import Counter
from pathlib import Path
from typing import List, Optional, Literal

from matplotlib.colors import is_color_like

from contigviz import (
    Category, CategoryRule, Database, ProteinID, PseudoID, FamilyID, ContigID, 
    fetch_contig, assign_category, order_by_clustering, plot_contigs,
)
from contigviz.schemas import SCHEMAS
from contigviz.models import ID_CLASSES


logger = logging.getLogger(__name__)

CATEGORY_COLOR = "#d62728"


def parse_arguments():
    parser = argparse.ArgumentParser(description="Plot contigs from one database with pygenomeviz, optionally coloring the features or contigs of one category")
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
    category_group = parser.add_argument_group("Catategory", "Color the features or contigs based on their IDs")
    category_group.add_argument(
        "--category", type=Path, metavar="FILE",
        help="A .sql query file that returns one column of IDs from the whole database, or a file with one ID per line"
    )
    category_group.add_argument(
        "--id-type", choices=sorted(ID_CLASSES),
        help="What the IDs refer to; a family ID colors every feature in that protein family (required with --category)"
    )
    category_group.add_argument(
        "--label", help="Name of the category in the legend (default: the file name without extension)"
    )
    category_group.add_argument(
        "--color", help=f"Color of the category (default: {CATEGORY_COLOR})"
    )
    category_group.add_argument(
        "--cache", type=Path, metavar="FILE",
        help="File to cache the IDs of a .sql category in, so the query only runs again when the database changes (default: None)"
    )
    category_group.add_argument(
        "--refresh-cache", action="store_true",
        help="Rerun the query and overwrite the IDs in --cache, even if they are up to date"
    )

    args = parser.parse_args()
    if args.category is None:
        given = [flag for flag, value in (("--id-type", args.id_type), ("--label", args.label), ("--color", args.color), ("--cache", args.cache)) if value is not None]
        if given:
            parser.error(f"{', '.join(given)} only works together with --category")
    elif args.id_type is None:
        parser.error("--category needs --id-type to know what its IDs refer to")
    elif args.cache is not None and args.category.suffix != ".sql":
        parser.error("--cache only works with a .sql category; a file of IDs is read directly")
    if args.refresh_cache and args.cache is None:
        parser.error("--refresh-cache only works together with --cache")
    if args.color is not None and not is_color_like(args.color):
        parser.error(f"argument --color: {args.color!r} is not a color")
    return args


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
    logging.basicConfig(level=logging.WARNING, handlers=[console_handler, file_handler])
    for name in (__name__, "contigviz"):
        logging.getLogger(name).setLevel(logging.DEBUG)


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


def build_rule(database: Database, category: Category, id_type: Literal["protein", "pseudo", "family", "contig"], path: Path, cache: Optional[Path] = None, refresh: bool = False) -> CategoryRule:
    """
    Make the rule of one category from a file: a .sql file is run as a query on the database, any other file is read as a list of IDs.
    Args:
        database (Database): The database to run the query on.
        category (Category): The category to assign.
        id_type (str): What the IDs refer to, one of the keys of ID_CLASSES, e.g. "protein".
        path (Path): The .sql query or the file with one ID per line.
        cache (Optional[Path]): The file to save the IDs of a query in, or None to always run the query.
        refresh (bool): Whether to run the query even if its saved IDs are up to date.
    Returns:
        CategoryRule: The rule.
    """
    id_class = ID_CLASSES[id_type]
    if path.suffix != ".sql":
        return CategoryRule.from_file(category, id_class, path)
    return CategoryRule.from_query(database, category, id_class, path.read_text(), cache, refresh)


def main():
    # ---------- Parse command line arguments and setup logging ----------
    args = parse_arguments()
    args.output_file.parent.mkdir(parents=True, exist_ok=True)
    setup_logging(args.output_file.with_suffix(".log"))
    contig_ids = read_contig_file(args.contig_file) if args.contig_file else args.contig_ids

    # ---------- Load the contigs and build the category rule ----------
    with Database(SCHEMAS[args.database]) as database:
        contigs = [fetch_contig(database,ContigID(contig_id)) for contig_id in contig_ids]
        rule = None
        if args.category is not None:
            category = Category(args.label or args.category.stem, args.color or CATEGORY_COLOR)
            rule = build_rule(database, category, args.id_type, args.category, args.cache, args.refresh_cache)

    # ---------- Assign the category to contigs ----------
    if rule is not None:
        assign_category(contigs, rule)
    for contig in contigs:
        category_counts = Counter(feature.category.label for feature in contig.features if feature.category is not None)
        logger.info(f"Loaded contig {contig.id} from {args.database}: {contig.length} bp, {len(contig.features)} features" + "".join(f", {count} {label}" for label, count in category_counts.items()))

    # ---------- Order contigs ----------
    if args.order_by_similarity:
        contigs = order_by_clustering(contigs)
        logger.info(f"Reordered contigs by protein family similarity: {', '.join(str(contig.id) for contig in contigs)}")

    # ---------- Plot contigs ----------
    plot_contigs(contigs, args.output_file)
    logger.info(f"Saved plot to {args.output_file}")


if __name__ == "__main__":
    main()
