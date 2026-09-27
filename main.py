import argparse
import os
from pathlib import Path
from time import time
from typing import Dict, List, Optional, Tuple

from plot_plasmid.utils import setup_logging
from plot_plasmid.config import load_config
from plot_plasmid.database import ContigLoader
from plot_plasmid.models import Contig

CONTIG_IDS = ["3889", "9353", "10599"]
DEFAULT_DATABASE = "kes2021"
DEFAULT_CONFIG_FILE = Path(__file__).parent / "config.yaml"
PLOT_DIR = "/zfshome/sunam274/compare_plasmids/results/plots"

LOG_DIR = "/zfshome/sunam274/compare_plasmids/results/plotting_data/logs"
os.makedirs(LOG_DIR, exist_ok=True)

def read_contig_file(file_path: str) -> Tuple[List[str], Dict[str, str]]:
    """
    Read contig IDs, and optionally a category per contig, from a file.
    Args:
        file_path (str): The path to a file with one contig ID per line, optionally followed by a tab and a category.
    Returns:
        Tuple[List[str], Dict[str, str]]: The contig IDs in file order, and a dictionary mapping contig ID to category for the contigs that have one.
    """
    contig_ids = []
    contig_categories = {}
    with open(file_path) as fh:
        for line in fh:
            line = line.rstrip("\n")
            if not line.strip():
                continue
            parts = line.split("\t")
            contig_id = parts[0].strip()
            contig_ids.append(contig_id)
            if len(parts) > 1 and parts[1].strip():
                contig_categories[contig_id] = parts[1].strip()
    return contig_ids, contig_categories

def parse_arguments():
    parser = argparse.ArgumentParser(description="Plot plasmid contigs with PyGenomeViz")
    input_group = parser.add_mutually_exclusive_group()
    input_group.add_argument(
        "-c", "--contig_ids", type=str, nargs="+", help="List of contig IDs to plot (space-separated)"
    )
    input_group.add_argument(
        "-f", "--contig_file", type=str, help="Path to a file with one contig ID per line, optionally followed by a tab and a category for label coloring"
    )
    parser.add_argument(
        "-d", "--database", type=str, default=DEFAULT_DATABASE, help=f"Name of the database in the config file to load the contigs from (default: {DEFAULT_DATABASE})"
    )
    parser.add_argument(
        "--config", type=str, default=str(DEFAULT_CONFIG_FILE), help="Path to the YAML config file (default: config.yaml next to this script)"
    )
    parser.add_argument(
        "-o", "--output_file", type=str, default=None, help="Output file path for the plot (e.g. contigs.svg)"
    )
    parser.add_argument(
        "-s", "--order-by-similarity", action="store_true",
        help="Reorder contigs by similarity of their protein families instead of the input order"
    )
    return parser.parse_args()

def run_plot(
        contig_ids: List[str],
        output_file: str,
        database: str = DEFAULT_DATABASE,
        config_file: str = str(DEFAULT_CONFIG_FILE),
        order_by_similarity: bool = False,
        contig_categories: Optional[Dict[str, str]] = None,
    ) -> List[Contig]:
    """
    Load the given contigs from the database and plot them.

    This is `main()`'s body minus argument parsing, so it can be called directly (e.g. from
    another script) with explicit arguments instead of going through argparse/`sys.argv`.
    Args:
        contig_ids (List[str]): A list of contig IDs to plot.
        output_file (str): The path to the output file.
        database (str): The name of the database in the config file to load the contigs from.
        config_file (str): The path to the YAML config file.
        order_by_similarity (bool): Whether to reorder contigs by similarity of their protein families.
        contig_categories (Optional[Dict[str, str]]): A dictionary mapping contig ID to a category label for track-label coloring.
    Returns:
        List[Contig]: The loaded contigs.
    Raises:
        ValueError: If the database is not in the config file.
    """
    contig_categories = contig_categories or {}
    logger = setup_logging(os.path.join(LOG_DIR, f"plot_plasmid_{int(time())}.log"))
    config = load_config(config_file)
    if database not in config.databases:
        raise ValueError(f"Unknown database '{database}'. Choose from: {', '.join(config.databases)} or add it to {config_file}.")

    # ---------- Get contig data ----------
    contigs = []
    with ContigLoader(database, config.databases[database], config.host, config.user) as loader:
        for contig_id in contig_ids:
            contig = loader.fetch_contig(contig_id)
            contig.category = contig_categories.get(contig.id)
            logger.info(f"Loaded contig {contig.id} from {database}: {contig.length} bp, {len(contig.features)} features")
            contigs.append(contig)
    return contigs

def main():
    args = parse_arguments()
    contig_categories: Dict[str, str] = {}
    if args.contig_file:
        contig_ids, contig_categories = read_contig_file(args.contig_file)
    elif args.contig_ids:
        contig_ids = args.contig_ids
    else:
        contig_ids = CONTIG_IDS
    if args.output_file:
        output_file = args.output_file
    else:
        if args.contig_file:
            base = os.path.splitext(os.path.basename(args.contig_file))[0]
        else:
            base = '_'.join(contig_ids)
        output_file = os.path.join(PLOT_DIR, f"{base}.svg")
    run_plot(contig_ids, output_file, args.database, args.config, args.order_by_similarity, contig_categories)

if __name__ == "__main__":
    main()
