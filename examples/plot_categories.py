"""
Plot a hybrid plasmid between its two parents, with the three ways to make a category rule: from IDs in the code,
from a file of IDs and from an SQL query.

The plot-contigs command colors one category; with the contigviz API you can combine as many as you like. Download
this script, replication_families.txt and amr_proteins.sql into one folder and run it there on the HPC:

    python plot_categories.py
"""
import logging
from pathlib import Path

from contigviz import (
    Category, CategoryRule, ContigID, Database, FamilyID, ProteinID, assign_category, fetch_contig, plot_contigs,
)
from contigviz.schemas import kes2021

# Show what contigviz does, e.g. how many IDs each rule has and whether the query ran or was read from the cache.
logging.basicConfig(format="%(message)s")
logging.getLogger("contigviz").setLevel(logging.INFO)

CONTIGS = "1991", "2082", "4338"  # parent, hybrid, parent: plotted from top to bottom, so the hybrid is linked to both
OUTPUT = Path("hybrid_plasmid.svg")

# ---------- 1. A rule from IDs in the code: color the contig names by parent or hybrid ----------
parent = CategoryRule(Category("Parent", "#1f77b4"), {ContigID("1991"), ContigID("4338")})
hybrid = CategoryRule(Category("Hybrid", "#d62728"), {ContigID("2082")})

# ---------- 2. A rule from a file of IDs: the protein families of replication initiation proteins ----------
replication = CategoryRule.from_file(Category("Replication", "#5e2bff"), FamilyID, Path("replication_families.txt"))

# ---------- 3. A rule from an SQL query: the proteins with an AMR annotation ----------
with Database(kes2021) as database:
    amr = CategoryRule.from_query(
        database, Category("AMR", "#33ffff"), ProteinID, Path("amr_proteins.sql").read_text(), cache=Path("amr_proteins.txt"),
    )
# --------- 4. Fetch the contigs from the database ----------
    contigs = [fetch_contig(database, ContigID(contig_id)) for contig_id in CONTIGS]

# ---------- 5. Assign the categories to the contigs ----------
for rule in (parent, hybrid, replication, amr):
    assign_category(contigs, rule)

# ---------- 6. Plot the contigs ----------
plot_contigs(contigs, OUTPUT)
print(f"Saved the plot to {OUTPUT}")
