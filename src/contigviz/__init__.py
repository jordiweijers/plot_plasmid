"""Load contigs and their features from the Dagan lab databases, assign categories to them and plot them with pygenomeviz."""
from importlib import import_module
from importlib.metadata import version
from typing import TYPE_CHECKING

__version__ = version("dagan-contigviz")

# The module each export lives in. It is only imported when the export is first used, so e.g. importing contigviz.models does not load pygenomeviz, scipy and the mysql connector.
_EXPORTS = {
    "CategoryRule": "contigviz.category",
    "assign_category": "contigviz.category",
    "Database": "contigviz.database",
    "fetch_contig": "contigviz.loader",
    "Category": "contigviz.models",
    "Contig": "contigviz.models",
    "ContigID": "contigviz.models",
    "FamilyID": "contigviz.models",
    "Feature": "contigviz.models",
    "PartialContig": "contigviz.models",
    "PartialContigID": "contigviz.models",
    "ProteinID": "contigviz.models",
    "PseudoID": "contigviz.models",
    "Region": "contigviz.models",
    "order_by_clustering": "contigviz.order",
    "plot_contigs": "contigviz.plot",
}

__all__ = list(_EXPORTS)


def __getattr__(name: str):
    if name not in _EXPORTS:
        raise AttributeError(f"module {__name__!r} has no attribute {name!r}")
    value = getattr(import_module(_EXPORTS[name]), name)
    globals()[name] = value  # so __getattr__ is only called once per name
    return value


if TYPE_CHECKING:  # never true when running, but lets IDEs and type checkers see the exports
    from contigviz.category import CategoryRule, assign_category
    from contigviz.database import Database
    from contigviz.loader import fetch_contig
    from contigviz.models import (
        Category, Contig, ContigID, FamilyID, Feature, PartialContig, PartialContigID, ProteinID, PseudoID, Region,
    )
    from contigviz.order import order_by_clustering
    from contigviz.plot import plot_contigs
