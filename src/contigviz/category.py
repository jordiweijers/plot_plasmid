import logging
from dataclasses import dataclass
from pathlib import Path
from typing import FrozenSet, Iterable, Optional, Type, Union

from contigviz.cache import load_ids, read_ids
from contigviz.database import Database
from contigviz.models import ID_CLASSES, Category, Contig, ContigID, FamilyID, ProteinID, PseudoID

logger = logging.getLogger(__name__)


@dataclass(frozen=True)
class CategoryRule:
    """
    Assigns a category to the contigs and features with the given IDs. To get the IDs from an SQL query, make the rule with CategoryRule.from_query, and from a file with CategoryRule.from_file.
    Attributes:
        category (Category): The category to assign. Several rules can share one category.
        ids (FrozenSet[Union[ProteinID, PseudoID, ContigID, FamilyID]]): The IDs, which may be of different kinds: a ProteinID or PseudoID matches a feature, a FamilyID every feature in that protein family, and a ContigID a contig. A list or set is turned into a frozenset, so looking up an ID stays fast.
    """
    category: Category
    ids: FrozenSet[Union[ProteinID, PseudoID, ContigID, FamilyID]]

    def __post_init__(self):
        # A frozen dataclass can only set its fields through object.__setattr__.
        object.__setattr__(self, "ids", frozenset(self.ids))
        wrong = [id for id in self.ids if type(id) not in ID_CLASSES.values()]
        if wrong:
            raise TypeError(
                f"The IDs of category rule '{self.category.label}' must be one of {', '.join(id_class.__name__ for id_class in ID_CLASSES.values())}, but include {wrong[0]!r}."
            )

    @classmethod
    def from_query(
            cls,
            database: Database,
            category: Category,
            id_class: Type[Union[ProteinID, PseudoID, ContigID, FamilyID]],
            query: str,
            cache: Optional[Path] = None,
            refresh: bool = False,
        ) -> "CategoryRule":
        """
        Make a rule from the IDs an SQL query returns.
        Args:
            database (Database): The database to run the query on.
            category (Category): The category to assign.
            id_class (Type[Union[ProteinID, PseudoID, ContigID, FamilyID]]): What the IDs the query returns refer to, e.g. ProteinID.
            query (str): An SQL query that returns one column of IDs from the whole database.
            cache (Optional[Path]): The Path where to save the IDs in, or None to always run the query.
            refresh (bool): Whether to run the query even if its saved IDs are up to date.
        Returns:
            CategoryRule: The rule.
        Raises:
            ValueError: If id_class is not one of the ID classes, the query uses @contig_id or does not return one column.
            RuntimeError: If MySQL fails to run the query.
        """
        if id_class not in ID_CLASSES.values():
            raise ValueError(f"Category '{category.label}' has IDs of class {id_class!r}, which must be one of {', '.join(id_class.__name__ for id_class in ID_CLASSES.values())}.")
        try:
            source_ids = load_ids(database, query, cache, refresh)
        except (RuntimeError, ValueError) as e:
            e.add_note(f"In category '{category.label}'")
            raise
        if not source_ids:
            logger.warning(f"Category '{category.label}': The query returned no IDs, so the rule matches nothing")
        rule = cls(category=category, ids=frozenset(id_class(source_id) for source_id in source_ids))
        logger.info(f"Category '{category.label}': {len(rule.ids)} IDs from its query")
        return rule

    @classmethod
    def from_file(
            cls,
            category: Category,
            id_class: Type[Union[ProteinID, PseudoID, ContigID, FamilyID]],
            path: Path,
        ) -> "CategoryRule":
        """
        Make a rule from the IDs in a file.
        Args:
            category (Category): The category to assign.
            id_class (Type[Union[ProteinID, PseudoID, ContigID, FamilyID]]): What the IDs in the file refer to, e.g. ProteinID.
            path (Path): A file with one ID per line. Blank lines and lines starting with # are skipped, so a cache
                file of from_query can be read as well.
        Returns:
            CategoryRule: The rule.
        Raises:
            ValueError: If id_class is not one of the ID classes.
            FileNotFoundError: If the file does not exist.
        """
        if id_class not in ID_CLASSES.values():
            raise ValueError(f"Category '{category.label}' has IDs of class {id_class!r}, which must be one of {', '.join(id_class.__name__ for id_class in ID_CLASSES.values())}.")
        source_ids = read_ids(path)
        if not source_ids:
            logger.warning(f"Category '{category.label}': {path} contains no IDs, so the rule matches nothing")
        rule = cls(category=category, ids=frozenset(id_class(source_id) for source_id in source_ids))
        logger.info(f"Category '{category.label}': {len(rule.ids)} IDs from {path}")
        return rule


def assign_category(contigs: Iterable[Contig], rule: CategoryRule) -> None:
    """
    Set the category of contigs and features according to the rule. If a contig or a feature already has a category, it will be overwritten with a warning.
    Args:
        contigs (Iterable[Contig]): The contigs, changed in place. Assign before slicing, so blocks keep the categories.
        rule (CategoryRule): The rule to assign.
    """
    overwritten = []
    for contig in contigs:
        items = [(contig, (contig.id,))] + [(feature, (feature.id, feature.protein_family)) for feature in contig.features]
        for item, keys in items:
            if any(key in rule.ids for key in keys):
                if item.category is not None and item.category != rule.category:
                    overwritten.append(f"{item.id} ({item.category.label})")
                item.category = rule.category
    if overwritten:
        logger.warning(
            f"Category '{rule.category.label}' overwrote the category of {len(overwritten)} contigs or features: "
            f"{', '.join(overwritten[:5])}{', ...' if len(overwritten) > 5 else ''}"
        )

