import logging
from collections import defaultdict
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, FrozenSet, Iterable, Optional, Sequence, Set, Tuple, Type, Union

from contigviz.cache import load_ids, read_ids
from contigviz.database import Database
from contigviz.models import ID_CLASSES, Category, Contig, ContigID, FamilyID, Feature, ProteinID, PseudoID

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


def assign_categories(contigs: Iterable[Contig], rules: Sequence[CategoryRule]) -> None:
    """
    Set the category of each contig and feature that a rule matches; the first matching rule wins. Contigs and features that no rule matches keep their category. When a later rule with another category also matches, one warning per pair of categories is logged.
    Args:
        contigs (Iterable[Contig]): The contigs, changed in place. Assign before slicing, so blocks keep the categories.
        rules (Sequence[CategoryRule]): The rules, in priority order.
    """
    conflicts: Dict[Tuple[Category, Category], Set[Union[ProteinID, PseudoID, FamilyID, ContigID]]] = defaultdict(set)
    for contig in contigs:
        _assign(contig, (contig.id,), rules, conflicts)
        for feature in contig.features:
            _assign(feature, (feature.id, feature.protein_family), rules, conflicts)
    for (winner, loser), ids in conflicts.items():
        examples = ", ".join(str(id) for id in sorted(ids, key=str)[:5])
        # Categories with the same label differ in color, so show it to tell them apart.
        first, second = (f"'{c.label}' ({c.color})" if winner.label == loser.label else f"'{c.label}'" for c in (winner, loser))
        logger.warning(
            f"{len(ids)} contigs or features match both {first} and {second}; kept {first}, which comes first: "
            f"{examples}{', ...' if len(ids) > 5 else ''}"
        )


def _assign(
        item: Union[Contig, Feature],
        keys: Tuple[Optional[Union[ProteinID, PseudoID, FamilyID, ContigID]], ...],
        rules: Sequence[CategoryRule],
        conflicts: Dict[Tuple[Category, Category], Set[Union[ProteinID, PseudoID, FamilyID, ContigID]]],
    ) -> None:
    """
    Set the category of one contig or feature to that of the first rule whose IDs contain one of its keys.
    Args:
        item (Union[Contig, Feature]): The contig or feature, changed in place.
        keys (Tuple[Optional[Union[ProteinID, PseudoID, FamilyID, ContigID]], ...]): Its IDs to look up: a contig's ID, or a feature's ID and protein family.
        rules (Sequence[CategoryRule]): The rules, in priority order.
        conflicts (Dict[Tuple[Category, Category], Set[Union[ProteinID, PseudoID, FamilyID, ContigID]]]): The conflicts found so far, by pair of categories; added to.
    """
    winner = None
    for rule in rules:
        if any(key in rule.ids for key in keys):
            if winner is None:
                winner = rule.category
            elif rule.category != winner:
                conflicts[(winner, rule.category)].add(item.id)
    if winner is not None:
        item.category = winner

