from __future__ import annotations

import logging
from dataclasses import dataclass, field
from typing import List, Optional, Literal
from matplotlib.colors import is_color_like

logger = logging.getLogger(__name__)


@dataclass(frozen=True, order=True)
class FeatureID:
    source_id: str
    type: Literal["protein", "pseudo"]

    def __str__(self) -> str:
        return f"{self.source_id}:{self.type}"


@dataclass(frozen=True, order=True)
class ContigID:
    source_id: str
    start: Optional[int] = None
    end: Optional[int] = None

    def __str__(self) -> str:
        if self.start is None:
            return self.source_id
        return f"{self.source_id}:{self.start}-{self.end}"


@dataclass(frozen=True, order=True)
class FamilyID:
    source_id: str

    def __str__(self) -> str:
        return self.source_id
    

@dataclass
class Feature:
    source_id: str
    type: Literal["protein", "pseudo"]
    start: int
    end: int
    strand: Literal["+", "-"]
    protein_family: Optional[FamilyID] = None
    category: Optional[Category] = None

    @property
    def id(self) -> FeatureID:
        return FeatureID(source_id=self.source_id, type=self.type)

    @property
    def wraps(self) -> bool:
        return self.start > self.end


@dataclass(frozen=True)
class Category:
    name: str
    color: str

    def __post_init__(self):
        if not is_color_like(self.color):
            raise ValueError(f"Invalid color code: {self.color}")

@dataclass
class Region:
    id: str
    start: int
    end: int


@dataclass(kw_only=True)
class Contig:
    source_id: str
    start: int
    end: int
    topology: Literal["linear", "circular"]
    features: List[Feature]
    regions: Optional[List[Region]] = None
    category: Optional[Category] = None

    def __post_init__(self):
        if self.start < 1 or self.end < 1:
            raise ValueError(f"Contig {self.id} start and end positions must be positive.")
        if self.start > self.end:
            raise ValueError(f"Contig {self.id} must start before it ends, but runs from {self.start} to {self.end}.")
        if self.topology == "linear" and any(feature.wraps for feature in self.features):
            raise ValueError(f"Linear contig {self.id} cannot have features that wrap around the origin.")

    @property
    def id(self) -> ContigID:
        return ContigID(source_id=self.source_id)

    @property
    def length(self) -> int:
        return self.end - self.start + 1

    def slice(self, start: int, end: int) -> PartialContig:
        """
        Cut out the block from start to end as a linear contig, keeping the original coordinates. Features entirely
        inside the block are kept and features entirely outside it are dropped. Features that cross a border of the
        block are dropped with a warning and kept in border_features. Features that wrap around the origin of a
        circular contig are cut by the origin, so they count as border features whenever they overlap the block.
        Args:
            start (int): The start of the block, at or after the start of this contig.
            end (int): The end of the block, at or after start and at or before the end of this contig.
        Returns:
            PartialContig: The block, with the category of this contig and its features inside the block.
        Raises:
            ValueError: If the block crosses the origin (start after end) or does not lie within this contig.
        """
        if start > end:
            raise ValueError(f"Cannot cut {start}-{end} from contig {self.id}: blocks that cross the origin are not supported yet.")
        if not self.start <= start <= end <= self.end:
            raise ValueError(f"Cannot cut {start}-{end} from contig {self.id}, which runs from {self.start} to {self.end}.")

        features, border_features = [], []
        for feature in self.features:
            if feature.wraps:
                inside = False
                overlaps = feature.start <= end or feature.end >= start
            else:
                inside = start <= feature.start and feature.end <= end
                overlaps = feature.start <= end and feature.end >= start
            if inside:
                features.append(feature)
            elif overlaps:
                border_features.append(feature)

        partial = PartialContig(
            source_id=self.source_id,
            start=start,
            end=end,
            topology="linear",
            features=features,
            regions=[region for region in self.regions if start <= region.start <= region.end <= end] if self.regions is not None else None,
            category=self.category,
            border_features=border_features,
        )
        if border_features:
            logger.warning(
                f"Dropped {len(border_features)} features of contig {self.id} that cross the borders of {partial.id}: "
                + ", ".join(str(feature.id) for feature in border_features)
            )
        return partial

@dataclass(kw_only=True)
class PartialContig(Contig):
    border_features: List[Feature] =field(default_factory=list)

    @property
    def id(self) -> ContigID:
        return ContigID(source_id=self.source_id, start=self.start, end=self.end)
        