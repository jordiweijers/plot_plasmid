from __future__ import annotations

from dataclasses import dataclass
from typing import List, Optional, Literal, Tuple


@dataclass
class Feature:
    id: str
    type: Literal["protein", "pseudo"]
    start: int
    end: int
    strand: Literal["+", "-"]
    protein_family: Optional[str] = None
    category: Optional[str] = None

    @property
    def key(self) -> Tuple[str, str]:
        return (self.id, self.type)


@dataclass
class Region:
    id: str
    start: int
    end: int


@dataclass
class Contig:
    id: str
    length: int
    features: List[Feature]
    regions: Optional[List[Region]] = None
    category: Optional[str] = None
