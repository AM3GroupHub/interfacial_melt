"""Type↔species mapping (vendored dataclass, self-contained).

We do not re-detect species: data/state_index.json already carries the audited
`type_to_species` per state. Only the lookup container is needed here.
"""
from __future__ import annotations
from dataclasses import dataclass


@dataclass
class SpeciesMap:
    type_to_species: dict[int, str]    # {1:"B", 2:"Fe"} or {1:"Fe"} for pure
    mean_v2: dict[int, float]
    v2_ratio_B_over_Fe: float | None
    expected_v2_ratio: float
    method: str

    def fe_type(self) -> int | None:
        for t, s in self.type_to_species.items():
            if s == "Fe":
                return t
        return None

    def b_type(self) -> int | None:
        for t, s in self.type_to_species.items():
            if s == "B":
                return t
        return None
