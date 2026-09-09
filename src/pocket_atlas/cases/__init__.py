from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

import yaml

from pocket_atlas.paths import CASES


@dataclass(frozen=True)
class Case:
    name: str
    raw: dict

    @property
    def chain(self) -> str:
        return str(self.raw.get("chain", "A"))

    @property
    def cryptic_residues(self) -> frozenset[int]:
        return frozenset(int(r) for r in self.raw["cryptic_site"]["residues"])

    @property
    def nmr_residues(self) -> frozenset[int]:
        residues = self.raw.get("nmr_exchange", {}).get("residues", [])
        return frozenset(int(r) for r in residues)

    def residue_set(self, key: str) -> frozenset[int]:
        block = self.raw.get(key, {})
        residues = block.get("residues", []) if isinstance(block, dict) else []
        return frozenset(int(r) for r in residues)


def load_case(name: str) -> Case:
    path = CASES / f"{name}.yaml"
    if not path.exists():
        available = ", ".join(p.stem for p in CASES.glob("*.yaml"))
        raise FileNotFoundError(f"Unknown case {name!r}. Available: {available}")
    with path.open() as handle:
        raw = yaml.safe_load(handle)
    return Case(name=raw["name"], raw=raw)


def list_cases() -> list[str]:
    return sorted(p.stem for p in CASES.glob("*.yaml"))


def case_path(name: str) -> Path:
    return CASES / f"{name}.yaml"
