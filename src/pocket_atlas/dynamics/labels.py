"""Residue priors: literature YAML, RelaxDB-CPMG, or cached Dyna-1."""

from __future__ import annotations

from pathlib import Path

import yaml

from pocket_atlas.cases import Case
from pocket_atlas.dynamics.dyna1 import high_exchange_residues, load_cached_scores

PRIORS = ("literature", "relaxdb", "dyna1")
RESOURCES = Path(__file__).resolve().parents[1] / "resources" / "relaxdb_cpmg.yaml"


def load_relaxdb_bundle() -> dict:
    with RESOURCES.open() as handle:
        return yaml.safe_load(handle)


def cpmg_entry(entry_id: str) -> dict:
    bundle = load_relaxdb_bundle()
    entry = bundle["entries"].get(entry_id)
    if entry is None:
        raise KeyError(f"unknown CPMG entry {entry_id!r}")
    return entry


def cpmg_exchange_residues(entry_id: str) -> frozenset[int]:
    """1-based indices into the deposited CPMG sequence (X/Y)."""
    return frozenset(int(r) for r in cpmg_entry(entry_id)["exchange_residues"])


def list_cpmg_entries() -> list[str]:
    return sorted(load_relaxdb_bundle()["entries"])


def relaxdb_residues(case_name: str) -> frozenset[int] | None:
    """Official RelaxDB-CPMG exchange residues, or None if this case has none.

    TEM-1 is not in RelaxDB-CPMG. The BLAC entry is Mtb BlaC (P9WKD3).
    """
    bundle = load_relaxdb_bundle()
    entry_id = bundle.get("case_map", {}).get(case_name)
    if not entry_id:
        return None
    entry = bundle["entries"][entry_id]
    return frozenset(int(r) for r in entry["exchange_residues"])


def resolve_prior(
    case: Case,
    prior: str = "literature",
    prefer_dyna1: bool = False,
) -> tuple[set[int], str]:
    """Return (residues, source). Unknown / missing priors fall back to literature."""
    if prior not in PRIORS:
        raise ValueError(f"prior must be one of {PRIORS}, got {prior!r}")

    if prefer_dyna1 or prior == "dyna1":
        scores = load_cached_scores(case.name)
        if scores:
            return high_exchange_residues(scores), "dyna1"
        if prior == "dyna1":
            # explicit dyna1 request with no cache → literature, not silent relaxdb
            return set(case.nmr_residues), "nmr_literature"

    if prior == "relaxdb":
        residues = relaxdb_residues(case.name)
        if residues is not None:
            return set(residues), "relaxdb_cpmg"
        return set(case.nmr_residues), "nmr_literature"

    return set(case.nmr_residues), "nmr_literature"
