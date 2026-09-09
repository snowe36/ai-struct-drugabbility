from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from pocket_atlas.cases import Case


@dataclass
class OverlapResult:
    n_protein: int
    n_cryptic: int
    n_nmr: int
    n_control: int
    cryptic_and_nmr: int
    control_and_nmr: int
    cryptic_enrichment: float
    control_enrichment: float
    odds_ratio: float
    source: str

    def as_dict(self) -> dict:
        return {
            "n_protein": self.n_protein,
            "n_cryptic": self.n_cryptic,
            "n_nmr": self.n_nmr,
            "n_control": self.n_control,
            "cryptic_and_nmr": self.cryptic_and_nmr,
            "control_and_nmr": self.control_and_nmr,
            "cryptic_enrichment": round(self.cryptic_enrichment, 3),
            "control_enrichment": round(self.control_enrichment, 3),
            "odds_ratio": round(self.odds_ratio, 3) if np.isfinite(self.odds_ratio) else None,
            "source": self.source,
        }


def _fraction(part: int, whole: int) -> float:
    if whole == 0:
        return 0.0
    return part / whole


def overlap_case(
    case: Case,
    protein_residues: set[int],
    nmr_residues: set[int] | None = None,
    control_key: str = "catalytic_site",
    source: str = "nmr_literature",
) -> OverlapResult:
    """Ask whether NMR-dynamic residues enrich in the cryptic lining.

    Control is the catalytic or nucleotide site — the obvious 'functional'
    residues that are *not* the cryptic pocket. Enrichment > 1 means the
    NMR prior lights the cryptic lining more than that control.
    """
    cryptic = set(case.cryptic_residues) & protein_residues
    if nmr_residues is None:
        nmr_residues = set(case.nmr_residues)
    nmr = set(nmr_residues) & protein_residues
    control = set(case.residue_set(control_key)) & protein_residues
    if not control:
        # KRAS uses nucleotide_site rather than catalytic_site.
        for key in ("nucleotide_site", "catalytic_site"):
            control = set(case.residue_set(key)) & protein_residues
            if control:
                break

    n = len(protein_residues)
    cryptic_and_nmr = len(cryptic & nmr)
    control_and_nmr = len(control & nmr)
    cryptic_frac = _fraction(cryptic_and_nmr, len(cryptic))
    control_frac = _fraction(control_and_nmr, len(control))
    background = _fraction(len(nmr), n)
    cryptic_enr = cryptic_frac / background if background else 0.0
    control_enr = control_frac / background if background else 0.0

    # Odds of NMR-positive given cryptic vs given control.
    a = cryptic_and_nmr
    b = len(cryptic) - cryptic_and_nmr
    c = control_and_nmr
    d = len(control) - control_and_nmr
    if min(a + b, c + d) == 0 or (b * c) == 0:
        odds = float("inf") if a > 0 and c == 0 else 0.0
    else:
        odds = (a / b) / (c / d) if b and d else float("inf")

    return OverlapResult(
        n_protein=n,
        n_cryptic=len(cryptic),
        n_nmr=len(nmr),
        n_control=len(control),
        cryptic_and_nmr=cryptic_and_nmr,
        control_and_nmr=control_and_nmr,
        cryptic_enrichment=cryptic_enr,
        control_enrichment=control_enr,
        odds_ratio=float(odds),
        source=source,
    )


def nmr_scores_from_labels(residues: set[int], protein_residues: set[int]) -> dict[int, float]:
    """Binary literature labels as a stand-in when Dyna-1 weights are absent.

    This is *not* Dyna-1. Captions must say so.
    """
    return {r: (1.0 if r in residues else 0.0) for r in protein_residues}
