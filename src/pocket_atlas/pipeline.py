from __future__ import annotations

import json
from dataclasses import dataclass, field
from pathlib import Path

from pocket_atlas.cases import Case, load_case
from pocket_atlas.dynamics import OverlapResult, overlap_case
from pocket_atlas.dynamics.dyna1 import high_exchange_residues, load_cached_scores
from pocket_atlas.io.pdb import Structure
from pocket_atlas.io.rcsb import load_structure
from pocket_atlas.paths import PROCESSED, ensure_dirs
from pocket_atlas.pockets import Pocket, detect_pockets, pocket_near_residues
from pocket_atlas.prepare import prepare_structure, write_pdb


@dataclass
class ArmResult:
    case: str
    tag: str  # apo | holo
    pdb_id: str
    n_residues: int
    pockets: list[Pocket]
    site_pocket: Pocket | None
    site_found: bool

    def as_dict(self) -> dict:
        return {
            "case": self.case,
            "tag": self.tag,
            "pdb_id": self.pdb_id,
            "n_residues": self.n_residues,
            "n_pockets": len(self.pockets),
            "site_found": self.site_found,
            "site_pocket": None if self.site_pocket is None else self.site_pocket.as_dict(),
            "top_pockets": [p.as_dict() for p in self.pockets[:5]],
        }


@dataclass
class Campaign:
    case: Case
    apo: ArmResult
    holo: ArmResult
    overlap: OverlapResult
    cryptic_in_apo: bool
    cryptic_in_holo: bool
    scores_source: str
    extra: dict = field(default_factory=dict)

    def as_dict(self) -> dict:
        return {
            "case": self.case.name,
            "title": self.case.raw.get("title"),
            "apo": self.apo.as_dict(),
            "holo": self.holo.as_dict(),
            "overlap": self.overlap.as_dict(),
            "cryptic_in_apo": self.cryptic_in_apo,
            "cryptic_in_holo": self.cryptic_in_holo,
            "scores_source": self.scores_source,
        }


def _run_arm(case: Case, tag: str) -> tuple[ArmResult, Structure]:
    spec = case.raw["structures"][tag]
    pdb_id = spec["pdb_id"]
    chain = case.chain
    ligand = spec.get("ligand_resname") if tag == "holo" else None
    structure = load_structure(pdb_id)
    prepared = prepare_structure(structure, chain=chain, keep_ligand=ligand)
    out = PROCESSED / f"{case.name}_{tag}_{pdb_id}.pdb"
    write_pdb(prepared, out)
    pockets = detect_pockets(prepared, chain=chain)
    site = pocket_near_residues(pockets, set(case.cryptic_residues))
    return (
        ArmResult(
            case=case.name,
            tag=tag,
            pdb_id=pdb_id,
            n_residues=len(prepared.residue_numbers(chain=chain)),
            pockets=pockets,
            site_pocket=site,
            site_found=site is not None,
        ),
        prepared,
    )


def run_campaign(case_name: str, prefer_dyna1: bool = True) -> Campaign:
    ensure_dirs()
    case = load_case(case_name)
    apo, apo_struct = _run_arm(case, "apo")
    holo, holo_struct = _run_arm(case, "holo")

    protein = set(holo_struct.residue_numbers(chain=case.chain))
    dyna_scores = load_cached_scores(case.name) if prefer_dyna1 else None
    if dyna_scores:
        nmr = high_exchange_residues(dyna_scores)
        source = "dyna1_cached"
    else:
        nmr = set(case.nmr_residues)
        source = "nmr_literature"

    overlap = overlap_case(case, protein_residues=protein, nmr_residues=nmr, source=source)
    return Campaign(
        case=case,
        apo=apo,
        holo=holo,
        overlap=overlap,
        cryptic_in_apo=apo.site_found,
        cryptic_in_holo=holo.site_found,
        scores_source=source,
        extra={"nmr_residues": sorted(nmr)},
    )


def write_campaign_json(campaign: Campaign, path: Path | None = None) -> Path:
    ensure_dirs()
    path = path or (PROCESSED / f"{campaign.case.name}_campaign.json")
    path.write_text(json.dumps(campaign.as_dict(), indent=2))
    return path
