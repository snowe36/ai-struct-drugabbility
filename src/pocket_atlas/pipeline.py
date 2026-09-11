from __future__ import annotations

import json
from dataclasses import dataclass, field
from pathlib import Path

from pocket_atlas.cases import Case, load_case
from pocket_atlas.dynamics import OverlapResult, overlap_case
from pocket_atlas.dynamics.labels import resolve_prior
from pocket_atlas.io.pdb import Structure
from pocket_atlas.io.rcsb import load_structure
from pocket_atlas.paths import PROCESSED, ensure_dirs
from pocket_atlas.pockets import DETECTOR_VERSION, Pocket, detect_pockets, pocket_near_residues
from pocket_atlas.prepare import prepare_structure, write_pdb


def _site_clearance(pocket: Pocket | None) -> float:
    if pocket is None:
        return 0.0
    return float(pocket.extra.get("seed_clearance", 0.0))


@dataclass
class ArmResult:
    case: str
    tag: str  # apo | holo
    pdb_id: str
    n_residues: int
    pockets: list[Pocket]
    site_pocket: Pocket | None
    site_found: bool
    ligand_mode: str = "exclude"

    def as_dict(self) -> dict:
        return {
            "case": self.case,
            "tag": self.tag,
            "pdb_id": self.pdb_id,
            "n_residues": self.n_residues,
            "n_pockets": len(self.pockets),
            "site_found": self.site_found,
            "ligand_mode": self.ligand_mode,
            "detector_version": DETECTOR_VERSION,
            "site_clearance": round(_site_clearance(self.site_pocket), 3),
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
            "detector_version": DETECTOR_VERSION,
            "apo": self.apo.as_dict(),
            "holo": self.holo.as_dict(),
            "overlap": self.overlap.as_dict(),
            "cryptic_in_apo": self.cryptic_in_apo,
            "cryptic_in_holo": self.cryptic_in_holo,
            "scores_source": self.scores_source,
            "extra": {k: v for k, v in self.extra.items() if k != "nmr_residues"},
        }


def _run_arm(case: Case, tag: str, ligand_mode: str = "exclude") -> tuple[ArmResult, Structure]:
    spec = case.raw["structures"][tag]
    pdb_id = spec["pdb_id"]
    chain = case.chain
    ligand = spec.get("ligand_resname") if tag == "holo" else None
    structure = load_structure(pdb_id)
    prepared = prepare_structure(structure, chain=chain, keep_ligand=ligand)
    out = PROCESSED / f"{case.name}_{tag}_{pdb_id}.pdb"
    write_pdb(prepared, out)
    mode = ligand_mode if tag == "holo" else "exclude"
    pockets = detect_pockets(prepared, chain=chain, ligand_mode=mode)
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
            ligand_mode=mode,
        ),
        prepared,
    )


def run_campaign(
    case_name: str,
    prefer_dyna1: bool = True,
    prior: str = "literature",
    holo_ligand_mode: str = "exclude",
    both_holo_modes: bool = True,
) -> Campaign:
    ensure_dirs()
    case = load_case(case_name)
    apo, apo_struct = _run_arm(case, "apo", ligand_mode="exclude")
    holo, holo_struct = _run_arm(case, "holo", ligand_mode=holo_ligand_mode)

    protein = set(holo_struct.residue_numbers(chain=case.chain))
    nmr, source = resolve_prior(case, prior=prior, prefer_dyna1=prefer_dyna1)
    overlap = overlap_case(case, protein_residues=protein, nmr_residues=nmr, source=source)

    extra: dict = {
        "nmr_residues": sorted(nmr),
        "detector_version": DETECTOR_VERSION,
        "prior": prior,
        "apo_clearance": _site_clearance(apo.site_pocket),
        "holo_clearance": _site_clearance(holo.site_pocket),
    }
    if both_holo_modes:
        other = "include" if holo_ligand_mode == "exclude" else "exclude"
        alt, _ = _run_arm(case, "holo", ligand_mode=other)
        extra["holo_modes"] = {
            holo.ligand_mode: holo.as_dict()["site_pocket"],
            alt.ligand_mode: alt.as_dict()["site_pocket"],
        }

    return Campaign(
        case=case,
        apo=apo,
        holo=holo,
        overlap=overlap,
        cryptic_in_apo=apo.site_found,
        cryptic_in_holo=holo.site_found,
        scores_source=source,
        extra=extra,
    )


def write_campaign_json(campaign: Campaign, path: Path | None = None) -> Path:
    ensure_dirs()
    path = path or (PROCESSED / f"{campaign.case.name}_campaign.json")
    path.write_text(json.dumps(campaign.as_dict(), indent=2))
    return path


def prepare_crystal(case_name: str) -> Path:
    """Write a prepared PDB for crystal-only cases (VP35 3FKE)."""
    ensure_dirs()
    case = load_case(case_name)
    structures = case.raw["structures"]
    key = "crystal" if "crystal" in structures else "apo"
    spec = structures[key]
    prepared = prepare_structure(load_structure(spec["pdb_id"]), chain=case.chain)
    out = PROCESSED / f"{case.name}_{key}_{spec['pdb_id']}.pdb"
    write_pdb(prepared, out)
    return out
