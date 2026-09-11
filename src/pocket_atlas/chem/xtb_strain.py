"""Ligand-only GFN2-xTB strain. Not a binding energy. Fail closed without xtb/tblite."""

from __future__ import annotations

import json
import re
import shutil
import subprocess
from pathlib import Path

from pocket_atlas.cases import load_case
from pocket_atlas.dock import HARTREE_TO_KCAL, ligand_atoms, write_ligand_xyz
from pocket_atlas.io.rcsb import load_structure
from pocket_atlas.paths import PROCESSED, REPORTS, ensure_dirs

HARTREE_RE = re.compile(r"TOTAL ENERGY\s+([-+0-9.]+)\s+Eh", re.I)


class XtbUnavailable(RuntimeError):
    pass


def find_xtb() -> str:
    path = shutil.which("xtb")
    if path is None:
        raise XtbUnavailable(
            "xtb not on PATH. Install Grimme xtb (GFN2) or skip strain. "
            "This is ligand strain, not protein DFT."
        )
    return path


def parse_xtb_energy(text: str) -> float:
    match = HARTREE_RE.search(text)
    if match is None:
        raise XtbUnavailable("xtb stdout had no TOTAL ENERGY")
    return float(match.group(1))


def strain_kcal(e_pose: float, e_opt: float) -> float:
    return (e_pose - e_opt) * HARTREE_TO_KCAL


def _xtb_energy(xyz: Path, opt: bool = False, charge: int = 0) -> float:
    xtb = find_xtb()
    cmd = [xtb, str(xyz), "--gfn", "2", "--chrg", str(charge), "--norestart"]
    if opt:
        cmd.extend(["--opt", "crude"])
    result = subprocess.run(cmd, check=True, capture_output=True, text=True, cwd=xyz.parent)
    text = (result.stdout or "") + "\n" + (result.stderr or "")
    return parse_xtb_energy(text)


def strain_xyz(xyz: Path, charge: int = 0) -> dict:
    """Single-point posed ligand vs GFN2-xTB opt. kcal/mol."""
    ensure_dirs()
    e_pose = _xtb_energy(xyz, opt=False, charge=charge)
    e_opt = _xtb_energy(xyz, opt=True, charge=charge)
    kcal = strain_kcal(e_pose, e_opt)
    return {
        "xyz": str(xyz),
        "e_pose_hartree": e_pose,
        "e_opt_hartree": e_opt,
        "strain_kcal": round(kcal, 3),
        "note": "ligand-only GFN2-xTB; not ΔG",
    }


def strain_case(case_name: str, charge: int = 0) -> dict:
    find_xtb()
    case = load_case(case_name)
    holo = case.raw["structures"]["holo"]
    structure = load_structure(holo["pdb_id"])
    atoms = ligand_atoms(structure, holo["ligand_resname"], chain=case.chain)
    dest = PROCESSED / "dock" / case_name
    dest.mkdir(parents=True, exist_ok=True)
    crystal = write_ligand_xyz(atoms, dest / f"{holo['ligand_resname']}_crystal.xyz")
    rows = {"crystal": strain_xyz(crystal, charge=charge)}
    for sdf in sorted(dest.glob("*.xyz")):
        if sdf == crystal:
            continue
        rows[sdf.stem] = strain_xyz(sdf, charge=charge)
    out = dest / "strain.json"
    out.write_text(json.dumps(rows, indent=2))
    (PROCESSED / f"{case_name}_strain.json").write_text(json.dumps(rows, indent=2))
    return {"case": case_name, "path": str(out), "rows": rows}


def write_strain_report(payloads: list[dict], path: Path | None = None) -> Path:
    ensure_dirs()
    path = path or (REPORTS / "strain.md")
    lines = [
        "# Ligand GFN2-xTB strain",
        "",
        "Pose single-point minus gas-phase opt. Not a binding free energy.",
        "",
        "| Case | Pose | Strain (kcal/mol) |",
        "|------|------|-------------------|",
    ]
    for payload in payloads:
        case = payload.get("case", "")
        for name, row in payload.get("rows", {}).items():
            lines.append(f"| {case} | {name} | {row.get('strain_kcal')} |")
    path.write_text("\n".join(lines) + "\n")
    return path
