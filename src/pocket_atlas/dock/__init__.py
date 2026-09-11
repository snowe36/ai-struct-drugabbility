"""Residue-box GNINA docking. Fail closed if gnina is missing."""

from __future__ import annotations

import json
import math
import shutil
import subprocess
from dataclasses import dataclass, field
from pathlib import Path

import numpy as np

from pocket_atlas.cases import Case, load_case
from pocket_atlas.dynamics.labels import resolve_prior
from pocket_atlas.io.pdb import Atom, Structure
from pocket_atlas.io.rcsb import load_structure
from pocket_atlas.paths import PROCESSED, REPORTS, ensure_dirs
from pocket_atlas.prepare import prepare_structure, write_pdb

BOX_PAD = 8.0
BOX_MIN = 16.0
HARTREE_TO_KCAL = 627.509


class DockUnavailable(RuntimeError):
    pass


@dataclass
class Box:
    name: str
    center: tuple[float, float, float]
    size: tuple[float, float, float]
    residues: list[int]

    def as_dict(self) -> dict:
        return {
            "name": self.name,
            "center": [round(c, 3) for c in self.center],
            "size": [round(s, 3) for s in self.size],
            "residues": self.residues,
        }


@dataclass
class Pose:
    rank: int
    vina: float | None
    cnn_score: float | None
    cnn_affinity: float | None
    rmsd: float | None
    coords: np.ndarray | None = None

    def as_dict(self) -> dict:
        return {
            "rank": self.rank,
            "vina": _round(self.vina),
            "cnn_score": _round(self.cnn_score),
            "cnn_affinity": _round(self.cnn_affinity),
            "rmsd": _round(self.rmsd),
        }


@dataclass
class DockJob:
    case: str
    receptor: str
    box: str
    ligand: str
    poses: list[Pose] = field(default_factory=list)
    sdf: str | None = None

    def best(self) -> Pose | None:
        return self.poses[0] if self.poses else None

    def as_dict(self) -> dict:
        best = self.best()
        return {
            "case": self.case,
            "receptor": self.receptor,
            "box": self.box,
            "ligand": self.ligand,
            "best": None if best is None else best.as_dict(),
            "n_poses": len(self.poses),
            "sdf": self.sdf,
        }


def _round(value: float | None, n: int = 3) -> float | None:
    if value is None or (isinstance(value, float) and not math.isfinite(value)):
        return None
    return round(float(value), n)


def control_key(case: Case) -> str:
    if case.residue_set("catalytic_site"):
        return "catalytic_site"
    return "nucleotide_site"


def residue_box(
    structure: Structure,
    residues: set[int] | frozenset[int],
    name: str,
    chain: str | None = None,
    pad: float = BOX_PAD,
    min_size: float = BOX_MIN,
) -> Box:
    cas = structure.residue_ca(chain=chain)
    pts = np.vstack([cas[r] for r in sorted(residues) if r in cas])
    lo = pts.min(axis=0)
    hi = pts.max(axis=0)
    center = tuple(float(x) for x in (lo + hi) / 2.0)
    extent = hi - lo + 2.0 * pad
    size = tuple(float(max(min_size, x)) for x in extent)
    kept = [r for r in sorted(residues) if r in cas]
    return Box(name=name, center=center, size=size, residues=kept)


def case_boxes(case: Case, structure: Structure, prior: str = "literature") -> dict[str, Box]:
    chain = case.chain
    nmr, _ = resolve_prior(case, prior=prior, prefer_dyna1=prior == "dyna1")
    if not nmr:
        nmr = set(case.nmr_residues)
    boxes = {
        "cryptic": residue_box(structure, case.cryptic_residues, "cryptic", chain=chain),
        "control": residue_box(structure, case.residue_set(control_key(case)), "control", chain=chain),
    }
    if nmr:
        boxes["nmr"] = residue_box(structure, nmr, "nmr", chain=chain)
    return boxes


def ligand_atoms(structure: Structure, resname: str, chain: str | None = None, resseq: int | None = None) -> list[Atom]:
    atoms = structure.ligand_atoms(resname, chain=chain)
    if not atoms:
        atoms = structure.ligand_atoms(resname, chain=None)
    if resseq is not None:
        atoms = [a for a in atoms if a.resseq == resseq]
    if not atoms:
        raise DockUnavailable(f"no ligand {resname} in {structure.pdb_id or structure.path}")
    return atoms


def write_ligand_pdb(atoms: list[Atom], path: Path) -> Path:
    path.parent.mkdir(parents=True, exist_ok=True)
    lines = []
    for i, atom in enumerate(atoms, start=1):
        lines.append(
            f"HETATM{i:5d} {atom.name:>4s} {atom.resname:>3s} {atom.chain:1s}"
            f"{atom.resseq:4d}    {atom.coord[0]:8.3f}{atom.coord[1]:8.3f}{atom.coord[2]:8.3f}"
            f"  1.00  0.00          {atom.element:>2s}"
        )
    lines.append("END")
    path.write_text("\n".join(lines) + "\n")
    return path


def write_ligand_xyz(atoms: list[Atom], path: Path, title: str = "ligand") -> Path:
    path.parent.mkdir(parents=True, exist_ok=True)
    lines = [str(len(atoms)), title]
    for atom in atoms:
        x, y, z = atom.coord
        lines.append(f"{atom.element:2s} {x:12.6f} {y:12.6f} {z:12.6f}")
    path.write_text("\n".join(lines) + "\n")
    return path


def parse_gnina_sdf(text: str) -> list[Pose]:
    """Parse GNINA multi-SDF: vina / CNN tags after each mol block."""
    poses: list[Pose] = []
    blocks = text.split("$$$$")
    for block in blocks:
        block = block.strip()
        if not block:
            continue
        vina = _sdf_tag(block, "minimizedAffinity")
        cnn_score = _sdf_tag(block, "CNNscore")
        cnn_aff = _sdf_tag(block, "CNNaffinity")
        coords = _sdf_coords(block)
        poses.append(
            Pose(
                rank=len(poses) + 1,
                vina=vina,
                cnn_score=cnn_score,
                cnn_affinity=cnn_aff,
                rmsd=None,
                coords=coords,
            )
        )
    return poses


def _sdf_tag(block: str, name: str) -> float | None:
    needle = f"> <{name}>"
    lines = block.splitlines()
    for i, line in enumerate(lines):
        if line.strip() == needle or line.strip().startswith(f"> <{name}"):
            if i + 1 < len(lines):
                try:
                    return float(lines[i + 1].strip())
                except ValueError:
                    return None
    return None


def _sdf_coords(block: str) -> np.ndarray | None:
    lines = block.splitlines()
    counts_i = next((i for i, ln in enumerate(lines) if "V2000" in ln or "V3000" in ln), None)
    if counts_i is None:
        return None
    try:
        n_atoms = int(lines[counts_i][0:3])
    except ValueError:
        return None
    xyz = []
    for line in lines[counts_i + 1 : counts_i + 1 + n_atoms]:
        if len(line) < 31:
            continue
        try:
            xyz.append([float(line[0:10]), float(line[10:20]), float(line[20:30])])
        except ValueError:
            continue
    if not xyz:
        return None
    return np.asarray(xyz, dtype=float)


def rmsd_heavy(ref: np.ndarray, pose: np.ndarray) -> float:
    if ref.size == 0 or pose.size == 0:
        return float("nan")
    n = min(len(ref), len(pose))
    delta = ref[:n] - pose[:n]
    return float(np.sqrt(np.mean(np.sum(delta * delta, axis=1))))


def find_gnina() -> str:
    path = shutil.which("gnina")
    if path is None:
        raise DockUnavailable(
            "gnina not on PATH. On GPU: download the gnina binary from "
            "https://github.com/gnina/gnina/releases and retry."
        )
    return path


def run_gnina(
    receptor: Path,
    ligand: Path,
    box: Box,
    out_sdf: Path,
    exhaustiveness: int = 8,
    num_modes: int = 9,
    seed: int = 1,
) -> list[Pose]:
    gnina = find_gnina()
    out_sdf.parent.mkdir(parents=True, exist_ok=True)
    cmd = [
        gnina,
        "-r", str(receptor),
        "-l", str(ligand),
        "--center_x", f"{box.center[0]:.3f}",
        "--center_y", f"{box.center[1]:.3f}",
        "--center_z", f"{box.center[2]:.3f}",
        "--size_x", f"{box.size[0]:.3f}",
        "--size_y", f"{box.size[1]:.3f}",
        "--size_z", f"{box.size[2]:.3f}",
        "--exhaustiveness", str(exhaustiveness),
        "--num_modes", str(num_modes),
        "--seed", str(seed),
        "-o", str(out_sdf),
    ]
    subprocess.run(cmd, check=True)
    poses = parse_gnina_sdf(out_sdf.read_text())
    if not poses:
        raise DockUnavailable(f"gnina wrote no poses to {out_sdf}")
    return poses


def _receptor_and_ligand(case: Case, tag: str) -> tuple[Structure, list[Atom], str]:
    spec = case.raw["structures"][tag]
    ligand_name = case.raw["structures"]["holo"]["ligand_resname"]
    structure = load_structure(spec["pdb_id"])
    holo_spec = case.raw["structures"]["holo"]
    holo_raw = load_structure(holo_spec["pdb_id"])
    lig = ligand_atoms(
        holo_raw,
        ligand_name,
        chain=case.chain,
        resseq=holo_spec.get("ligand_resseq"),
    )
    receptor = prepare_structure(structure, chain=case.chain, keep_ligand=None)
    return receptor, lig, ligand_name


def dock_case(
    case_name: str,
    prior: str = "literature",
    exhaustiveness: int = 8,
) -> list[DockJob]:
    ensure_dirs()
    find_gnina()
    case = load_case(case_name)
    if "holo" not in case.raw.get("structures", {}):
        raise DockUnavailable(f"{case_name} has no holo ligand (skip VP35)")
    dest = PROCESSED / "dock" / case_name
    dest.mkdir(parents=True, exist_ok=True)

    holo_struct, lig, ligand_name = _receptor_and_ligand(case, "holo")
    lig_pdb = write_ligand_pdb(lig, dest / f"{ligand_name}.pdb")
    write_ligand_xyz(lig, dest / f"{ligand_name}_crystal.xyz", title=f"{ligand_name} crystal")
    ref_xyz = np.vstack([a.coord for a in lig])

    jobs: list[DockJob] = []
    for tag in ("apo", "holo"):
        receptor, _, _ = _receptor_and_ligand(case, tag)
        rec_path = dest / f"{tag}_receptor.pdb"
        write_pdb(receptor, rec_path)
        boxes = case_boxes(case, receptor, prior=prior)
        (dest / f"{tag}_boxes.json").write_text(
            json.dumps({k: v.as_dict() for k, v in boxes.items()}, indent=2)
        )
        for box_name, box in boxes.items():
            sdf = dest / f"{tag}_{box_name}.sdf"
            poses = run_gnina(rec_path, lig_pdb, box, sdf, exhaustiveness=exhaustiveness)
            for pose in poses:
                if pose.coords is not None:
                    pose.rmsd = rmsd_heavy(ref_xyz, pose.coords)
            best = poses[0]
            if best.coords is not None:
                elements = [a.element for a in lig[: len(best.coords)]]
                _write_pose_xyz(
                    elements,
                    best.coords,
                    dest / f"{tag}_{box_name}_best.xyz",
                    title=f"{case_name} {tag} {box_name}",
                )
            job = DockJob(
                case=case_name,
                receptor=tag,
                box=box_name,
                ligand=ligand_name,
                poses=poses,
                sdf=str(sdf),
            )
            jobs.append(job)
    (PROCESSED / f"{case_name}_dock.json").write_text(
        json.dumps([j.as_dict() for j in jobs], indent=2)
    )
    return jobs


def write_dock_report(jobs: list[DockJob], path: Path | None = None) -> Path:
    ensure_dirs()
    path = path or (REPORTS / "dock.md")
    lines = [
        "# GNINA recovery docking",
        "",
        "Boxes are residue centroids (cryptic lining, NMR prior, catalytic/nucleotide control).",
        "Sotorasib is covalent to Cys12; RMSD to crystal is the metric, not ΔG.",
        "",
        "| Case | Receptor | Box | Ligand | Vina | CNN aff. | RMSD (Å) |",
        "|------|----------|-----|--------|------|----------|----------|",
    ]
    payload = []
    for job in jobs:
        best = job.best()
        vina = "" if best is None else best.vina
        cnn = "" if best is None else best.cnn_affinity
        rmsd = "" if best is None else best.rmsd
        lines.append(
            f"| {job.case} | {job.receptor} | {job.box} | {job.ligand} "
            f"| {vina} | {cnn} | {rmsd} |"
        )
        payload.append(job.as_dict())
    path.write_text("\n".join(lines) + "\n")
    (PROCESSED / "dock_summary.json").write_text(json.dumps(payload, indent=2))
    return path


def _write_pose_xyz(elements: list[str], coords: np.ndarray, path: Path, title: str) -> Path:
    path.parent.mkdir(parents=True, exist_ok=True)
    n = min(len(elements), len(coords))
    lines = [str(n), title]
    for element, xyz in zip(elements[:n], coords[:n], strict=True):
        lines.append(f"{element:2s} {xyz[0]:12.6f} {xyz[1]:12.6f} {xyz[2]:12.6f}")
    path.write_text("\n".join(lines) + "\n")
    return path
