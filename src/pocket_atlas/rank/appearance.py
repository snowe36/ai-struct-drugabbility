"""VP35 closed crystal vs a high-CV open frame. YAML lining stays out."""

from __future__ import annotations

from pathlib import Path

from pocket_atlas.cases import load_case
from pocket_atlas.dynamics.labels import resolve_prior
from pocket_atlas.io.pdb import read_pdb
from pocket_atlas.paths import PROCESSED, ensure_dirs
from pocket_atlas.pockets import detect_pockets
from pocket_atlas.prepare import prepare_structure, write_pdb
from pocket_atlas.rank import ApoRanking, rank_apo, rank_case
from pocket_atlas.sample import TrajectorySkipped, occupancy_vs_cv


class FrameUnavailable(RuntimeError):
    pass


def find_strided_traj() -> tuple[Path, Path] | None:
    dest = PROCESSED / "vp35"
    xtc = dest / "vp35_strided.xtc"
    pdb = dest / "vp35_strided.pdb"
    if xtc.exists() and pdb.exists():
        return xtc, pdb
    return None


def high_cv_frame_index(xtc: Path, pdb: Path) -> tuple[int, float]:
    try:
        import mdtraj as md
        import numpy as np
    except ImportError as exc:
        raise FrameUnavailable("mdtraj is required for the VP35 open-frame arm") from exc

    from pocket_atlas.sample import _ca_index

    case = load_case("vp35_iid")
    cv = case.raw["open_closed_cv"]
    traj = md.load(str(xtc), top=str(pdb))
    i = _ca_index(traj.topology, int(cv["residue_i"]))
    j = _ca_index(traj.topology, int(cv["residue_j"]))
    dist, _opened = occupancy_vs_cv(traj.xyz[:, i] * 10.0, traj.xyz[:, j] * 10.0)
    idx = int(np.argmax(dist))
    return idx, float(dist[idx])


def structure_from_frame(xtc: Path, pdb: Path, frame: int, chain: str = "A"):
    try:
        import mdtraj as md
    except ImportError as exc:
        raise FrameUnavailable("mdtraj is required for the VP35 open-frame arm") from exc

    ensure_dirs()
    traj = md.load(str(xtc), top=str(pdb))
    if frame < 0 or frame >= traj.n_frames:
        raise FrameUnavailable(f"frame {frame} out of range (n={traj.n_frames})")
    dest = PROCESSED / "vp35"
    dest.mkdir(parents=True, exist_ok=True)
    out = dest / f"vp35_open_frame_{frame}.pdb"
    traj[frame].save_pdb(str(out))
    structure = prepare_structure(read_pdb(out, pdb_id="open"), chain=chain)
    write_pdb(structure, dest / f"vp35_open_frame_{frame}_prep.pdb")
    return structure


def rank_open_frame(prior: str = "dyna1") -> tuple[ApoRanking, dict]:
    pair = find_strided_traj()
    if pair is None:
        raise FrameUnavailable(
            "no strided VP35 xtc at data/processed/vp35/. Keep the 3FKE crystal arm."
        )
    xtc, pdb = pair
    try:
        frame, cv_angstrom = high_cv_frame_index(xtc, pdb)
    except TrajectorySkipped as exc:
        raise FrameUnavailable(str(exc)) from exc
    case = load_case("vp35_iid")
    structure = structure_from_frame(xtc, pdb, frame, chain=case.chain)
    pockets = detect_pockets(structure, chain=case.chain, ligand_mode="exclude")
    protein = set(structure.residue_numbers(chain=case.chain))
    prior_residues, source = resolve_prior(case, prior=prior, prefer_dyna1=prior == "dyna1")
    ranking = rank_apo(
        pockets,
        protein_residues=protein,
        prior_residues=set(prior_residues),
        prior_source=source,
        case="vp35_iid",
        tag="open_frame",
        pdb_id=f"fah_frame_{frame}",
    )
    extra = {
        "xtc": str(xtc),
        "frame": frame,
        "cv_angstrom": round(cv_angstrom, 2),
        "n_pockets": len(pockets),
    }
    ranking.extra.update(extra)
    return ranking, extra


def rank_vp35(prior: str = "dyna1", include_open: bool = True) -> dict:
    """Crystal 3FKE ranking; open-frame arm if the strided xtc is on disk."""
    crystal = rank_case("vp35_iid", prior=prior)
    payload: dict = {
        "crystal": crystal,
        "open": None,
        "open_skipped": None,
    }
    if not include_open:
        payload["open_skipped"] = "open-frame arm disabled"
        return payload
    try:
        opened, extra = rank_open_frame(prior=prior)
        payload["open"] = opened
        payload["open_extra"] = extra
    except FrameUnavailable as exc:
        payload["open_skipped"] = str(exc)
    return payload
