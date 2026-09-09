"""VP35 FAST+FAH: download, extract with a hard stride, score open vs closed.

Not imported by pocket-demo or CI. Do not commit trajectories.
"""

from __future__ import annotations

import json
import tarfile
from pathlib import Path

import requests

from pocket_atlas.cases import load_case
from pocket_atlas.dynamics import overlap_case
from pocket_atlas.dynamics.dyna1 import high_exchange_residues, load_cached_scores
from pocket_atlas.dynamics.labels import resolve_prior
from pocket_atlas.paths import PROCESSED, RAW, ensure_dirs

ZENODO_API = "https://zenodo.org/api/records/{record}"


class TrajectorySkipped(RuntimeError):
    pass


def zenodo_files(record: str) -> list[dict]:
    response = requests.get(ZENODO_API.format(record=record), timeout=60)
    response.raise_for_status()
    return response.json().get("files", [])


def fetch_vp35(
    which: str | None = None,
    dest_dir: Path | None = None,
    force: bool = False,
) -> Path:
    case = load_case("vp35_iid")
    meta = case.raw["zenodo"]
    key = which or meta["default"]
    filename = meta["files"][key]
    ensure_dirs()
    out_dir = dest_dir or (RAW / "vp35")
    out_dir.mkdir(parents=True, exist_ok=True)
    dest = out_dir / filename
    if dest.exists() and not force:
        return dest

    files = zenodo_files(str(meta["record"]))
    match = next((f for f in files if f.get("key") == filename or filename in f.get("key", "")), None)
    if match is None:
        names = [f.get("key") for f in files]
        raise TrajectorySkipped(f"{filename} not in Zenodo {meta['record']}. Available: {names}")
    url = match.get("links", {}).get("self") or match.get("links", {}).get("download")
    if not url:
        raise TrajectorySkipped("Zenodo file has no download URL")
    with requests.get(url, stream=True, timeout=120) as response:
        response.raise_for_status()
        with dest.open("wb") as handle:
            for chunk in response.iter_content(chunk_size=1024 * 1024):
                if chunk:
                    handle.write(chunk)
    return dest


def occupancy_vs_cv(ca_i, ca_j, pocket_open_cutoff: float = 8.0):
    import numpy as np

    i = np.asarray(ca_i)
    j = np.asarray(ca_j)
    dist = np.linalg.norm(i - j, axis=-1)
    return dist, dist > pocket_open_cutoff


def _traj_members(archive: Path) -> list[str]:
    with tarfile.open(archive, "r:*") as tar:
        return [m.name for m in tar.getmembers() if m.isfile()]


def extract_vp35(
    archive: Path,
    dest_dir: Path | None = None,
    stride: int | None = None,
    max_frames: int | None = None,
) -> dict:
    """Extract the tarball and write a strided xtc/pdb pair via mdtraj.

    Hard stride: keep every `stride`-th frame, cap at `max_frames`.
    """
    try:
        import mdtraj as md
    except ImportError as exc:
        raise TrajectorySkipped("mdtraj is required for VP35 extract (`pip install -e '.[md]'`)") from exc

    case = load_case("vp35_iid")
    meta = case.raw["zenodo"]
    stride = int(stride if stride is not None else meta.get("stride", 50))
    max_frames = int(max_frames if max_frames is not None else meta.get("max_frames", 400))
    dest = dest_dir or (PROCESSED / "vp35")
    dest.mkdir(parents=True, exist_ok=True)

    names = _traj_members(archive)
    traj_name = next((n for n in names if n.endswith((".xtc", ".dcd", ".trr"))), None)
    top_name = next((n for n in names if n.endswith((".pdb", ".gro"))), None)
    if traj_name is None or top_name is None:
        raise TrajectorySkipped(f"no traj+topology in {archive.name}. files={names[:20]}")

    extract_root = dest / "extract"
    extract_root.mkdir(parents=True, exist_ok=True)
    wanted = {traj_name, top_name}
    with tarfile.open(archive, "r:*") as tar:
        for member in tar.getmembers():
            if member.name in wanted:
                tar.extract(member, path=extract_root)

    top_path = extract_root / top_name
    traj_path = extract_root / traj_name
    # mdtraj stride is applied at load; then cap frames.
    traj = md.load(str(traj_path), top=str(top_path), stride=stride)
    if len(traj) > max_frames:
        traj = traj[:max_frames]
    out_xtc = dest / "vp35_strided.xtc"
    out_pdb = dest / "vp35_strided.pdb"
    traj.save_xtc(str(out_xtc))
    traj[0].save_pdb(str(out_pdb))
    manifest = {
        "archive": str(archive),
        "stride": stride,
        "max_frames": max_frames,
        "n_frames": int(traj.n_frames),
        "xtc": str(out_xtc),
        "pdb": str(out_pdb),
        "source_traj": traj_name,
        "source_top": top_name,
    }
    (dest / "extract_manifest.json").write_text(json.dumps(manifest, indent=2))
    return manifest


def analyze_vp35(
    xtc: Path | None = None,
    pdb: Path | None = None,
    pocket_open_cutoff: float = 8.0,
) -> dict:
    """CA 225–295 occupancy, open vs closed frames, lining ∩ Dyna-1 if cached."""
    try:
        import mdtraj as md
        import numpy as np
    except ImportError as exc:
        raise TrajectorySkipped("mdtraj is required for VP35 analyze") from exc

    dest = PROCESSED / "vp35"
    xtc = Path(xtc) if xtc else dest / "vp35_strided.xtc"
    pdb = Path(pdb) if pdb else dest / "vp35_strided.pdb"
    if not xtc.exists() or not pdb.exists():
        raise TrajectorySkipped(f"missing strided traj {xtc} or {pdb}. Run extract first.")

    case = load_case("vp35_iid")
    cv = case.raw["open_closed_cv"]
    traj = md.load(str(xtc), top=str(pdb))
    top = traj.topology
    i = _ca_index(top, int(cv["residue_i"]))
    j = _ca_index(top, int(cv["residue_j"]))
    dist, opened = occupancy_vs_cv(
        traj.xyz[:, i] * 10.0,
        traj.xyz[:, j] * 10.0,
        pocket_open_cutoff=pocket_open_cutoff,
    )
    n_open = int(np.sum(opened))
    n_closed = int(np.sum(~opened))

    scores = load_cached_scores(case.name)
    if scores:
        nmr = high_exchange_residues(scores)
        source = "dyna1"
    else:
        nmr, source = resolve_prior(case, prior="dyna1", prefer_dyna1=True)

    protein = {res.resSeq for res in top.residues if res.is_protein}
    overlap = overlap_case(case, protein_residues=protein, nmr_residues=nmr, source=source)
    # mdtraj xyz is nm; occupancy_vs_cv here is called on Å coords.
    result = {
        "n_frames": int(traj.n_frames),
        "cv_residues": [int(cv["residue_i"]), int(cv["residue_j"])],
        "open_cutoff_angstrom": pocket_open_cutoff,
        "n_open": n_open,
        "n_closed": n_closed,
        "open_fraction": float(n_open / len(dist)) if len(dist) else 0.0,
        "cv_mean_angstrom": float(np.mean(dist)),
        "cv_min_angstrom": float(np.min(dist)),
        "cv_max_angstrom": float(np.max(dist)),
        "lining_residues": sorted(case.cryptic_residues),
        "overlap": overlap.as_dict(),
        "prior": source,
        "note": "mdtraj stores nm; CV cutoff default 8 Å. No RelaxDB CPMG for VP35.",
    }
    out = dest / "vp35_occupancy.json"
    dest.mkdir(parents=True, exist_ok=True)
    out.write_text(json.dumps(result, indent=2))
    result["path"] = str(out)
    return result


def _ca_index(topology, resseq: int) -> int:
    for atom in topology.atoms:
        if atom.name == "CA" and atom.residue.resSeq == resseq:
            return atom.index
    raise TrajectorySkipped(f"no CA for residue {resseq} in VP35 topology")
