"""Download public VP35 FAH/FAST trajectories from Zenodo.

These archives are multi-gigabyte. CI and pocket-demo never call this.
Stride hard; never commit tarballs.
"""

from __future__ import annotations

from pathlib import Path

import requests

from pocket_atlas.cases import load_case
from pocket_atlas.paths import RAW, ensure_dirs

ZENODO_API = "https://zenodo.org/api/records/{record}"


class TrajectorySkipped(RuntimeError):
    pass


def zenodo_files(record: str) -> list[dict]:
    response = requests.get(ZENODO_API.format(record=record), timeout=60)
    response.raise_for_status()
    files = response.json().get("files", [])
    return files


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
        raise TrajectorySkipped(
            f"{filename} not in Zenodo {meta['record']}. Available: {names}"
        )
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


def occupancy_vs_cv(
    ca_i,
    ca_j,
    pocket_open_cutoff: float = 8.0,
):
    """Per-frame open/closed from the 225–295 CA distance.

    Mallimadugula et al. use this CV. We only compute it when mdtraj has
    already stacked the coordinates — this helper stays numpy-only.
    """
    import numpy as np

    i = np.asarray(ca_i)
    j = np.asarray(ca_j)
    dist = np.linalg.norm(i - j, axis=-1)
    return dist, dist > pocket_open_cutoff
