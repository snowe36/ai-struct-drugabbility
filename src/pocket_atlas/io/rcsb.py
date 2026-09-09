from __future__ import annotations

from pathlib import Path

import requests

from pocket_atlas.io.pdb import Structure, parse_pdb, read_pdb
from pocket_atlas.paths import RAW, ensure_dirs

RCSB_PDB = "https://files.rcsb.org/download/{pdb_id}.pdb"


def fetch_pdb(pdb_id: str, dest_dir: Path | None = None, force: bool = False) -> Path:
    """Download a PDB from RCSB. Cached under data/raw/."""
    ensure_dirs()
    pdb_id = pdb_id.upper()
    out_dir = dest_dir or RAW
    out_dir.mkdir(parents=True, exist_ok=True)
    path = out_dir / f"{pdb_id}.pdb"
    if path.exists() and not force:
        return path
    url = RCSB_PDB.format(pdb_id=pdb_id)
    response = requests.get(url, timeout=60)
    response.raise_for_status()
    if "HEADER" not in response.text[:200] and "ATOM" not in response.text:
        raise RuntimeError(f"RCSB did not return a PDB for {pdb_id}")
    path.write_text(response.text)
    return path


def load_structure(pdb_id: str, dest_dir: Path | None = None) -> Structure:
    path = fetch_pdb(pdb_id, dest_dir=dest_dir)
    return read_pdb(path, pdb_id=pdb_id)


def parse_pdb_text(text: str, pdb_id: str = "") -> Structure:
    return parse_pdb(text, pdb_id=pdb_id)
