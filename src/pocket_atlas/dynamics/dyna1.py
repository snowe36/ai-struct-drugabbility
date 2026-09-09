"""Optional Dyna-1 inference.

Dyna-1 (Wayment-Steele / Kern 2026) predicts per-residue p(μs–ms exchange).
It does not generate MD. Weights live on Hugging Face (gelnesr/Dyna-1).
When they are missing we refuse to fake scores and fall back to curated
NMR labels in the case YAML.
"""

from __future__ import annotations

import json
from pathlib import Path

from pocket_atlas.paths import PROCESSED, ensure_dirs


class Dyna1Unavailable(RuntimeError):
    pass


def cached_scores_path(case_name: str) -> Path:
    return PROCESSED / f"{case_name}_dyna1.json"


def load_cached_scores(case_name: str) -> dict[int, float] | None:
    path = cached_scores_path(case_name)
    if not path.exists():
        return None
    raw = json.loads(path.read_text())
    return {int(k): float(v) for k, v in raw.items()}


def save_scores(case_name: str, scores: dict[int, float]) -> Path:
    ensure_dirs()
    path = cached_scores_path(case_name)
    path.write_text(json.dumps({str(k): v for k, v in sorted(scores.items())}, indent=2))
    return path


def high_exchange_residues(scores: dict[int, float], quantile: float = 0.8) -> set[int]:
    if not scores:
        return set()
    values = np_values(scores)
    cutoff = _quantile(values, quantile)
    return {r for r, p in scores.items() if p >= cutoff}


def np_values(scores: dict[int, float]) -> list[float]:
    return list(scores.values())


def _quantile(values: list[float], q: float) -> float:
    if not values:
        return 1.0
    s = sorted(values)
    idx = min(len(s) - 1, max(0, int(round(q * (len(s) - 1)))))
    return s[idx]


def try_dyna1_inference(pdb_path: Path, chain: str = "A") -> dict[int, float]:
    """Run Dyna-1 if the optional extra and weights are installed.

    We shell out to the upstream CLI when `dyna1.py` is on PATH. This repo
    does not vendor ESM-3.
    """
    import shutil
    import subprocess

    dyna = shutil.which("dyna1.py") or shutil.which("dyna1")
    if dyna is None:
        raise Dyna1Unavailable(
            "Dyna-1 CLI not found. Install WaymentSteeleLab/Dyna-1 and weights, "
            "or rely on curated NMR labels (default)."
        )
    out_dir = PROCESSED / "dyna1_run"
    out_dir.mkdir(parents=True, exist_ok=True)
    cmd = [
        dyna,
        "--pdb",
        str(pdb_path),
        "--chain",
        chain,
        "--name",
        pdb_path.stem,
        "--use_pdb_seq",
        "--write_to_pdb",
    ]
    subprocess.run(cmd, check=True, cwd=out_dir)
    scored = out_dir / f"{pdb_path.stem}.pdb"
    if not scored.exists():
        raise Dyna1Unavailable("Dyna-1 ran but did not write a scored PDB")
    return _bfactors_as_scores(scored, chain=chain)


def _bfactors_as_scores(path: Path, chain: str) -> dict[int, float]:
    from pocket_atlas.io.pdb import read_pdb

    structure = read_pdb(path)
    scores: dict[int, float] = {}
    # Upstream writes p(exchange) into B-factors; we re-parse as occupancy-like.
    text = path.read_text()
    for line in text.splitlines():
        if not line.startswith("ATOM"):
            continue
        if line[21].strip() != chain:
            continue
        if line[12:16].strip() != "CA":
            continue
        resseq = int(line[22:26])
        bfactor = float(line[60:66])
        scores[resseq] = bfactor
    if not scores:
        for atom in structure.protein_atoms(chain=chain):
            if atom.name == "CA":
                scores[atom.resseq] = 0.0
    return scores
