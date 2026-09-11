"""Optional Dyna-1 inference.

Upstream CLI (WaymentSteeleLab/Dyna-1 `dyna1.py`):
  python dyna1.py --pdb PATH --chain A --name NAME --use_pdb_seq --write_to_pdb --save_dir DIR

Writes `{name}-Dyna1.csv` (position, residue, p_exchange) and `{name}-Dyna1.pdb`.
Must run with cwd = the Dyna-1 repo so configs/esm3.yml and model/weights resolve.
Set DYNA1_ROOT if the CLI is not already launched from that tree.
"""

from __future__ import annotations

import json
import os
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


def remap_index_scores(scores: dict[int, float], pdb_path: Path, chain: str = "A") -> dict[int, float]:
    """Map 1..N Dyna-1 positions onto PDB resseq (3FKE is 218–340, not 1–123)."""
    from pocket_atlas.io.pdb import read_pdb

    resseqs = sorted(read_pdb(pdb_path).residue_ca(chain=chain))
    ordered = [scores[i] for i in sorted(scores)]
    n = min(len(resseqs), len(ordered))
    return {resseqs[i]: ordered[i] for i in range(n)}


def high_exchange_residues(scores: dict[int, float], quantile: float = 0.8) -> set[int]:
    if not scores:
        return set()
    values = list(scores.values())
    cutoff = _quantile(values, quantile)
    return {r for r, p in scores.items() if p >= cutoff}


def _quantile(values: list[float], q: float) -> float:
    if not values:
        return 1.0
    s = sorted(values)
    idx = min(len(s) - 1, max(0, int(round(q * (len(s) - 1)))))
    return s[idx]


def try_dyna1_inference(pdb_path: Path, chain: str = "A", name: str | None = None) -> dict[int, float]:
    import shutil
    import subprocess

    dyna = shutil.which("dyna1.py") or shutil.which("dyna1")
    if dyna is None:
        raise Dyna1Unavailable(
            "Dyna-1 CLI not found. Install WaymentSteeleLab/Dyna-1 and gelnesr/Dyna-1 "
            "weights, or use --prior literature / relaxdb."
        )
    name = name or pdb_path.stem
    out_dir = PROCESSED / "dyna1_run"
    out_dir.mkdir(parents=True, exist_ok=True)
    cmd = [
        dyna,
        "--pdb",
        str(Path(pdb_path).resolve()),
        "--chain",
        chain,
        "--name",
        name,
        "--use_pdb_seq",
        "--write_to_pdb",
        "--save_dir",
        str(out_dir),
    ]
    cwd = os.environ.get("DYNA1_ROOT") or str(out_dir)
    subprocess.run(cmd, check=True, cwd=cwd)
    csv_path = out_dir / f"{name}-Dyna1.csv"
    if csv_path.exists():
        return _scores_from_csv(csv_path)
    scored = out_dir / f"{name}-Dyna1.pdb"
    if scored.exists():
        return _bfactors_as_scores(scored, chain=chain)
    raise Dyna1Unavailable(f"Dyna-1 ran but did not write {name}-Dyna1.csv or .pdb")


def _scores_from_csv(path: Path) -> dict[int, float]:
    scores: dict[int, float] = {}
    for i, line in enumerate(path.read_text().splitlines()):
        if i == 0 and "p_exchange" in line:
            continue
        parts = [p.strip() for p in line.split(",")]
        if len(parts) < 3:
            continue
        scores[int(float(parts[0]))] = float(parts[2])
    if not scores:
        raise Dyna1Unavailable(f"no p_exchange rows in {path}")
    return scores


def _bfactors_as_scores(path: Path, chain: str) -> dict[int, float]:
    scores: dict[int, float] = {}
    for line in path.read_text().splitlines():
        if not line.startswith("ATOM"):
            continue
        if line[21].strip() != chain:
            continue
        if line[12:16].strip() != "CA":
            continue
        scores[int(line[22:26])] = float(line[60:66])
    return scores
