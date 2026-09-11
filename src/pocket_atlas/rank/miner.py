"""PocketMiner lining-mean baseline. Fail closed; do not train a GNN here."""

from __future__ import annotations

import json
import os
import shutil
from pathlib import Path

from pocket_atlas.paths import PROCESSED


class MinerUnavailable(RuntimeError):
    pass


def find_miner() -> str:
    path = os.environ.get("POCKETMINER") or shutil.which("pocketminer")
    if path is None:
        raise MinerUnavailable(
            "PocketMiner not on PATH. Set POCKETMINER or cache "
            "data/processed/{case}_miner.json. This repo does not train a GNN."
        )
    return path


def scores_from_bfactor_pdb(path: Path, chain: str = "A") -> dict[int, float]:
    scores: dict[int, float] = {}
    for line in Path(path).read_text().splitlines():
        if not line.startswith("ATOM"):
            continue
        if line[21].strip() != chain:
            continue
        if line[12:16].strip() != "CA":
            continue
        if len(line) < 66:
            continue
        scores[int(line[22:26])] = float(line[60:66])
    if not scores:
        raise MinerUnavailable(f"no CA B-factors in {path}")
    return scores


def load_cached_miner(case_name: str) -> dict[int, float] | None:
    path = PROCESSED / f"{case_name}_miner.json"
    if not path.exists():
        return None
    raw = json.loads(path.read_text())
    return {int(k): float(v) for k, v in raw.items()}


def try_miner_scores(case_name: str, pdb_path: Path | None = None) -> dict[int, float]:
    cached = load_cached_miner(case_name)
    if cached:
        return cached
    if pdb_path is not None and Path(pdb_path).exists():
        return scores_from_bfactor_pdb(Path(pdb_path))
    find_miner()
    raise MinerUnavailable(
        "PocketMiner CLI is present but this wrapper does not shell out. "
        "Cache residue probabilities as data/processed/{case}_miner.json."
    )
