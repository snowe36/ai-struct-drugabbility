from __future__ import annotations

from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
DATA = ROOT / "data"
RAW = DATA / "raw"
PROCESSED = DATA / "processed"
REPORTS = ROOT / "out"
FIGURES = REPORTS / "figures"
DEMO_OUT = ROOT / "demo" / "outputs"
CASES = Path(__file__).resolve().parent / "cases"


def ensure_dirs() -> None:
    for path in (RAW, PROCESSED, FIGURES, DEMO_OUT):
        path.mkdir(parents=True, exist_ok=True)
