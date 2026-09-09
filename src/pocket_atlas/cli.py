from __future__ import annotations

import argparse
import json

from pocket_atlas.cases import list_cases, load_case
from pocket_atlas.dynamics.dyna1 import Dyna1Unavailable, save_scores, try_dyna1_inference
from pocket_atlas.io.rcsb import fetch_pdb
from pocket_atlas.paths import ensure_dirs
from pocket_atlas.pipeline import run_campaign, write_campaign_json
from pocket_atlas.report import write_markdown
from pocket_atlas.sample import TrajectorySkipped, fetch_vp35
from pocket_atlas.viz.figures import write_all_figures


def _case_arg(parser: argparse.ArgumentParser) -> None:
    parser.add_argument("--case", default="tem1_horn", choices=list_cases())


def fetch_main(argv: list[str] | None = None) -> None:
    parser = argparse.ArgumentParser(description="Download case PDBs from RCSB")
    _case_arg(parser)
    args = parser.parse_args(argv)
    case = load_case(args.case)
    ensure_dirs()
    for spec in case.raw.get("structures", {}).values():
        pdb_id = spec.get("pdb_id")
        if pdb_id:
            path = fetch_pdb(pdb_id)
            print(f"fetched {pdb_id} -> {path}")


def prepare_main(argv: list[str] | None = None) -> None:
    parser = argparse.ArgumentParser(description="Run apo/holo prepare + detect")
    _case_arg(parser)
    args = parser.parse_args(argv)
    camp = run_campaign(args.case)
    path = write_campaign_json(camp)
    print(json.dumps(camp.as_dict(), indent=2))
    print(f"wrote {path}")


def detect_main(argv: list[str] | None = None) -> None:
    prepare_main(argv)


def dyna_main(argv: list[str] | None = None) -> None:
    parser = argparse.ArgumentParser(description="Optional Dyna-1 inference")
    _case_arg(parser)
    parser.add_argument("--pdb", required=True)
    args = parser.parse_args(argv)
    from pathlib import Path

    case = load_case(args.case)
    try:
        scores = try_dyna1_inference(Path(args.pdb), chain=case.chain)
    except Dyna1Unavailable as exc:
        print(exc)
        return
    path = save_scores(args.case, scores)
    print(f"wrote {len(scores)} residue scores -> {path}")


def overlap_main(argv: list[str] | None = None) -> None:
    parser = argparse.ArgumentParser(description="NMR/Dyna-1 overlap with cryptic lining")
    parser.add_argument("--cases", nargs="+", default=["tem1_horn", "kras_switch2"])
    args = parser.parse_args(argv)
    camps = [run_campaign(name) for name in args.cases]
    md = write_markdown(camps)
    print(md.read_text())


def fetch_md_main(argv: list[str] | None = None) -> None:
    parser = argparse.ArgumentParser(description="Download VP35 Zenodo trajectory (multi-GB)")
    parser.add_argument("--which", default=None)
    parser.add_argument("--yes", action="store_true", help="Confirm the large download")
    args = parser.parse_args(argv)
    if not args.yes:
        print("Refusing to download multi-GB FAH archives without --yes")
        print("See src/pocket_atlas/cases/vp35_iid.yaml for Zenodo 15854842")
        return
    try:
        path = fetch_vp35(which=args.which)
    except TrajectorySkipped as exc:
        print(exc)
        return
    print(f"downloaded {path} ({path.stat().st_size / 1e9:.2f} GB)")


def report_main(argv: list[str] | None = None) -> None:
    parser = argparse.ArgumentParser(description="Write markdown reports")
    parser.add_argument("--cases", nargs="+", default=["tem1_horn", "kras_switch2"])
    args = parser.parse_args(argv)
    camps = [run_campaign(name) for name in args.cases]
    for camp in camps:
        write_campaign_json(camp)
    path = write_markdown(camps)
    print(f"wrote {path}")


def figures_main(argv: list[str] | None = None) -> None:
    parser = argparse.ArgumentParser(description="Write matplotlib figures")
    parser.add_argument("--cases", nargs="+", default=["tem1_horn", "kras_switch2"])
    args = parser.parse_args(argv)
    camps = [run_campaign(name) for name in args.cases]
    paths = write_all_figures(camps)
    for path in paths:
        print(path)
