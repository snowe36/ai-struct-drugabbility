from __future__ import annotations

import argparse
import json
from pathlib import Path

from pocket_atlas.cases import list_cases, load_case
from pocket_atlas.dynamics.dyna1 import Dyna1Unavailable, save_scores, try_dyna1_inference
from pocket_atlas.dynamics.labels import PRIORS
from pocket_atlas.io.rcsb import fetch_pdb
from pocket_atlas.paths import ensure_dirs
from pocket_atlas.pipeline import run_campaign, write_campaign_json
from pocket_atlas.report import write_markdown
from pocket_atlas.sample import TrajectorySkipped, analyze_vp35, extract_vp35, fetch_vp35
from pocket_atlas.viz.figures import write_all_figures


def _case_arg(parser: argparse.ArgumentParser) -> None:
    parser.add_argument("--case", default="tem1_horn", choices=list_cases())


def _prior_arg(parser: argparse.ArgumentParser) -> None:
    parser.add_argument(
        "--prior",
        default="literature",
        choices=PRIORS,
        help="literature = YAML (CI/demo). relaxdb = RelaxDB-CPMG. dyna1 = cached scores.",
    )


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
    _prior_arg(parser)
    parser.add_argument("--ligand-mode", default="exclude", choices=("exclude", "include"))
    parser.add_argument("--no-both-modes", action="store_true")
    args = parser.parse_args(argv)
    camp = run_campaign(
        args.case,
        prefer_dyna1=args.prior == "dyna1",
        prior=args.prior,
        holo_ligand_mode=args.ligand_mode,
        both_holo_modes=not args.no_both_modes,
    )
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
    case = load_case(args.case)
    try:
        scores = try_dyna1_inference(Path(args.pdb), chain=case.chain, name=case.name)
    except Dyna1Unavailable as exc:
        print(exc)
        return
    path = save_scores(args.case, scores)
    print(f"wrote {len(scores)} residue scores -> {path}")


def overlap_main(argv: list[str] | None = None) -> None:
    parser = argparse.ArgumentParser(description="NMR/Dyna-1 overlap with cryptic lining")
    parser.add_argument("--cases", nargs="+", default=["tem1_horn", "kras_switch2"])
    _prior_arg(parser)
    parser.add_argument("--prefer-dyna1", action="store_true")
    args = parser.parse_args(argv)
    camps = [
        run_campaign(name, prefer_dyna1=args.prefer_dyna1, prior=args.prior)
        for name in args.cases
    ]
    md = write_markdown(camps)
    print(md.read_text())


def fetch_md_main(argv: list[str] | None = None) -> None:
    parser = argparse.ArgumentParser(description="Download VP35 Zenodo trajectory (multi-GB)")
    parser.add_argument("--which", default=None)
    parser.add_argument("--yes", action="store_true", help="Confirm the large download")
    parser.add_argument("--extract", action="store_true", help="Extract + stride after download")
    parser.add_argument("--analyze", action="store_true", help="Score CA 225–295 occupancy")
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
    if args.extract:
        manifest = extract_vp35(path)
        print(json.dumps(manifest, indent=2))
    if args.analyze:
        result = analyze_vp35()
        print(json.dumps({k: v for k, v in result.items() if k != "path"}, indent=2))


def md_analyze_main(argv: list[str] | None = None) -> None:
    parser = argparse.ArgumentParser(description="Extract/analyze VP35 FAH (not in demo/CI)")
    parser.add_argument("--archive", default=None, help="Path to downloaded tarball")
    parser.add_argument("--extract", action="store_true")
    parser.add_argument("--analyze", action="store_true")
    args = parser.parse_args(argv)
    if not args.extract and not args.analyze:
        parser.error("pass --extract and/or --analyze")
    try:
        if args.extract:
            if not args.archive:
                parser.error("--extract needs --archive")
            print(json.dumps(extract_vp35(Path(args.archive)), indent=2))
        if args.analyze:
            print(json.dumps(analyze_vp35(), indent=2))
    except TrajectorySkipped as exc:
        print(exc)


def fetch_nmr_main(argv: list[str] | None = None) -> None:
    parser = argparse.ArgumentParser(
        description="Refresh RelaxDB-CPMG labels from the Dyna-1 repo zip (optional)"
    )
    parser.add_argument("--url", default=(
        "https://raw.githubusercontent.com/WaymentSteeleLab/Dyna-1/main/"
        "data/RelaxDB_datasets/RelaxDB_CPMG_4jun2026.json.zip"
    ))
    args = parser.parse_args(argv)
    from pocket_atlas.dynamics.labels import RESOURCES, load_relaxdb_bundle

    print(f"Vendored labels live at {RESOURCES}")
    bundle = load_relaxdb_bundle()
    print(f"source: {bundle.get('source')}")
    print(f"note: {bundle.get('note')}")
    print("CI uses --prior literature; this file is the RelaxDB-CPMG fallback.")
    print(f"(re-fetch URL, not run: {args.url})")


def report_main(argv: list[str] | None = None) -> None:
    parser = argparse.ArgumentParser(description="Write markdown reports")
    parser.add_argument("--cases", nargs="+", default=["tem1_horn", "kras_switch2"])
    _prior_arg(parser)
    parser.add_argument("--prefer-dyna1", action="store_true")
    args = parser.parse_args(argv)
    camps = [
        run_campaign(name, prefer_dyna1=args.prefer_dyna1, prior=args.prior)
        for name in args.cases
    ]
    for camp in camps:
        write_campaign_json(camp)
    path = write_markdown(camps)
    print(f"wrote {path}")


def figures_main(argv: list[str] | None = None) -> None:
    parser = argparse.ArgumentParser(description="Write matplotlib figures")
    parser.add_argument("--cases", nargs="+", default=["tem1_horn", "kras_switch2"])
    _prior_arg(parser)
    parser.add_argument("--prefer-dyna1", action="store_true")
    args = parser.parse_args(argv)
    camps = [
        run_campaign(name, prefer_dyna1=args.prefer_dyna1, prior=args.prior)
        for name in args.cases
    ]
    paths = write_all_figures(camps)
    for path in paths:
        print(path)


def dock_main(argv: list[str] | None = None) -> None:
    parser = argparse.ArgumentParser(description="GNINA residue-box docking (fail closed)")
    parser.add_argument("--cases", nargs="+", default=["tem1_horn", "kras_switch2"])
    _prior_arg(parser)
    parser.add_argument("--exhaustiveness", type=int, default=8)
    args = parser.parse_args(argv)
    from pocket_atlas.dock import DockUnavailable, dock_case, write_dock_report

    jobs = []
    try:
        for name in args.cases:
            jobs.extend(dock_case(name, prior=args.prior, exhaustiveness=args.exhaustiveness))
    except DockUnavailable as exc:
        print(exc)
        return
    path = write_dock_report(jobs)
    print(path.read_text())


def strain_main(argv: list[str] | None = None) -> None:
    parser = argparse.ArgumentParser(description="Ligand-only GFN2-xTB strain (fail closed)")
    parser.add_argument("--cases", nargs="+", default=["tem1_horn", "kras_switch2"])
    args = parser.parse_args(argv)
    from pocket_atlas.chem.xtb_strain import XtbUnavailable, strain_case, write_strain_report

    payloads = []
    try:
        for name in args.cases:
            payloads.append(strain_case(name))
    except XtbUnavailable as exc:
        print(exc)
        return
    path = write_strain_report(payloads)
    print(path.read_text())


def prepare_crystal_main(argv: list[str] | None = None) -> None:
    parser = argparse.ArgumentParser(description="Prepare a crystal-only case (VP35 3FKE)")
    _case_arg(parser)
    args = parser.parse_args(argv)
    from pocket_atlas.pipeline import prepare_crystal

    path = prepare_crystal(args.case)
    print(path)


def rank_main(argv: list[str] | None = None) -> None:
    parser = argparse.ArgumentParser(
        description="Rank apo cavities without YAML site linings, then eval ranks"
    )
    parser.add_argument("--case", default=None, choices=list_cases())
    parser.add_argument(
        "--cases",
        nargs="+",
        default=None,
        choices=list_cases(),
    )
    _prior_arg(parser)
    parser.add_argument(
        "--open-frame",
        action="store_true",
        help="VP35: also rank a high-CV FAH frame if the strided xtc exists",
    )
    args = parser.parse_args(argv)
    names = args.cases or ([args.case] if args.case else ["tem1_horn", "kras_switch2"])
    from pocket_atlas.rank import rank_case
    from pocket_atlas.rank.eval import (
        eval_labeled_sites,
        write_combined_rank_report,
        write_rank_report,
    )
    from pocket_atlas.rank.miner import MinerUnavailable, try_miner_scores
    from pocket_atlas.viz.figures import fig_rank_bars

    combined = []
    miner_skip_printed = False
    for name in names:
        miner_scores = None
        try:
            miner_scores = try_miner_scores(name)
        except MinerUnavailable as exc:
            if not miner_skip_printed:
                print(f"PocketMiner skipped: {exc}")
                miner_skip_printed = True
        ranking = rank_case(name, prior=args.prior, miner_scores=miner_scores)
        sites = eval_labeled_sites(ranking)
        write_rank_report(ranking, sites)
        combined.append((ranking, sites))
        print(json.dumps({"ranking": ranking.as_dict(), "sites": [s.as_dict() for s in sites]}, indent=2))

    if "vp35_iid" in names and args.open_frame:
        from pocket_atlas.rank.appearance import FrameUnavailable, rank_open_frame

        try:
            opened, extra = rank_open_frame(prior=args.prior)
            open_sites = eval_labeled_sites(opened)
            write_rank_report(opened, open_sites)
            combined.append((opened, open_sites))
            print(json.dumps({
                "open": opened.as_dict(),
                "sites": [s.as_dict() for s in open_sites],
                "extra": extra,
            }, indent=2))
        except FrameUnavailable as exc:
            print(f"VP35 open-frame skipped: {exc}")

    path = write_combined_rank_report(combined)
    fig = fig_rank_bars(combined)
    print(f"wrote {path}")
    print(f"wrote {fig}")


def plane_main(argv: list[str] | None = None) -> None:
    parser = argparse.ArgumentParser(
        description="Map experimental CPMG/literature exchange onto apo clearance maxima"
    )
    parser.parse_args(argv)
    from pocket_atlas.rank.plane import run_exchange_plane, write_plane_report
    from pocket_atlas.viz.figures import fig_exchange_plane

    rows = run_exchange_plane()
    report = write_plane_report(rows)
    fig = fig_exchange_plane(rows)
    for ranking in rows:
        top = max(ranking.scored, key=lambda s: s.enrichment) if ranking.scored else None
        extra = ranking.extra
        print(
            f"{ranking.case:12s} {ranking.pdb_id}  pockets={len(ranking.scored):2d}  "
            f"prior={ranking.n_prior:3d}  mapped={extra.get('mapped_frac', '')}  "
            f"top_enr={None if top is None else round(top.enrichment, 2)}"
        )
    print(f"wrote {report}")
    print(f"wrote {fig}")
