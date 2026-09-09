from __future__ import annotations

from pocket_atlas.paths import DEMO_OUT, FIGURES, REPORTS, ensure_dirs
from pocket_atlas.pipeline import run_campaign, write_campaign_json
from pocket_atlas.report import write_markdown
from pocket_atlas.viz.figures import write_all_figures


def run_demo() -> None:
    """Five-minute path: crystals + literature NMR labels. No FAH, no Dyna-1 weights."""
    ensure_dirs()
    cases = ["tem1_horn", "kras_switch2"]
    campaigns = []
    for name in cases:
        print(f"running {name} …")
        camp = run_campaign(name, prefer_dyna1=False)
        write_campaign_json(camp)
        campaigns.append(camp)
        print(
            f"  apo site={camp.cryptic_in_apo}  holo site={camp.cryptic_in_holo}  "
            f"cryptic enrichment={camp.overlap.cryptic_enrichment:.2f}  "
            f"control enrichment={camp.overlap.control_enrichment:.2f}"
        )
    md = write_markdown(campaigns, path=REPORTS / "overlap.md")
    figs = write_all_figures(campaigns)
    (DEMO_OUT / "README.txt").write_text(
        "Demo used literature NMR labels, not Dyna-1 weights, and did not download VP35 FAH.\n"
        f"Report: {md}\nFigures:\n" + "\n".join(str(p) for p in figs) + "\n"
    )
    print(f"wrote {md}")
    print(f"figures in {FIGURES}")


def run_demo_main(argv: list[str] | None = None) -> None:
    run_demo()
