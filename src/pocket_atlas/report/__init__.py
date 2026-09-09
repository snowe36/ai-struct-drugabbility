from __future__ import annotations

from pathlib import Path

from pocket_atlas.paths import REPORTS, ensure_dirs
from pocket_atlas.pipeline import Campaign


def write_markdown(campaigns: list[Campaign], path: Path | None = None) -> Path:
    ensure_dirs()
    path = path or (REPORTS / "overlap.md")
    chunks = [
        "# Pocket Atlas — NMR dynamics vs cryptic sites",
        "",
        "Dyna-1 predicts *where* micro–millisecond motion occurs. It does not",
        "generate molecular dynamics. 100 ns public atlases are the wrong",
        "ensemble for cryptic opening; VP35 FAH/FAST is the long-MD arm.",
        "",
    ]
    for camp in campaigns:
        ov = camp.overlap
        apo_vol = camp.apo.site_pocket.volume if camp.apo.site_pocket else 0.0
        holo_vol = camp.holo.site_pocket.volume if camp.holo.site_pocket else 0.0
        apo_clr = camp.extra.get("apo_clearance", 0.0)
        holo_clr = camp.extra.get("holo_clearance", 0.0)
        chunks += [
            f"## {camp.case.raw.get('title', camp.case.name)}",
            "",
            f"- Apo `{camp.apo.pdb_id}` site recovered: **{camp.cryptic_in_apo}** "
            f"(clearance {apo_clr:.2f} Å; volume {apo_vol:.0f} Å³ is secondary)",
            f"- Holo `{camp.holo.pdb_id}` site recovered: **{camp.cryptic_in_holo}** "
            f"(clearance {holo_clr:.2f} Å; volume {holo_vol:.0f} Å³ is secondary)",
            f"- NMR ∩ cryptic lining: {ov.cryptic_and_nmr}/{ov.n_cryptic} (enrichment {ov.cryptic_enrichment:.2f})",
            f"- NMR ∩ control site: {ov.control_and_nmr}/{ov.n_control} (enrichment {ov.control_enrichment:.2f})",
            f"- Odds ratio (cryptic vs control): {ov.odds_ratio:.2f}",
            f"- Prior: `{ov.source}` (Dyna-1 vs literature when both exist)",
            "",
        ]
        if camp.case.name == "tem1_horn":
            chunks += [
                "TEM-1 horn lining is mostly buried hydrophobic core, not the Ω-loop.",
                "Savard/Gagné μs–ms exchange is reported at the Ω-loop and active-site",
                "vicinity. If cryptic enrichment is low, that is a result: NMR dynamics",
                "are a prior for *motion*, not a pocket oracle.",
                "",
            ]
        if camp.case.name == "kras_switch2":
            chunks += [
                "KRAS switch-I/II carry most published μs–ms NMR signal and also line",
                "the sotorasib site. Enrichment here is the expected positive control",
                "for the same question.",
                "",
            ]
    chunks += [
        "## Limitations",
        "",
        "- `--prior literature` is the curated YAML (CI/demo, offline).",
        "- `--prior relaxdb` uses RelaxDB-CPMG. KRAS is in that set; TEM-1 is not "
        "(BLAC_CPMG is Mtb BlaC P9WKD3, not TEM-1).",
        "- Dyna-1 weights are optional; captions must say Dyna-1 vs literature.",
        "- Headline geometry is seed clearance, not matched-site volume.",
        "- VP35 trajectories are not downloaded by the demo (multi-GB Zenodo archives).",
        "",
    ]
    path.write_text("\n".join(chunks))
    return path
