from __future__ import annotations

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.patches import FancyBboxPatch

from pocket_atlas.paths import FIGURES, ensure_dirs
from pocket_atlas.pipeline import Campaign
from pocket_atlas.viz.palette import (
    ACCENT,
    CARD_BG,
    CORAL,
    FACE,
    GRID,
    MUSTARD,
    MUTED,
    SAGE,
    TEAL,
    TEXT,
)


def _style(ax) -> None:
    ax.set_facecolor(FACE)
    for spine in ax.spines.values():
        spine.set_color(GRID)
    ax.tick_params(colors=MUTED)
    ax.xaxis.label.set_color(TEXT)
    ax.yaxis.label.set_color(TEXT)
    ax.title.set_color(TEXT)


def _prior_label(source: str) -> str:
    if source == "dyna1":
        return "Dyna-1"
    if source == "relaxdb_cpmg":
        return "RelaxDB-CPMG"
    return "literature NMR"


def _case_title(camp: Campaign) -> str:
    raw = camp.case.raw.get("title") or camp.case.name
    if camp.case.name == "tem1_horn":
        return "TEM-1 horn"
    if camp.case.name == "kras_switch2":
        return "KRAS switch-II"
    return raw


def fig_overlap_enrichment(campaigns: list[Campaign], path: Path | None = None) -> Path:
    """Hero: one panel per case so a 0.00 bar and an 8× control are both readable."""
    ensure_dirs()
    path = path or (FIGURES / "fig1_overlap_enrichment.png")
    n = max(len(campaigns), 1)
    fig, axes = plt.subplots(1, n, figsize=(3.54 * n, 3.9), dpi=200, sharey=True)
    fig.patch.set_facecolor(FACE)
    if n == 1:
        axes = [axes]

    ymax = 1.4
    for camp in campaigns:
        ymax = max(ymax, camp.overlap.cryptic_enrichment, camp.overlap.control_enrichment)
    ymax = ymax * 1.28

    sources = []
    for ax, camp in zip(axes, campaigns, strict=True):
        ov = camp.overlap
        sources.append(ov.source)
        vals = [ov.cryptic_enrichment, ov.control_enrichment]
        counts = [
            f"{ov.cryptic_and_nmr}/{ov.n_cryptic}",
            f"{ov.control_and_nmr}/{ov.n_control}",
        ]
        names = ["Cryptic lining", "Catalytic / nucleotide"]
        colors = [TEAL, ACCENT]
        ax.bar(range(2), vals, color=colors, width=0.62, zorder=2)
        ax.axhline(1.0, color=MUTED, ls="--", lw=1, zorder=1)
        trans = ax.get_xaxis_transform()  # data x, axes y
        for i, (val, count) in enumerate(zip(vals, counts, strict=True)):
            if val < 0.35:
                ax.text(
                    i,
                    -0.14,
                    f"{val:.2f} ({count})",
                    transform=trans,
                    ha="center",
                    va="top",
                    color=TEXT,
                    fontsize=8,
                    clip_on=False,
                )
            else:
                pad = 0.38 if abs(val - 1.0) < 0.25 else 0.2
                ax.text(
                    i,
                    val + pad,
                    f"{val:.2f} ({count})",
                    ha="center",
                    va="bottom",
                    color=TEXT,
                    fontsize=8,
                    zorder=3,
                )
        ax.set_xticks(range(2))
        ax.set_xticklabels(names)
        ax.set_ylim(0, ymax)
        ax.set_title(_case_title(camp), color=TEXT, fontsize=11)
        _style(ax)
        ax.tick_params(axis="x", labelsize=8, pad=10)

    axes[0].set_ylabel("Enrichment vs protein background")
    fig.suptitle("Do NMR-timescale residues sit on cryptic sites?", color=TEXT, fontsize=12, y=1.02)
    prior = _prior_label(sources[0] if sources else "")
    extras = ""
    if len(set(sources)) > 1:
        extras = "  ·  " + "; ".join(
            f"{_case_title(c)}: {_prior_label(c.overlap.source)}" for c in campaigns
        )
    fig.text(
        0.5,
        -0.04,
        f"Prior: {prior}{extras}. Dashed line = no enrichment.",
        ha="center",
        color=MUTED,
        fontsize=8,
    )
    fig.tight_layout()
    fig.subplots_adjust(bottom=0.22)
    fig.savefig(path, bbox_inches="tight", facecolor=FACE)
    plt.close(fig)
    return path


def fig_dyna1_vs_literature(
    literature: list[Campaign],
    dyna1: list[Campaign],
    path: Path | None = None,
) -> Path:
    """Cryptic enrichment: published NMR vs Dyna-1 top-quintile p(exchange)."""
    ensure_dirs()
    path = path or (FIGURES / "fig5_dyna1_vs_literature.png")
    by_dyna = {c.case.name: c for c in dyna1}
    pairs = [(lit, by_dyna[lit.case.name]) for lit in literature if lit.case.name in by_dyna]
    n = max(len(pairs), 1)
    fig, axes = plt.subplots(1, n, figsize=(3.54 * n, 3.9), dpi=200, sharey=True)
    fig.patch.set_facecolor(FACE)
    if n == 1:
        axes = [axes]

    ymax = 1.4
    for lit, dyn in pairs:
        ymax = max(ymax, lit.overlap.cryptic_enrichment, dyn.overlap.cryptic_enrichment)
    ymax = ymax * 1.28

    for ax, (lit, dyn) in zip(axes, pairs, strict=True):
        vals = [lit.overlap.cryptic_enrichment, dyn.overlap.cryptic_enrichment]
        counts = [
            f"{lit.overlap.cryptic_and_nmr}/{lit.overlap.n_cryptic}",
            f"{dyn.overlap.cryptic_and_nmr}/{dyn.overlap.n_cryptic}",
        ]
        ax.bar(range(2), vals, color=[SAGE, MUSTARD], width=0.62, zorder=2)
        ax.axhline(1.0, color=MUTED, ls="--", lw=1, zorder=1)
        trans = ax.get_xaxis_transform()
        for i, (val, count) in enumerate(zip(vals, counts, strict=True)):
            if val < 0.35:
                ax.text(
                    i, -0.14, f"{val:.2f} ({count})",
                    transform=trans, ha="center", va="top", color=TEXT, fontsize=8, clip_on=False,
                )
            else:
                ax.text(
                    i, val + 0.12, f"{val:.2f} ({count})",
                    ha="center", va="bottom", color=TEXT, fontsize=8, zorder=3,
                )
        ax.set_xticks(range(2))
        ax.set_xticklabels(["Literature NMR", "Dyna-1"])
        ax.set_ylim(0, ymax)
        ax.set_title(_case_title(lit), color=TEXT, fontsize=11)
        _style(ax)
        ax.tick_params(axis="x", labelsize=8, pad=10)

    axes[0].set_ylabel("Cryptic enrichment vs protein background")
    fig.suptitle("Does predicted exchange recover the NMR prior?", color=TEXT, fontsize=12, y=1.02)
    fig.text(
        0.5, -0.04,
        "Dyna-1 = top quintile p(μs–ms exchange). Dashed line = no enrichment.",
        ha="center", color=MUTED, fontsize=8,
    )
    fig.tight_layout()
    fig.subplots_adjust(bottom=0.22)
    fig.savefig(path, bbox_inches="tight", facecolor=FACE)
    plt.close(fig)
    return path


def fig_apo_holo_volumes(campaign: Campaign, path: Path | None = None) -> Path:
    """Apo vs holo matched-site clearance. Volume sits under the state, not on the bar."""
    ensure_dirs()
    slug = "tem1" if campaign.case.name == "tem1_horn" else (
        "kras" if campaign.case.name == "kras_switch2" else campaign.case.name
    )
    path = path or (FIGURES / f"fig_{slug}_clearance.png")
    fig, ax = plt.subplots(figsize=(3.54, 3.5), dpi=200)
    fig.patch.set_facecolor(FACE)
    clears = []
    names = []
    vols = []
    for arm in (campaign.apo, campaign.holo):
        extra = arm.site_pocket.extra if arm.site_pocket else {}
        clears.append(float(extra.get("seed_clearance", 0.0)))
        vols.append(arm.site_pocket.volume if arm.site_pocket else 0.0)
        names.append(arm.tag)
    ax.bar(names, clears, color=[MUTED, TEAL], width=0.55)
    ax.set_ylabel("Matched-site seed clearance (Å)")
    ax.set_title(_case_title(campaign))
    ticks = [f"{name}\n{vol:.0f} Å³" for name, vol in zip(names, vols, strict=True)]
    ax.set_xticks(range(len(names)))
    ax.set_xticklabels(ticks)
    top = max(clears) if clears else 1.0
    ax.set_ylim(0, top * 1.2)
    _style(ax)
    fig.tight_layout()
    fig.savefig(path, bbox_inches="tight", facecolor=FACE)
    plt.close(fig)
    return path


def fig_cartoon_trace(campaign: Campaign, which: str = "holo", path: Path | None = None) -> Path:
    """Cα trace: cryptic lining vs NMR-dynamic vs the rest. Not PyMOL."""
    from pocket_atlas.io.rcsb import load_structure
    from pocket_atlas.prepare import prepare_structure

    ensure_dirs()
    path = path or (FIGURES / f"fig_{campaign.case.name}_{which}_trace.png")
    spec = campaign.case.raw["structures"][which]
    structure = prepare_structure(load_structure(spec["pdb_id"]), chain=campaign.case.chain)
    cas = structure.residue_ca(chain=campaign.case.chain)
    res = sorted(cas)
    xyz = np.vstack([cas[r] for r in res])
    # Simple PCA projection for a readable 2D trace.
    xyz = xyz - xyz.mean(axis=0)
    _, _, vt = np.linalg.svd(xyz, full_matrices=False)
    xy = xyz @ vt[:2].T

    cryptic = campaign.case.cryptic_residues
    nmr = set(campaign.extra.get("nmr_residues", [])) or set(campaign.case.nmr_residues)
    colors = []
    sizes = []
    for r in res:
        if r in cryptic and r in nmr:
            colors.append(CORAL)
            sizes.append(36)
        elif r in cryptic:
            colors.append(TEAL)
            sizes.append(28)
        elif r in nmr:
            colors.append(ACCENT)
            sizes.append(22)
        else:
            colors.append(GRID)
            sizes.append(8)

    fig, ax = plt.subplots(figsize=(3.54, 3.6), dpi=200)
    fig.patch.set_facecolor(FACE)
    ax.plot(xy[:, 0], xy[:, 1], color=GRID, lw=1.4, zorder=0)
    ax.scatter(xy[:, 0], xy[:, 1], c=colors, s=sizes, linewidths=0, zorder=1)
    ax.set_aspect("equal")
    ax.axis("off")
    ax.set_title(f"{campaign.case.name} {which}  ·  teal=cryptic  peach=NMR  coral=both")
    fig.tight_layout()
    fig.savefig(path, bbox_inches="tight", facecolor=FACE)
    plt.close(fig)
    return path


def fig_tractability_card(campaign: Campaign, path: Path | None = None) -> Path:
    ensure_dirs()
    path = path or (FIGURES / f"fig_{campaign.case.name}_card.png")
    ov = campaign.overlap
    apo_vol = campaign.apo.site_pocket.volume if campaign.apo.site_pocket else 0.0
    holo_vol = campaign.holo.site_pocket.volume if campaign.holo.site_pocket else 0.0
    apo_clr = campaign.extra.get("apo_clearance", 0.0)
    holo_clr = campaign.extra.get("holo_clearance", 0.0)
    lines = [
        campaign.case.raw.get("title", campaign.case.name),
        f"Apo site: {'yes' if campaign.cryptic_in_apo else 'no'}  "
        f"clearance {apo_clr:.2f} Å  ({apo_vol:.0f} Å³)",
        f"Holo site: {'yes' if campaign.cryptic_in_holo else 'no'}  "
        f"clearance {holo_clr:.2f} Å  ({holo_vol:.0f} Å³)",
        f"NMR ∩ cryptic lining: {ov.cryptic_and_nmr}/{ov.n_cryptic}",
        f"NMR ∩ control site:   {ov.control_and_nmr}/{ov.n_control}",
        f"Cryptic enrichment: {ov.cryptic_enrichment:.2f}   control: {ov.control_enrichment:.2f}",
        f"Prior: {ov.source}  ·  Dyna-1 vs literature",
    ]
    fig, ax = plt.subplots(figsize=(3.54, 3.2), dpi=200)
    fig.patch.set_facecolor(FACE)
    ax.axis("off")
    ax.add_patch(
        FancyBboxPatch(
            (0.02, 0.04), 0.96, 0.92,
            boxstyle="round,pad=0.02,rounding_size=0.04",
            facecolor=CARD_BG, edgecolor=GRID, transform=ax.transAxes, linewidth=1,
        )
    )
    ax.text(0.08, 0.88, lines[0], transform=ax.transAxes, color=TEXT, fontsize=10, fontweight="bold")
    ax.text(0.08, 0.72, "\n".join(lines[1:]), transform=ax.transAxes, color=MUTED, fontsize=8, family="monospace", va="top")
    fig.savefig(path, bbox_inches="tight", facecolor=FACE)
    plt.close(fig)
    return path


def write_all_figures(campaigns: list[Campaign]) -> list[Path]:
    ensure_dirs()
    paths = [fig_overlap_enrichment(campaigns, path=FIGURES / "fig1_overlap_enrichment.png")]
    clearance_names = {
        "tem1_horn": FIGURES / "fig2_tem1_clearance.png",
        "kras_switch2": FIGURES / "fig3_kras_clearance.png",
    }
    for camp in campaigns:
        dest = clearance_names.get(camp.case.name)
        paths.append(fig_apo_holo_volumes(camp, path=dest))
        paths.append(fig_cartoon_trace(camp, which="holo"))
        paths.append(fig_tractability_card(camp))
    return paths


def fig_dock_heatmap(jobs: list[dict], path: Path | None = None) -> Path:
    """Best-pose CNN affinity by case × box × apo/holo."""
    ensure_dirs()
    path = path or (FIGURES / "fig6_dock_cnn.png")
    cases = []
    for job in jobs:
        if job["case"] not in cases:
            cases.append(job["case"])
    boxes = ["cryptic", "nmr", "control"]
    box_color = {"cryptic": TEAL, "nmr": SAGE, "control": ACCENT}
    fig, axes = plt.subplots(1, max(len(cases), 1), figsize=(3.7 * max(len(cases), 1), 3.9), dpi=200, sharey=True)
    fig.patch.set_facecolor(FACE)
    if len(cases) <= 1:
        axes = [axes]
    for ax, case in zip(axes, cases, strict=False):
        x = np.arange(2)
        width = 0.25
        for k, box in enumerate(boxes):
            vals = []
            for rec in ("apo", "holo"):
                hit = next(
                    (
                        j["best"].get("cnn_affinity")
                        for j in jobs
                        if j["case"] == case and j["receptor"] == rec and j["box"] == box and j.get("best")
                    ),
                    None,
                )
                vals.append(np.nan if hit is None else float(hit))
            ax.bar(x + (k - 1) * width, vals, width=width, color=box_color[box], label=box, zorder=2)
        ax.set_xticks(x)
        ax.set_xticklabels(["apo", "holo"])
        title = "TEM-1 horn" if case == "tem1_horn" else (
            "KRAS switch-II" if case == "kras_switch2" else case
        )
        ax.set_title(title, color=TEXT, fontsize=11)
        _style(ax)
    axes[0].set_ylabel("GNINA CNN affinity")
    axes[-1].legend(frameon=False, fontsize=8)
    fig.suptitle("Does the NMR box recover the crystal ligand?", color=TEXT, fontsize=12, y=1.02)
    fig.text(
        0.5, -0.04,
        "Higher is better. Sotorasib is covalent to Cys12; this is a pose score, not ΔG.",
        ha="center", color=MUTED, fontsize=8,
    )
    fig.tight_layout()
    fig.savefig(path, bbox_inches="tight", facecolor=FACE)
    plt.close(fig)
    return path


def fig_xtb_strain(payloads: list[dict], path: Path | None = None) -> Path:
    ensure_dirs()
    path = path or (FIGURES / "fig7_xtb_strain.png")
    labels, vals, colors = [], [], []
    palette = [SAGE, TEAL, MUSTARD, ACCENT, CORAL]
    k = 0
    for payload in payloads:
        case = payload.get("case", "")
        short = "TEM-1" if "tem1" in case else ("KRAS" if "kras" in case else case)
        for name, row in payload.get("rows", {}).items():
            labels.append(f"{short}\n{name}")
            vals.append(float(row.get("strain_kcal") or 0.0))
            colors.append(palette[k % len(palette)])
            k += 1
    fig, ax = plt.subplots(figsize=(max(3.54, 0.9 * max(len(vals), 1)), 3.6), dpi=200)
    fig.patch.set_facecolor(FACE)
    ax.bar(range(len(vals)), vals, color=colors or [TEAL], width=0.65)
    ax.set_xticks(range(len(vals)))
    ax.set_xticklabels(labels or [""], fontsize=7)
    ax.set_ylabel("GFN2-xTB ligand strain (kcal/mol)")
    ax.set_title("Crystal vs docked ligand strain")
    _style(ax)
    fig.text(0.5, -0.06, "Ligand only. Not a binding free energy.", ha="center", color=MUTED, fontsize=8)
    fig.tight_layout()
    fig.savefig(path, bbox_inches="tight", facecolor=FACE)
    plt.close(fig)
    return path


def fig_vp35_occupancy(result: dict, path: Path | None = None) -> Path:
    ensure_dirs()
    path = path or (FIGURES / "fig8_vp35_occupancy.png")
    fig, ax = plt.subplots(figsize=(3.54, 3.6), dpi=200)
    fig.patch.set_facecolor(FACE)
    n_open = int(result.get("n_open", 0))
    n_closed = int(result.get("n_closed", 0))
    ax.bar(["closed", "open"], [n_closed, n_open], color=[MUTED, TEAL], width=0.55)
    ov = result.get("overlap") or {}
    ax.set_ylabel("Strided frames")
    ax.set_title("VP35 CA 225–295 occupancy")
    frac = result.get("open_fraction", 0.0)
    ax.text(
        0.5, -0.18,
        f"open {frac:.2f}  ·  lining ∩ prior {ov.get('cryptic_and_nmr', 0)}/{ov.get('n_cryptic', 0)}"
        f"  ({ov.get('source', 'none')})",
        transform=ax.transAxes, ha="center", color=MUTED, fontsize=8,
    )
    _style(ax)
    fig.tight_layout()
    fig.savefig(path, bbox_inches="tight", facecolor=FACE)
    plt.close(fig)
    return path


def fig_rank_bars(rows: list, path: Path | None = None) -> Path:
    """Geometry vs hybrid rank of labeled sites. Lower rank is better."""
    ensure_dirs()
    path = path or (FIGURES / "fig9_discovery_ranks.png")
    panels = [(ranking, sites) for ranking, sites in rows if sites]
    n = max(len(panels), 1)
    fig, axes = plt.subplots(1, n, figsize=(3.7 * n, 3.9), dpi=200)
    fig.patch.set_facecolor(FACE)
    if n == 1:
        axes = [axes]
    titles = {
        "tem1_horn": "TEM-1",
        "kras_switch2": "KRAS",
        "vp35_iid": "VP35",
    }
    site_labels = {
        "omega_loop": "Ω-loop",
        "switch2": "switch-II",
        "nucleotide": "nucleotide",
        "horn": "horn",
        "catalytic": "catalytic",
        "cryptic": "cryptic",
    }
    ymax = 1
    for _ranking, sites in panels:
        for row in sites:
            if row.found and row.geometry_rank:
                ymax = max(ymax, row.geometry_rank)
            if row.found and row.hybrid_rank:
                ymax = max(ymax, row.hybrid_rank)
    for ax, (ranking, sites) in zip(axes, panels, strict=False):
        labels = [site_labels.get(row.site, row.site) for row in sites]
        x = np.arange(len(labels))
        geo = [row.geometry_rank if row.found and row.geometry_rank else np.nan for row in sites]
        hyb = [row.hybrid_rank if row.found and row.hybrid_rank else np.nan for row in sites]
        width = 0.36
        ax.bar(x - width / 2, geo, width=width, color=TEAL, label="geometry", zorder=2)
        ax.bar(x + width / 2, hyb, width=width, color=SAGE, label="hybrid", zorder=2)
        ax.set_xticks(x)
        ax.set_xticklabels(labels, fontsize=8)
        ax.set_ylabel("Rank (lower is better)")
        ax.set_ylim(0.5, ymax + 0.8)
        ax.invert_yaxis()
        ax.set_title(titles.get(ranking.case, ranking.case), color=TEXT, fontsize=11)
        _style(ax)
    if panels:
        axes[-1].legend(frameon=False, fontsize=8)
    fig.suptitle("Apo cavity rank of labeled sites", color=TEXT, fontsize=12, y=1.02)
    fig.text(
        0.5, -0.05,
        "YAML linings used only after ranking. Hybrid = dscore × (1 + NMR enrichment).",
        ha="center", color=MUTED, fontsize=8,
    )
    fig.tight_layout()
    fig.savefig(path, bbox_inches="tight", facecolor=FACE)
    plt.close(fig)
    return path


def fig_exchange_plane(rankings: list, path: Path | None = None) -> Path:
    """Every apo cavity: lining exchange enrichment vs ligand-scale dscore."""
    from pocket_atlas.rank.plane import SHORT_LABELS

    ensure_dirs()
    path = path or (FIGURES / "fig10_exchange_clearance.png")
    n = len(rankings) + 1
    ncols = 4
    nrows = int(np.ceil(n / ncols))
    fig, axes = plt.subplots(nrows, ncols, figsize=(9.2, 2.35 * nrows), dpi=200, sharex=False, sharey=False)
    fig.patch.set_facecolor(FACE)
    axes = np.atleast_1d(axes).ravel()
    for ax, ranking in zip(axes, rankings, strict=False):
        xs = [s.enrichment for s in ranking.scored]
        ys = [s.pocket.dscore for s in ranking.scored]
        ax.axvline(1.0, color=MUTED, ls="--", lw=0.8, zorder=1)
        if ys:
            ax.axhline(float(np.median(ys)), color=GRID, ls=":", lw=0.8, zorder=1)
        ax.scatter(xs, ys, s=16, color=TEAL, alpha=0.8, zorder=2, edgecolors="none")
        if ranking.scored:
            top = max(ranking.scored, key=lambda s: s.enrichment)
            ax.scatter(
                [top.enrichment],
                [top.pocket.dscore],
                s=56,
                facecolors="none",
                edgecolors=CORAL,
                linewidths=1.3,
                zorder=3,
            )
        ax.set_title(SHORT_LABELS.get(ranking.case, ranking.case), color=TEXT, fontsize=9)
        _style(ax)
    key = axes[len(rankings)]
    key.set_xlim(0, 1)
    key.set_ylim(0, 1)
    key.set_xticks([])
    key.set_yticks([])
    key.text(0.25, 0.76, "static hole", ha="center", va="center", color=MUTED, fontsize=7, clip_on=False)
    key.text(0.75, 0.76, "dynamic site", ha="center", va="center", color=SAGE, fontsize=7, clip_on=False)
    key.text(0.25, 0.24, "noise", ha="center", va="center", color=MUTED, fontsize=7, clip_on=False)
    key.text(0.75, 0.24, "catalytic /\ngating", ha="center", va="center", color=CORAL, fontsize=7, clip_on=False)
    key.axvline(0.5, color=GRID, lw=0.8)
    key.axhline(0.5, color=GRID, lw=0.8)
    key.set_xlabel("enrichment →", fontsize=8)
    key.set_ylabel("dscore →", fontsize=8)
    key.set_title("quadrants", color=TEXT, fontsize=9)
    _style(key)
    for ax in axes[len(rankings) + 1 :]:
        ax.axis("off")
    fig.supxlabel("Exchange enrichment in cavity lining", color=TEXT, fontsize=10)
    fig.supylabel("Clearance-maxima dscore", color=TEXT, fontsize=10)
    fig.tight_layout()
    fig.savefig(path, bbox_inches="tight", facecolor=FACE)
    plt.close(fig)
    return path
