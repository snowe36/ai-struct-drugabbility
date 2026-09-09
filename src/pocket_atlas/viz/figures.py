from __future__ import annotations

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.patches import FancyBboxPatch

from pocket_atlas.paths import FIGURES, ensure_dirs
from pocket_atlas.pipeline import Campaign
from pocket_atlas.viz.palette import ACCENT, CARD_BG, CORAL, FACE, GRID, MUTED, TEAL, TEXT


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
