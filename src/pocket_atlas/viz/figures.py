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


def fig_overlap_enrichment(campaigns: list[Campaign], path: Path | None = None) -> Path:
    """Hero figure: NMR/Dyna-1 enrichment in cryptic lining vs catalytic control."""
    ensure_dirs()
    path = path or (FIGURES / "fig4_overlap_enrichment.png")
    labels = []
    cryptic = []
    control = []
    sources = []
    for camp in campaigns:
        labels.append(camp.case.name.replace("_", " "))
        cryptic.append(camp.overlap.cryptic_enrichment)
        control.append(camp.overlap.control_enrichment)
        sources.append(camp.scores_source)

    fig, ax = plt.subplots(figsize=(7.08, 3.6), dpi=200)
    fig.patch.set_facecolor(FACE)
    x = np.arange(len(labels))
    width = 0.36
    ax.bar(x - width / 2, cryptic, width, color=TEAL, label="Cryptic lining")
    ax.bar(x + width / 2, control, width, color=ACCENT, label="Catalytic / nucleotide site")
    ax.axhline(1.0, color=MUTED, ls="--", lw=1, label="No enrichment")
    ax.set_xticks(x)
    ax.set_xticklabels(labels)
    ax.set_ylabel("Enrichment vs protein background")
    ax.set_title("Do NMR-timescale residues sit on cryptic sites?")
    ax.legend(frameon=False, loc="upper right")
    _style(ax)
    source = sources[0] if sources else ""
    ax.text(
        0.0,
        -0.22,
        f"Prior: {source}. Enrichment = (NMR ∩ set) / |set|, divided by NMR fraction in the chain.",
        transform=ax.transAxes,
        color=MUTED,
        fontsize=8,
    )
    fig.tight_layout()
    fig.savefig(path, bbox_inches="tight", facecolor=FACE)
    plt.close(fig)
    return path


def fig_apo_holo_volumes(campaign: Campaign, path: Path | None = None) -> Path:
    ensure_dirs()
    path = path or (FIGURES / f"fig_{campaign.case.name}_apo_holo.png")
    fig, ax = plt.subplots(figsize=(3.54, 3.4), dpi=200)
    fig.patch.set_facecolor(FACE)
    vols = []
    names = []
    for arm, _color in ((campaign.apo, MUTED), (campaign.holo, TEAL)):
        vol = arm.site_pocket.volume if arm.site_pocket else 0.0
        vols.append(vol)
        names.append(arm.tag)
    ax.bar(names, vols, color=[MUTED, TEAL])
    ax.set_ylabel("Matched-site volume (Å³)")
    ax.set_title(campaign.case.raw.get("title", campaign.case.name))
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
    lines = [
        campaign.case.raw.get("title", campaign.case.name),
        f"Apo site recovered: {'yes' if campaign.cryptic_in_apo else 'no'}  ({apo_vol:.0f} Å³)",
        f"Holo site recovered: {'yes' if campaign.cryptic_in_holo else 'no'}  ({holo_vol:.0f} Å³)",
        f"NMR ∩ cryptic lining: {ov.cryptic_and_nmr}/{ov.n_cryptic}",
        f"NMR ∩ control site:   {ov.control_and_nmr}/{ov.n_control}",
        f"Cryptic enrichment: {ov.cryptic_enrichment:.2f}   control: {ov.control_enrichment:.2f}",
        f"Prior: {ov.source}",
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
    paths = [fig_overlap_enrichment(campaigns)]
    for camp in campaigns:
        paths.append(fig_apo_holo_volumes(camp))
        paths.append(fig_cartoon_trace(camp, which="holo"))
        paths.append(fig_tractability_card(camp))
    return paths
