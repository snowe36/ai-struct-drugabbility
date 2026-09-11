"""Post-hoc labeled-site ranks. YAML linings are allowed here only."""

from __future__ import annotations

import json
from dataclasses import dataclass
from pathlib import Path

from pocket_atlas.cases import load_case
from pocket_atlas.paths import PROCESSED, REPORTS, ensure_dirs
from pocket_atlas.pockets import pocket_near_residues
from pocket_atlas.rank import ApoRanking, ScoredPocket

# Evaluation labels only. Ranking must not import this table.
LABELED_SITES = {
    "tem1_horn": (
        ("horn", "cryptic_site"),
        ("omega_loop", "omega_loop"),
        ("catalytic", "catalytic_site"),
    ),
    "kras_switch2": (
        ("switch2", "cryptic_site"),
        ("nucleotide", "nucleotide_site"),
    ),
    "vp35_iid": (
        ("cryptic", "cryptic_site"),
    ),
}


@dataclass
class SiteEval:
    site: str
    found: bool
    pocket_index: int | None
    geometry_rank: int | None
    prior_rank: int | None
    hybrid_rank: int | None
    miner_rank: int | None
    n_overlap: int
    dscore: float | None
    enrichment: float | None
    hybrid: float | None

    def as_dict(self) -> dict:
        return {
            "site": self.site,
            "found": self.found,
            "pocket_index": self.pocket_index,
            "geometry_rank": self.geometry_rank,
            "prior_rank": self.prior_rank,
            "hybrid_rank": self.hybrid_rank,
            "miner_rank": self.miner_rank,
            "n_overlap": self.n_overlap,
            "dscore": None if self.dscore is None else round(self.dscore, 3),
            "enrichment": None if self.enrichment is None else round(self.enrichment, 3),
            "hybrid": None if self.hybrid is None else round(self.hybrid, 3),
        }


def _scored_for_pocket(ranking: ApoRanking, pocket_index: int) -> ScoredPocket | None:
    for item in ranking.scored:
        if item.pocket.index == pocket_index:
            return item
    return None


def eval_labeled_sites(ranking: ApoRanking, min_overlap: int = 3) -> list[SiteEval]:
    case = load_case(ranking.case)
    rows: list[SiteEval] = []
    for site, key in LABELED_SITES.get(ranking.case, ()):
        residues = set(case.cryptic_residues if key == "cryptic_site" else case.residue_set(key))
        pocket = pocket_near_residues(ranking.pockets, residues, min_overlap=min_overlap)
        if pocket is None:
            rows.append(
                SiteEval(
                    site=site,
                    found=False,
                    pocket_index=None,
                    geometry_rank=None,
                    prior_rank=None,
                    hybrid_rank=None,
                    miner_rank=None,
                    n_overlap=0,
                    dscore=None,
                    enrichment=None,
                    hybrid=None,
                )
            )
            continue
        scored = _scored_for_pocket(ranking, pocket.index)
        n_overlap = len(set(pocket.lining_residues) & residues)
        rows.append(
            SiteEval(
                site=site,
                found=True,
                pocket_index=pocket.index,
                geometry_rank=None if scored is None else scored.geometry_rank,
                prior_rank=None if scored is None else scored.prior_rank,
                hybrid_rank=None if scored is None else scored.hybrid_rank,
                miner_rank=None if scored is None else scored.miner_rank,
                n_overlap=n_overlap,
                dscore=None if scored is None else scored.pocket.dscore,
                enrichment=None if scored is None else scored.enrichment,
                hybrid=None if scored is None else scored.hybrid,
            )
        )
    return rows


def write_rank_report(
    ranking: ApoRanking,
    sites: list[SiteEval],
    path: Path | None = None,
) -> Path:
    ensure_dirs()
    slug = f"{ranking.case}_{ranking.tag}"
    path = path or (REPORTS / f"rank_{slug}.md")
    lines = [
        f"# Apo discovery ranks — {ranking.case}",
        "",
        f"Ranker {ranking.extra.get('ranker_version', '')} on {ranking.tag} `{ranking.pdb_id}`. "
        f"Prior: {ranking.prior_source}. YAML linings used only in this table.",
        "",
        "| Site | Found | Geometry rank | Prior rank | Hybrid rank | Miner rank | Enrichment |",
        "|------|-------|---------------|------------|-------------|------------|------------|",
    ]
    for row in sites:
        miner = "" if row.miner_rank is None else row.miner_rank
        enr = "" if row.enrichment is None else f"{row.enrichment:.2f}"
        lines.append(
            f"| {row.site} | {row.found} | {row.geometry_rank or ''} | {row.prior_rank or ''} "
            f"| {row.hybrid_rank or ''} | {miner} | {enr} |"
        )
    path.write_text("\n".join(lines) + "\n")
    payload = {"ranking": ranking.as_dict(), "sites": [r.as_dict() for r in sites]}
    (PROCESSED / f"{slug}_rank.json").write_text(json.dumps(payload, indent=2))
    return path


def write_combined_rank_report(
    rows: list[tuple[ApoRanking, list[SiteEval]]],
    path: Path | None = None,
) -> Path:
    ensure_dirs()
    path = path or (REPORTS / "rank.md")
    lines = [
        "# Apo discovery ranks",
        "",
        "YAML linings used only after ranks were frozen. Lower rank is better.",
        "",
        "| Case | Site | Geometry | Prior | Hybrid | Miner | Enrichment |",
        "|------|------|----------|-------|--------|-------|------------|",
    ]
    for ranking, sites in rows:
        for row in sites:
            miner = "" if row.miner_rank is None else row.miner_rank
            enr = "" if row.enrichment is None else f"{row.enrichment:.2f}"
            geo = "" if row.geometry_rank is None else row.geometry_rank
            pri = "" if row.prior_rank is None else row.prior_rank
            hyb = "" if row.hybrid_rank is None else row.hybrid_rank
            lines.append(
                f"| {ranking.case} | {row.site} | {geo} | {pri} | {hyb} | {miner} | {enr} |"
            )
    path.write_text("\n".join(lines) + "\n")
    return path
