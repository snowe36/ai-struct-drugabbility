"""Leak-free apo cavity ranking. YAML cryptic linings stay out of this module."""

from __future__ import annotations

from dataclasses import dataclass, field

from pocket_atlas.cases import load_case
from pocket_atlas.dynamics.labels import resolve_prior
from pocket_atlas.io.pdb import Structure
from pocket_atlas.io.rcsb import load_structure
from pocket_atlas.paths import PROCESSED, ensure_dirs
from pocket_atlas.pockets import DETECTOR_VERSION, Pocket, detect_pockets
from pocket_atlas.prepare import prepare_structure, write_pdb

RANKER_VERSION = "0.1.0"


@dataclass
class ScoredPocket:
    pocket: Pocket
    enrichment: float
    hybrid: float
    geometry_rank: int = 0
    prior_rank: int = 0
    hybrid_rank: int = 0
    miner_mean: float | None = None
    miner_hybrid: float | None = None
    miner_rank: int | None = None

    def as_dict(self) -> dict:
        payload = {
            "index": self.pocket.index,
            "dscore": round(self.pocket.dscore, 3),
            "enrichment": round(self.enrichment, 3),
            "hybrid": round(self.hybrid, 3),
            "geometry_rank": self.geometry_rank,
            "prior_rank": self.prior_rank,
            "hybrid_rank": self.hybrid_rank,
            "volume": round(self.pocket.volume, 1),
            "seed_clearance": round(float(self.pocket.extra.get("seed_clearance", 0.0)), 3),
            "lining_residues": self.pocket.lining_residues,
            "centroid": [round(float(x), 2) for x in self.pocket.centroid],
        }
        if self.miner_mean is not None:
            payload["miner_mean"] = round(self.miner_mean, 3)
            payload["miner_hybrid"] = round(self.miner_hybrid or 0.0, 3)
            payload["miner_rank"] = self.miner_rank
        return payload


@dataclass
class ApoRanking:
    case: str
    tag: str
    pdb_id: str
    prior_source: str
    n_protein: int
    n_prior: int
    scored: list[ScoredPocket] = field(default_factory=list)
    extra: dict = field(default_factory=dict)

    @property
    def pockets(self) -> list[Pocket]:
        return [s.pocket for s in self.scored]

    def as_dict(self) -> dict:
        return {
            "case": self.case,
            "tag": self.tag,
            "pdb_id": self.pdb_id,
            "ranker_version": RANKER_VERSION,
            "detector_version": DETECTOR_VERSION,
            "prior_source": self.prior_source,
            "n_protein": self.n_protein,
            "n_prior": self.n_prior,
            "n_pockets": len(self.scored),
            "pockets": [s.as_dict() for s in self.scored],
            "extra": self.extra,
        }


def lining_enrichment(
    lining: set[int],
    prior: set[int],
    protein: set[int],
) -> float:
    """(prior fraction in lining) / (prior fraction in protein)."""
    lining = lining & protein
    prior = prior & protein
    if not lining or not protein or not prior:
        return 0.0
    frac = len(lining & prior) / len(lining)
    background = len(prior) / len(protein)
    if background <= 0:
        return 0.0
    return frac / background


def hybrid_score(dscore: float, enrichment: float) -> float:
    return float(dscore) * (1.0 + float(enrichment))


def _assign_ranks(scored: list[ScoredPocket], key, attr: str) -> None:
    order = sorted(scored, key=key, reverse=True)
    for rank, item in enumerate(order, start=1):
        setattr(item, attr, rank)


def rank_apo(
    pockets: list[Pocket],
    protein_residues: set[int],
    prior_residues: set[int],
    prior_source: str,
    case: str = "",
    tag: str = "apo",
    pdb_id: str = "",
    miner_scores: dict[int, float] | None = None,
) -> ApoRanking:
    """Rank apo cavities. Does not take a Case or YAML site linings."""
    protein = set(protein_residues)
    prior = set(prior_residues) & protein
    scored: list[ScoredPocket] = []
    for pocket in pockets:
        lining = set(pocket.lining_residues) & protein
        enr = lining_enrichment(lining, prior, protein)
        hyb = hybrid_score(pocket.dscore, enr)
        miner_mean = None
        miner_hyb = None
        if miner_scores:
            vals = [miner_scores[r] for r in lining if r in miner_scores]
            miner_mean = float(sum(vals) / len(vals)) if vals else 0.0
            miner_hyb = hybrid_score(pocket.dscore, miner_mean)
        scored.append(
            ScoredPocket(
                pocket=pocket,
                enrichment=enr,
                hybrid=hyb,
                miner_mean=miner_mean,
                miner_hybrid=miner_hyb,
            )
        )
    _assign_ranks(scored, lambda s: s.pocket.dscore, "geometry_rank")
    _assign_ranks(scored, lambda s: s.enrichment, "prior_rank")
    _assign_ranks(scored, lambda s: s.hybrid, "hybrid_rank")
    if miner_scores:
        _assign_ranks(scored, lambda s: s.miner_hybrid or 0.0, "miner_rank")
    scored.sort(key=lambda s: s.hybrid_rank)
    return ApoRanking(
        case=case,
        tag=tag,
        pdb_id=pdb_id,
        prior_source=prior_source,
        n_protein=len(protein),
        n_prior=len(prior),
        scored=scored,
        extra={"ranker_version": RANKER_VERSION},
    )


def apo_structure_spec(raw_structures: dict) -> tuple[str, dict]:
    """Apo arm, or crystal for closed-only cases. No site linings."""
    if "apo" in raw_structures:
        return "apo", raw_structures["apo"]
    return "crystal", raw_structures["crystal"]


def load_apo_arm(case_name: str) -> tuple[str, str, Structure, list[Pocket], set[int]]:
    """Fetch/prepare/detect the apo (or crystal) arm only."""
    ensure_dirs()
    case = load_case(case_name)
    tag, spec = apo_structure_spec(case.raw["structures"])
    pdb_id = spec["pdb_id"]
    structure = prepare_structure(load_structure(pdb_id), chain=case.chain)
    write_pdb(structure, PROCESSED / f"{case.name}_{tag}_{pdb_id}.pdb")
    pockets = detect_pockets(structure, chain=case.chain, ligand_mode="exclude")
    protein = set(structure.residue_numbers(chain=case.chain))
    return tag, pdb_id, structure, pockets, protein


def rank_case(
    case_name: str,
    prior: str = "literature",
    miner_scores: dict[int, float] | None = None,
    use_miner: bool = False,
) -> ApoRanking:
    """Apo ranking for a named case. Uses NMR/Dyna-1 prior, not site YAML."""
    case = load_case(case_name)
    tag, pdb_id, _structure, pockets, protein = load_apo_arm(case_name)
    prior_residues, source = resolve_prior(case, prior=prior, prefer_dyna1=prior == "dyna1")
    miner_note = None
    if miner_scores is None and use_miner:
        from pocket_atlas.rank.miner import MinerUnavailable, try_miner_scores

        try:
            miner_scores = try_miner_scores(case_name)
        except MinerUnavailable as exc:
            miner_note = str(exc)
            miner_scores = None
    ranking = rank_apo(
        pockets,
        protein_residues=protein,
        prior_residues=set(prior_residues),
        prior_source=source,
        case=case_name,
        tag=tag,
        pdb_id=pdb_id,
        miner_scores=miner_scores,
    )
    if miner_note:
        ranking.extra["miner"] = miner_note
    return ranking
