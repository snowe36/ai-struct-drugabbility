from __future__ import annotations

import inspect
from pathlib import Path

import numpy as np
import pytest

from pocket_atlas.cases import load_case
from pocket_atlas.pockets import Pocket
from pocket_atlas.rank import RANKER_VERSION, hybrid_score, lining_enrichment, rank_apo
from pocket_atlas.rank.eval import eval_labeled_sites
from pocket_atlas.rank.miner import MinerUnavailable, scores_from_bfactor_pdb, try_miner_scores


def _pocket(index: int, lining: list[int], dscore: float) -> Pocket:
    return Pocket(
        index=index,
        points=np.zeros((5, 3)),
        volume=40.0,
        enclosure=0.5,
        hydrophobicity=0.0,
        polarity=0.2,
        centroid=np.zeros(3),
        lining_residues=lining,
        dscore=dscore,
    )


def test_ranker_version_frozen():
    assert RANKER_VERSION == "0.1.0"


def test_rank_apo_source_has_no_site_yaml():
    from pocket_atlas import rank as rank_mod

    text = Path(rank_mod.__file__).read_text()
    src = inspect.getsource(rank_apo)
    for token in ("cryptic_site", "omega_loop", "nucleotide_site", "catalytic_site"):
        assert token not in text
        assert token not in src
    appearance = Path(rank_mod.__file__).parent / "appearance.py"
    assert "cryptic_site" not in appearance.read_text()
    miner = Path(rank_mod.__file__).parent / "miner.py"
    assert "cryptic_site" not in miner.read_text()


def test_lining_enrichment_and_hybrid():
    protein = set(range(1, 11))
    prior = {1, 2}
    lining = {1, 2, 3, 4}
    # 2/4 over 2/10 = 2.5
    assert lining_enrichment(lining, prior, protein) == pytest.approx(2.5)
    assert hybrid_score(2.0, 2.5) == pytest.approx(7.0)
    assert lining_enrichment({9, 10}, prior, protein) == 0.0


def test_tem1_hybrid_promotes_omega_over_horn():
    case = load_case("tem1_horn")
    protein = set(range(1, 291))
    prior = set(case.nmr_residues)
    horn = _pocket(0, sorted(case.cryptic_residues), dscore=2.0)
    omega = _pocket(1, sorted(case.residue_set("omega_loop")), dscore=1.0)
    cat = _pocket(2, sorted(case.residue_set("catalytic_site")), dscore=1.4)
    ranking = rank_apo(
        [horn, omega, cat],
        protein_residues=protein,
        prior_residues=prior,
        prior_source="nmr_literature",
        case="tem1_horn",
        pdb_id="synthetic",
    )
    by_idx = {s.pocket.index: s for s in ranking.scored}
    assert by_idx[0].geometry_rank < by_idx[1].geometry_rank
    sites = {row.site: row for row in eval_labeled_sites(ranking)}
    assert sites["omega_loop"].hybrid_rank < sites["horn"].hybrid_rank
    assert sites["horn"].found and sites["omega_loop"].found


def test_kras_hybrid_promotes_switch2_over_nucleotide():
    case = load_case("kras_switch2")
    protein = set(range(1, 171))
    prior = set(case.nmr_residues)
    switch2 = _pocket(0, sorted(case.cryptic_residues), dscore=1.2)
    nucleotide = _pocket(1, sorted(case.residue_set("nucleotide_site")), dscore=1.8)
    ranking = rank_apo(
        [switch2, nucleotide],
        protein_residues=protein,
        prior_residues=prior,
        prior_source="nmr_literature",
        case="kras_switch2",
        pdb_id="synthetic",
    )
    by_idx = {s.pocket.index: s for s in ranking.scored}
    assert by_idx[1].geometry_rank < by_idx[0].geometry_rank
    sites = {row.site: row for row in eval_labeled_sites(ranking)}
    assert sites["switch2"].hybrid_rank < sites["nucleotide"].hybrid_rank


def test_miner_fail_closed_without_weights():
    with pytest.raises(MinerUnavailable):
        try_miner_scores("tem1_horn")


def test_miner_bfactor_pdb_and_lining_mean(tmp_path):
    pdb = tmp_path / "miner.pdb"
    pdb.write_text(
        "ATOM      1  CA  ALA A   1       0.000   0.000   0.000  1.00  0.90           C  \n"
        "ATOM      2  CA  ALA A   2       1.000   0.000   0.000  1.00  0.10           C  \n"
        "ATOM      3  CA  ALA A   3       2.000   0.000   0.000  1.00  0.50           C  \n"
    )
    scores = scores_from_bfactor_pdb(pdb)
    assert scores[1] == pytest.approx(0.90)
    pocket = _pocket(0, [1, 2], dscore=1.0)
    other = _pocket(1, [3], dscore=1.0)
    ranking = rank_apo(
        [pocket, other],
        protein_residues={1, 2, 3},
        prior_residues={1},
        prior_source="test",
        miner_scores=scores,
    )
    by_idx = {s.pocket.index: s for s in ranking.scored}
    assert by_idx[0].miner_mean == pytest.approx(0.5)
    assert by_idx[0].miner_rank == 1


def test_tem1_kras_live_ranks_skip_if_rcsb_down():
    """Record apo ranks on real crystals. Skip if RCSB is unreachable."""
    from pocket_atlas.rank import rank_case

    try:
        tem = rank_case("tem1_horn", prior="literature")
        kras = rank_case("kras_switch2", prior="literature")
    except Exception as exc:
        pytest.skip(f"RCSB apo PDBs unavailable: {exc}")

    tem_sites = {row.site: row for row in eval_labeled_sites(tem)}
    kras_sites = {row.site: row for row in eval_labeled_sites(kras)}
    assert tem.scored and kras.scored
    assert "horn" in tem_sites and "switch2" in kras_sites


def test_vp35_open_frame_skips_without_xtc():
    from pocket_atlas.rank.appearance import FrameUnavailable, find_strided_traj, rank_open_frame

    if find_strided_traj() is not None:
        pytest.skip("strided VP35 xtc is present")
    with pytest.raises(FrameUnavailable):
        rank_open_frame()
