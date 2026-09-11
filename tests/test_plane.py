from pathlib import Path

import numpy as np

from pocket_atlas.pockets import Pocket
from pocket_atlas.rank import rank_apo
from pocket_atlas.rank.plane import load_panel, plane_row_stats, write_plane_report


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


def test_aqadk_crystal_is_sequence_matched():
    panel = load_panel()
    aqadk = next(p for p in panel["proteins"] if p["id"] == "aqadk")
    assert aqadk["pdb_id"] == "2RH5"
    assert len(panel["proteins"]) == 10
    assert panel["literature_extra"][0]["id"] == "tem1_horn"


def test_plane_source_has_no_site_yaml():
    from pocket_atlas.rank import plane as plane_mod

    text = Path(plane_mod.__file__).read_text()
    for token in ("cryptic_site", "omega_loop", "nucleotide_site", "catalytic_site"):
        assert token not in text


def test_plane_row_stats_and_report(tmp_path):
    protein = set(range(1, 11))
    prior = {1, 2}
    ranking = rank_apo(
        [_pocket(0, [1, 2, 3, 4], 1.0), _pocket(1, [9, 10], 2.0)],
        protein_residues=protein,
        prior_residues=prior,
        prior_source="relaxdb_cpmg",
        case="cypa",
        pdb_id="2CPL",
    )
    stats = plane_row_stats(ranking)
    assert stats["max_enrichment"] == 2.5
    assert stats["dscore_at_max_enr"] == 1.0
    assert stats["enrichment_at_max_dscore"] == 0.0
    ranking.extra["mapped_frac"] = 0.99
    path = write_plane_report([ranking], path=tmp_path / "plane.md")
    text = path.read_text()
    assert "CypA" in text
    assert "2.50" in text
