from __future__ import annotations

import numpy as np

from pocket_atlas.dynamics import overlap_case
from pocket_atlas.io.pdb import Atom, Structure
from pocket_atlas.pockets import detect_pockets, pocket_near_residues


def _atom(resseq: int, x: float, y: float, z: float, name: str = "CA", resname: str = "LEU") -> Atom:
    return Atom(
        name=name,
        resname=resname,
        resseq=resseq,
        chain="A",
        coord=np.array([x, y, z], dtype=float),
        element=name[:1],
        is_het=False,
    )


def test_detects_synthetic_buried_void():
    """A hollow hydrophobic shell should produce one interior pocket."""
    atoms = []
    res = 1
    # Two concentric-ish layers of a 10 Å cube, empty 4 Å core.
    for x in range(-8, 9, 2):
        for y in range(-8, 9, 2):
            for z in range(-8, 9, 2):
                if max(abs(x), abs(y), abs(z)) < 4:
                    continue
                if max(abs(x), abs(y), abs(z)) > 8:
                    continue
                atoms.append(_atom(res, float(x), float(y), float(z)))
                res += 1
    structure = Structure(atoms=atoms, pdb_id="void")
    pockets = detect_pockets(structure, min_neighbors=12, min_points=6, grid=1.0)
    assert pockets, "expected a void in the hollow shell"
    assert pockets[0].volume >= 8.0


def test_open_vs_closed_shell():
    """Removing one face of the shell should shrink or abolish the void."""
    def shell(open_face: bool) -> Structure:
        atoms = []
        res = 1
        for x in range(-6, 7, 2):
            for y in range(-6, 7, 2):
                for z in range(-6, 7, 2):
                    if max(abs(x), abs(y), abs(z)) < 3:
                        continue
                    if open_face and x > 4:
                        continue
                    atoms.append(_atom(res, float(x), float(y), float(z)))
                    res += 1
        return Structure(atoms=atoms)

    closed = detect_pockets(shell(False), min_neighbors=10, min_points=5, grid=1.0)
    opened = detect_pockets(shell(True), min_neighbors=10, min_points=5, grid=1.0)
    closed_vol = closed[0].volume if closed else 0.0
    opened_vol = opened[0].volume if opened else 0.0
    assert closed_vol >= opened_vol


def test_overlap_enrichment_math():
    from pocket_atlas.cases import Case

    case = Case(
        name="toy",
        raw={
            "cryptic_site": {"residues": [1, 2, 3, 4]},
            "nmr_exchange": {"residues": [1, 2, 10]},
            "catalytic_site": {"residues": [20, 21]},
        },
    )
    protein = set(range(1, 31))
    result = overlap_case(case, protein_residues=protein, source="test")
    assert result.cryptic_and_nmr == 2
    assert result.control_and_nmr == 0
    assert result.cryptic_enrichment > result.control_enrichment
    assert result.odds_ratio > 1 or result.odds_ratio == float("inf")


def test_pocket_near_residues_requires_overlap():
    from pocket_atlas.pockets import Pocket

    pocket = Pocket(
        index=0,
        points=np.zeros((5, 3)),
        volume=20,
        enclosure=0.5,
        hydrophobicity=1.0,
        polarity=0.2,
        centroid=np.zeros(3),
        lining_residues=[10, 11, 12, 13],
        dscore=1.0,
    )
    other = Pocket(
        index=1,
        points=np.zeros((5, 3)),
        volume=80,
        enclosure=0.4,
        hydrophobicity=0.0,
        polarity=0.5,
        centroid=np.zeros(3),
        lining_residues=[10, 11],
        dscore=2.0,
    )
    assert pocket_near_residues([pocket], {10, 11, 12}, min_overlap=3) is pocket
    assert pocket_near_residues([other, pocket], {10, 11, 12}, min_overlap=2) is pocket
    assert pocket_near_residues([pocket], {99, 100}, min_overlap=3) is None


def test_pdb_parser_skips_altloc_b():
    from pocket_atlas.io.pdb import parse_pdb

    text = (
        "ATOM      1  CA  ALA A   1       0.000   0.000   0.000  1.00 10.00           C  \n"
        "ATOM      2  CA BALA A   1       1.000   0.000   0.000  1.00 10.00           C  \n"
    )
    structure = parse_pdb(text, pdb_id="x")
    cas = [a for a in structure.atoms if a.name == "CA"]
    assert len(cas) == 1
    assert cas[0].coord[0] == 0.0
