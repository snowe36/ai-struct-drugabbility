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


def test_ligand_include_shrinks_occupied_void():
    """A HETATM in the void should survive exclude and shrink/kill include."""
    atoms = []
    res = 1
    for x in range(-8, 9, 2):
        for y in range(-8, 9, 2):
            for z in range(-8, 9, 2):
                if max(abs(x), abs(y), abs(z)) < 4:
                    continue
                if max(abs(x), abs(y), abs(z)) > 8:
                    continue
                atoms.append(_atom(res, float(x), float(y), float(z)))
                res += 1
    ligands = []
    for x in range(-4, 5):
        for y in range(-4, 5):
            for z in range(-4, 5):
                ligands.append(
                    Atom(
                        name="C1",
                        resname="LIG",
                        resseq=999,
                        chain="A",
                        coord=np.array([float(x), float(y), float(z)], dtype=float),
                        element="C",
                        is_het=True,
                    )
                )
    structure = Structure(atoms=atoms + ligands, pdb_id="void")
    excluded = detect_pockets(structure, ligand_mode="exclude", min_neighbors=12, min_points=6)
    included = detect_pockets(structure, ligand_mode="include", min_neighbors=12, min_points=6)
    assert excluded, "exclude mode should still see the empty void"
    assert excluded[0].extra.get("ligand_mode") == "exclude"

    def _near_origin(pockets):
        for pocket in pockets:
            if float(np.linalg.norm(pocket.centroid)) < 3.0:
                return pocket
        return None

    assert _near_origin(excluded) is not None
    # Occupied core: include must not report a void at the origin.
    assert _near_origin(included) is None


def test_horn_cluster_exists_and_is_not_the_chain():
    """1BTL/1PZO: horn cluster exists, ligand-scale, not the whole protein."""
    import pytest

    from pocket_atlas.cases import load_case
    from pocket_atlas.io.rcsb import load_structure
    from pocket_atlas.pockets import DETECTOR_VERSION
    from pocket_atlas.prepare import prepare_structure

    assert DETECTOR_VERSION == "0.2.0"
    case = load_case("tem1_horn")
    try:
        apo = prepare_structure(load_structure("1BTL"), chain="A")
        holo = prepare_structure(load_structure("1PZO"), chain="A", keep_ligand="CBT")
    except Exception as exc:
        pytest.skip(f"RCSB 1BTL/1PZO unavailable: {exc}")

    n_res = len(apo.residue_numbers(chain="A"))
    assert n_res > 200
    for tag, structure, mode in (("apo", apo, "exclude"), ("holo", holo, "exclude")):
        pockets = detect_pockets(structure, chain="A", ligand_mode=mode)
        site = pocket_near_residues(pockets, set(case.cryptic_residues))
        assert site is not None, f"{tag}: horn cluster missing"
        assert 200 < site.volume < 8000, f"{tag}: volume {site.volume} is not ligand-scale"
        assert len(site.lining_residues) < n_res * 0.4, f"{tag}: lining is the chain"
        assert site.extra.get("seed_clearance", 0) >= 2.6
    included = detect_pockets(holo, chain="A", ligand_mode="include")
    horn_incl = pocket_near_residues(included, set(case.cryptic_residues))
    # Include is a different measurement, not a volume contest.
    assert horn_incl is None or horn_incl.extra.get("ligand_mode") == "include"


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
