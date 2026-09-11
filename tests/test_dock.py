from __future__ import annotations

import numpy as np
import pytest

from pocket_atlas.cases import load_case
from pocket_atlas.dock import (
    DockUnavailable,
    case_boxes,
    find_gnina,
    parse_gnina_sdf,
    residue_box,
    rmsd_heavy,
    write_ligand_pdb,
    write_ligand_xyz,
)
from pocket_atlas.io.pdb import Atom, Structure


def _ca(resseq: int, x: float, y: float, z: float) -> Atom:
    return Atom(
        name="CA",
        resname="ALA",
        resseq=resseq,
        chain="A",
        coord=np.array([x, y, z], dtype=float),
        element="C",
        is_het=False,
    )


def test_residue_box_centroid_and_padding():
    atoms = [_ca(1, 0, 0, 0), _ca(2, 10, 0, 0), _ca(3, 0, 10, 0)]
    structure = Structure(atoms=atoms, pdb_id="toy")
    box = residue_box(structure, {1, 2, 3}, "cryptic", pad=2.0, min_size=4.0)
    assert box.center == (5.0, 5.0, 0.0)
    assert box.size[0] == pytest.approx(14.0)
    assert box.size[1] == pytest.approx(14.0)
    assert box.size[2] >= 4.0
    assert box.residues == [1, 2, 3]


def test_case_boxes_have_cryptic_nmr_control():
    case = load_case("tem1_horn")
    residues = sorted(
        case.cryptic_residues | case.nmr_residues | case.residue_set("catalytic_site")
    )
    atoms = [_ca(r, float(i), 0.0, 0.0) for i, r in enumerate(residues)]
    structure = Structure(atoms=atoms, pdb_id="fake")
    boxes = case_boxes(case, structure, prior="literature")
    assert set(boxes) == {"cryptic", "nmr", "control"}
    assert boxes["cryptic"].residues
    assert boxes["control"].residues
    assert 70 in boxes["control"].residues


def test_write_ligand_pdb_and_xyz(tmp_path):
    atoms = [
        Atom("C1", "CBT", 300, "A", np.array([1.0, 2.0, 3.0]), "C", True),
        Atom("O1", "CBT", 300, "A", np.array([1.5, 2.0, 3.0]), "O", True),
    ]
    pdb = write_ligand_pdb(atoms, tmp_path / "lig.pdb")
    xyz = write_ligand_xyz(atoms, tmp_path / "lig.xyz")
    text = pdb.read_text()
    assert "CBT" in text and "HETATM" in text
    xyz_lines = xyz.read_text().splitlines()
    assert xyz_lines[0] == "2"
    assert xyz_lines[2].startswith("C")


def test_parse_gnina_sdf_and_rmsd():
    sdf = """
  fake
     RDKit          3D

  2  0  0  0  0  0  0  0  0  0999 V2000
    0.0000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    1.0000    0.0000    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
M  END
> <minimizedAffinity>
-7.50

> <CNNscore>
0.42

> <CNNaffinity>
5.10

$$$$
"""
    poses = parse_gnina_sdf(sdf)
    assert len(poses) == 1
    assert poses[0].vina == pytest.approx(-7.5)
    assert poses[0].cnn_affinity == pytest.approx(5.1)
    assert poses[0].coords is not None
    ref = np.array([[0.0, 0.0, 0.0], [1.0, 0.0, 0.0]])
    assert rmsd_heavy(ref, poses[0].coords) == pytest.approx(0.0)


def test_find_gnina_fails_closed():
    import pocket_atlas.dock as dock

    original = dock.shutil.which
    dock.shutil.which = lambda _name: None
    try:
        with pytest.raises(DockUnavailable, match="gnina"):
            find_gnina()
    finally:
        dock.shutil.which = original


def test_remap_dyna1_index_to_pdb_resseq(tmp_path):
    from pocket_atlas.dynamics.dyna1 import remap_index_scores
    from pocket_atlas.prepare import write_pdb

    atoms = [
        Atom("CA", "ASP", 218, "A", np.array([0.0, 0.0, 0.0]), "C"),
        Atom("CA", "ILE", 219, "A", np.array([1.0, 0.0, 0.0]), "C"),
        Atom("CA", "SER", 220, "A", np.array([2.0, 0.0, 0.0]), "C"),
    ]
    path = tmp_path / "toy.pdb"
    write_pdb(Structure(atoms=atoms), path)
    mapped = remap_index_scores({1: 0.1, 2: 0.9, 3: 0.2}, path, chain="A")
    assert mapped == {218: 0.1, 219: 0.9, 220: 0.2}
