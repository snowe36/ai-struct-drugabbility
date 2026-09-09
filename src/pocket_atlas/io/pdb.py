from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path

import numpy as np

PROTEIN_RESIDUES = {
    "ALA", "ARG", "ASN", "ASP", "CYS", "GLN", "GLU", "GLY", "HIS", "ILE",
    "LEU", "LYS", "MET", "PHE", "PRO", "SER", "THR", "TRP", "TYR", "VAL",
    "MSE", "SEC",
}


@dataclass
class Atom:
    name: str
    resname: str
    resseq: int
    chain: str
    coord: np.ndarray
    element: str
    is_het: bool = False


@dataclass
class Structure:
    atoms: list[Atom]
    pdb_id: str = ""
    path: str = ""
    extra: dict = field(default_factory=dict)

    def protein_atoms(self, chain: str | None = None, heavy: bool = True) -> list[Atom]:
        out = []
        for atom in self.atoms:
            if atom.is_het:
                continue
            if atom.resname not in PROTEIN_RESIDUES:
                continue
            if chain is not None and atom.chain != chain:
                continue
            if heavy and atom.element == "H":
                continue
            out.append(atom)
        return out

    def coords(self, atoms: list[Atom] | None = None) -> np.ndarray:
        chosen = atoms if atoms is not None else self.protein_atoms()
        if not chosen:
            return np.zeros((0, 3), dtype=float)
        return np.vstack([a.coord for a in chosen])

    def residue_ca(self, chain: str | None = None) -> dict[int, np.ndarray]:
        cas: dict[int, np.ndarray] = {}
        for atom in self.protein_atoms(chain=chain):
            if atom.name == "CA":
                cas[atom.resseq] = atom.coord
        return cas

    def ligand_atoms(self, resname: str, chain: str | None = None) -> list[Atom]:
        out = []
        for atom in self.atoms:
            if atom.resname != resname:
                continue
            if chain is not None and atom.chain != chain:
                continue
            if atom.element == "H":
                continue
            out.append(atom)
        return out

    def residue_numbers(self, chain: str | None = None) -> list[int]:
        return sorted({a.resseq for a in self.protein_atoms(chain=chain)})


def parse_pdb(text: str, pdb_id: str = "") -> Structure:
    atoms: list[Atom] = []
    for line in text.splitlines():
        rec = line[:6]
        if rec not in ("ATOM  ", "HETATM"):
            continue
        if len(line) < 54:
            continue
        alt = line[16]
        if alt not in (" ", "A"):
            continue
        name = line[12:16].strip()
        resname = line[17:20].strip()
        chain = line[21].strip() or "A"
        try:
            resseq = int(line[22:26])
        except ValueError:
            continue
        x = float(line[30:38])
        y = float(line[38:46])
        z = float(line[46:54])
        if len(line) >= 78 and line[76:78].strip():
            element = line[76:78].strip().title()
        else:
            element = "".join(c for c in name if c.isalpha())[:1] or "C"
        atoms.append(
            Atom(
                name=name,
                resname=resname,
                resseq=resseq,
                chain=chain,
                coord=np.array([x, y, z], dtype=float),
                element=element,
                is_het=rec == "HETATM" or resname not in PROTEIN_RESIDUES,
            )
        )
    return Structure(atoms=atoms, pdb_id=pdb_id)


def read_pdb(path: str | bytes, pdb_id: str = "") -> Structure:
    p = Path(path)
    text = p.read_text()
    structure = parse_pdb(text, pdb_id=pdb_id or p.stem)
    structure.path = str(p)
    return structure
