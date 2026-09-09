from __future__ import annotations

from pathlib import Path

from pocket_atlas.io.pdb import Structure
from pocket_atlas.paths import PROCESSED, ensure_dirs

DROP_RESNAMES = {"HOH", "WAT", "SO4", "PO4", "GOL", "EDO", "PEG", "CL", "NA"}


def prepare_structure(
    structure: Structure,
    chain: str = "A",
    keep_ligand: str | None = None,
) -> Structure:
    """Drop waters, unused chains, and non-ligand HETATMs. No physics yet.

    Hydrogens and minimization are optional OpenMM steps; the cryptic vs apo
    comparison is a heavy-atom geometric question and does not need them.
    """
    keep: list = []
    for atom in structure.atoms:
        if atom.chain != chain:
            continue
        if atom.resname in DROP_RESNAMES:
            continue
        if atom.is_het:
            if keep_ligand and atom.resname == keep_ligand:
                keep.append(atom)
            continue
        keep.append(atom)
    prepared = Structure(
        atoms=keep,
        pdb_id=structure.pdb_id,
        path=structure.path,
        extra={"chain": chain, "keep_ligand": keep_ligand},
    )
    return prepared


def write_pdb(structure: Structure, path: Path) -> Path:
    ensure_dirs()
    path.parent.mkdir(parents=True, exist_ok=True)
    lines = []
    for i, atom in enumerate(structure.atoms, start=1):
        rec = "HETATM" if atom.is_het else "ATOM  "
        lines.append(
            f"{rec}{i:5d} {atom.name:>4s} {atom.resname:>3s} {atom.chain:1s}"
            f"{atom.resseq:4d}    {atom.coord[0]:8.3f}{atom.coord[1]:8.3f}"
            f"{atom.coord[2]:8.3f}  1.00  0.00          {atom.element:>2s}"
        )
    lines.append("END")
    path.write_text("\n".join(lines) + "\n")
    return path


def prepared_path(pdb_id: str, tag: str = "prep") -> Path:
    return PROCESSED / f"{pdb_id}_{tag}.pdb"
