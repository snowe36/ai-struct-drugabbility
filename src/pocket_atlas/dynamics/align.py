"""Map experimental CPMG sequence indices onto crystal resseq."""

from __future__ import annotations

from pocket_atlas.io.pdb import Structure

AA3_TO_1 = {
    "ALA": "A", "ARG": "R", "ASN": "N", "ASP": "D", "CYS": "C",
    "GLN": "Q", "GLU": "E", "GLY": "G", "HIS": "H", "ILE": "I",
    "LEU": "L", "LYS": "K", "MET": "M", "PHE": "F", "PRO": "P",
    "SER": "S", "THR": "T", "TRP": "W", "TYR": "Y", "VAL": "V",
    "MSE": "M", "SEC": "C",
}


def chain_ca_sequence(structure: Structure, chain: str) -> tuple[list[int], str]:
    resseqs: list[int] = []
    letters: list[str] = []
    seen: set[int] = set()
    for atom in structure.protein_atoms(chain=chain, heavy=True):
        if atom.name != "CA" or atom.resseq in seen:
            continue
        seen.add(atom.resseq)
        resseqs.append(atom.resseq)
        letters.append(AA3_TO_1.get(atom.resname, "X"))
    return resseqs, "".join(letters)


def map_cpmg_index_to_resseq(cpmg_seq: str, resseqs: list[int], pdb_seq: str) -> dict[int, int]:
    """1-based CPMG sequence index → PDB resseq. Empty if the sequences do not share a seed."""
    if not cpmg_seq or not pdb_seq or not resseqs:
        return {}
    if cpmg_seq in pdb_seq:
        off = pdb_seq.index(cpmg_seq)
        return {i + 1: resseqs[off + i] for i in range(len(cpmg_seq)) if off + i < len(resseqs)}
    if pdb_seq in cpmg_seq:
        off = cpmg_seq.index(pdb_seq)
        return {off + i + 1: resseqs[i] for i in range(len(pdb_seq))}
    mapping: dict[int, int] = {}
    for seed in (20, 12, 8):
        for i in range(0, len(cpmg_seq) - seed + 1):
            j = pdb_seq.find(cpmg_seq[i : i + seed])
            if j < 0:
                continue
            a, b = i, j
            while a < len(cpmg_seq) and b < len(pdb_seq) and cpmg_seq[a] == pdb_seq[b]:
                mapping[a + 1] = resseqs[b]
                a += 1
                b += 1
            a, b = i - 1, j - 1
            while a >= 0 and b >= 0 and cpmg_seq[a] == pdb_seq[b]:
                mapping[a + 1] = resseqs[b]
                a -= 1
                b -= 1
            if len(mapping) >= 0.5 * min(len(cpmg_seq), len(pdb_seq)):
                return mapping
        mapping = {}
    return mapping
