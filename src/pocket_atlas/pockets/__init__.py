"""Ligand-scale pocket detection.

SiteMap-inspired geometry, not SiteMap. Packing interstices (empty space
~1.5 Å between atoms) percolate through a protein core if you 6-connect
them. Real ligandable cavities are *local maxima* of the clearance field
in a ligand-scale band (~2.6–5 Å). We keep those maxima, cluster them,
and gather nearby empty grid points. That is an fpocket-style alpha-sphere
idea without shipping fpocket.
"""

from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np
from scipy.spatial import cKDTree

from pocket_atlas.io.pdb import Atom, Structure

DETECTOR_VERSION = "0.2.0"

KYTE_DOOLITTLE = {
    "ILE": 4.5, "VAL": 4.2, "LEU": 3.8, "PHE": 2.8, "CYS": 2.5, "MET": 1.9,
    "ALA": 1.8, "GLY": -0.4, "THR": -0.7, "SER": -0.8, "TRP": -0.9,
    "TYR": -1.3, "PRO": -1.6, "HIS": -3.2, "GLU": -3.5, "GLN": -3.5,
    "ASP": -3.5, "ASN": -3.5, "LYS": -3.9, "ARG": -4.5, "MSE": 1.9,
}


@dataclass
class Pocket:
    index: int
    points: np.ndarray
    volume: float
    enclosure: float
    hydrophobicity: float
    polarity: float
    centroid: np.ndarray
    lining_residues: list[int]
    dscore: float
    extra: dict = field(default_factory=dict)

    def as_dict(self) -> dict:
        payload = {
            "index": self.index,
            "volume": round(self.volume, 1),
            "enclosure": round(self.enclosure, 3),
            "hydrophobicity": round(self.hydrophobicity, 3),
            "polarity": round(self.polarity, 3),
            "dscore": round(self.dscore, 3),
            "centroid": [round(float(x), 2) for x in self.centroid],
            "lining_residues": self.lining_residues,
            "n_points": int(len(self.points)),
        }
        if self.extra:
            payload["extra"] = {
                k: (round(v, 3) if isinstance(v, float) else v)
                for k, v in self.extra.items()
            }
        return payload


def _dscore(volume: float, enclosure: float, hydrophobicity: float) -> float:
    """Documented analogue, not SiteMap Dscore.

    A buried hydrophobic void (~150 Å³, enclosure 0.7, KD > 0) maps toward
    ~1.0. Solvent-exposed polar dimples sit well below 0.8.
    """
    vol_term = min(volume / 200.0, 1.5)
    hyd_term = (hydrophobicity + 4.5) / 9.0
    return float(0.45 * vol_term + 0.35 * enclosure + 0.20 * hyd_term)


def detect_pockets(
    structure: Structure,
    chain: str = "A",
    grid: float = 1.0,
    probe: float = 1.4,
    r_lo: float = 2.6,
    r_hi: float = 5.0,
    cluster_cut: float = 4.0,
    gather_radius: float = 5.5,
    neighbor_radius: float = 8.0,
    min_neighbors: int = 16,
    min_points: int = 10,
    lining_cutoff: float = 4.5,
    ligand_mode: str = "exclude",
) -> list[Pocket]:
    """ligand_mode: exclude = protein only; include = protein + HETATMs in the field."""
    if ligand_mode not in {"exclude", "include"}:
        raise ValueError(f"ligand_mode must be exclude|include, got {ligand_mode!r}")
    protein = structure.protein_atoms(chain=chain, heavy=True)
    if ligand_mode == "include":
        ligands = [
            a
            for a in structure.atoms
            if a.is_het and a.element != "H" and (chain is None or a.chain == chain)
        ]
        atoms = protein + ligands
    else:
        atoms = protein
    if len(atoms) < 10:
        return []
    coords = structure.coords(atoms)
    tree = cKDTree(coords)

    pad = 4.0
    mins = coords.min(axis=0) - pad
    maxs = coords.max(axis=0) + pad
    xs = np.arange(mins[0], maxs[0], grid)
    ys = np.arange(mins[1], maxs[1], grid)
    zs = np.arange(mins[2], maxs[2], grid)
    shape = (len(xs), len(ys), len(zs))
    gx, gy, gz = np.meshgrid(xs, ys, zs, indexing="ij")
    points = np.column_stack([gx.ravel(), gy.ravel(), gz.ravel()])

    dist, _ = tree.query(points, k=1, workers=-1)
    dist3 = dist.reshape(shape)

    maxima = _local_maxima(dist3, xs, ys, zs, r_lo=r_lo, r_hi=r_hi)
    if len(maxima) == 0:
        return []

    n_near = tree.query_ball_point(maxima, r=neighbor_radius, workers=-1, return_length=True)
    buried = maxima[n_near >= min_neighbors]
    if len(buried) == 0:
        return []

    clusters = _cluster_points(buried, cutoff=cluster_cut)
    pockets: list[Pocket] = []
    for members in clusters:
        pts, enclosure = _gather_cavity(
            points,
            dist,
            members,
            tree,
            probe=probe,
            r_hi=r_hi,
            gather_radius=gather_radius,
            neighbor_radius=neighbor_radius,
            min_neighbors=min_neighbors,
        )
        if len(pts) < min_points:
            continue
        seed_clearance = float(tree.query(members, k=1, workers=-1)[0].max())
        # A lone packing maximum (~2.6 Å) is not a ligandable cavity.
        if len(members) == 1 and seed_clearance < 3.0:
            continue
        volume = float(len(pts) * grid**3)
        centroid = pts.mean(axis=0)
        lining = _lining_residues(protein, pts, cutoff=lining_cutoff)
        hyd, pol = _lining_chemistry(protein, lining)
        r_max = float(np.max(np.linalg.norm(members - members.mean(axis=0), axis=1))) if len(members) > 1 else 0.0
        pockets.append(
            Pocket(
                index=0,
                points=pts,
                volume=volume,
                enclosure=enclosure,
                hydrophobicity=hyd,
                polarity=pol,
                centroid=centroid,
                lining_residues=lining,
                dscore=_dscore(volume, enclosure, hyd),
                extra={
                    "n_spheres": int(len(members)),
                    "seed_clearance": seed_clearance,
                    "cluster_span": r_max,
                    "ligand_mode": ligand_mode,
                    "detector_version": DETECTOR_VERSION,
                },
            )
        )
    pockets.sort(key=lambda p: p.dscore, reverse=True)
    for i, pocket in enumerate(pockets):
        pocket.index = i
    return pockets


def _local_maxima(
    dist3: np.ndarray,
    xs: np.ndarray,
    ys: np.ndarray,
    zs: np.ndarray,
    r_lo: float,
    r_hi: float,
) -> np.ndarray:
    """26-connected local maxima of clearance in the ligand-scale band."""
    band = (dist3 > r_lo) & (dist3 < r_hi)
    if not np.any(band):
        return np.zeros((0, 3))

    padded = np.pad(dist3, 1, mode="constant", constant_values=-np.inf)
    max_nbr = np.full_like(dist3, -np.inf)
    for dx in (-1, 0, 1):
        for dy in (-1, 0, 1):
            for dz in (-1, 0, 1):
                if dx == dy == dz == 0:
                    continue
                sl = padded[1 + dx : 1 + dx + dist3.shape[0], 1 + dy : 1 + dy + dist3.shape[1], 1 + dz : 1 + dz + dist3.shape[2]]
                max_nbr = np.maximum(max_nbr, sl)
    is_max = band & (dist3 >= max_nbr)
    ii, jj, kk = np.where(is_max)
    if len(ii) == 0:
        return np.zeros((0, 3))
    return np.column_stack([xs[ii], ys[jj], zs[kk]])


def _cluster_points(points: np.ndarray, cutoff: float) -> list[np.ndarray]:
    """Connected components with cutoff as the edge length."""
    order = np.lexsort((points[:, 0], points[:, 1], points[:, 2]))
    unused = np.ones(len(points), dtype=bool)
    tree = cKDTree(points)
    clusters: list[np.ndarray] = []
    for i in order:
        if not unused[i]:
            continue
        idxs = tree.query_ball_point(points[i], r=cutoff)
        # Grow once from the seed so a 4 Å chain still clusters.
        seen = set(idxs)
        queue = list(idxs)
        while queue:
            j = queue.pop()
            for k in tree.query_ball_point(points[j], r=cutoff):
                if k not in seen and unused[k]:
                    seen.add(k)
                    queue.append(k)
        members = np.array(sorted(seen))
        unused[members] = False
        clusters.append(points[members])
    return clusters


def _gather_cavity(
    points: np.ndarray,
    dist: np.ndarray,
    seeds: np.ndarray,
    atom_tree: cKDTree,
    probe: float,
    r_hi: float,
    gather_radius: float,
    neighbor_radius: float,
    min_neighbors: int,
) -> tuple[np.ndarray, float]:
    """Empty grid points near ligand-scale seeds, clipped so cores cannot percolate."""
    seed_tree = cKDTree(seeds)
    d_seed, _ = seed_tree.query(points, k=1, workers=-1)
    mask = (dist > probe) & (dist < r_hi) & (d_seed <= gather_radius)
    candidates = points[mask]
    if len(candidates) == 0:
        return np.zeros((0, 3)), 0.0
    n_near = atom_tree.query_ball_point(
        candidates, r=neighbor_radius, workers=-1, return_length=True
    )
    keep = n_near >= min_neighbors
    kept = candidates[keep]
    if len(kept) == 0:
        return np.zeros((0, 3)), 0.0
    enclosure = float(np.clip(n_near[keep] / 80.0, 0.0, 1.0).mean())
    return kept, enclosure


def _lining_residues(atoms: list[Atom], points: np.ndarray, cutoff: float) -> list[int]:
    coords = np.vstack([a.coord for a in atoms])
    tree = cKDTree(coords)
    hits = tree.query_ball_point(points, r=cutoff, workers=-1)
    residues = {atoms[j].resseq for idxs in hits for j in idxs}
    return sorted(residues)


def _lining_chemistry(atoms: list[Atom], lining: list[int]) -> tuple[float, float]:
    resnames: dict[int, str] = {}
    for atom in atoms:
        if atom.resseq in lining and atom.name == "CA":
            resnames[atom.resseq] = atom.resname
    if not resnames:
        return 0.0, 0.0
    kd = [KYTE_DOOLITTLE.get(name, 0.0) for name in resnames.values()]
    hyd = float(np.mean(kd))
    polar = {"SER", "THR", "ASN", "GLN", "ASP", "GLU", "LYS", "ARG", "HIS", "TYR"}
    polarity = float(sum(1 for n in resnames.values() if n in polar) / len(resnames))
    return hyd, polarity


def pocket_near_residues(
    pockets: list[Pocket],
    residues: set[int],
    min_overlap: int = 3,
) -> Pocket | None:
    """Pocket whose lining overlaps `residues` most; Dscore breaks ties."""
    best: Pocket | None = None
    best_overlap = 0
    for pocket in pockets:
        n_overlap = len(set(pocket.lining_residues) & residues)
        if n_overlap < min_overlap:
            continue
        better_overlap = n_overlap > best_overlap
        tie = n_overlap == best_overlap and best is not None and pocket.dscore > best.dscore
        if best is None or better_overlap or tie:
            best = pocket
            best_overlap = n_overlap
    return best
