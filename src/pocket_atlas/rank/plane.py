"""Experimental exchange on ligand-scale clearance maxima."""

from __future__ import annotations

import json
from pathlib import Path

import yaml

from pocket_atlas.dynamics.align import chain_ca_sequence, map_cpmg_index_to_resseq
from pocket_atlas.dynamics.labels import RESOURCES, cpmg_entry
from pocket_atlas.io.rcsb import load_structure
from pocket_atlas.paths import PROCESSED, REPORTS, ensure_dirs
from pocket_atlas.pockets import detect_pockets
from pocket_atlas.prepare import prepare_structure, write_pdb
from pocket_atlas.rank import ApoRanking, rank_apo, rank_case

PANEL = RESOURCES.parent / "cpmg_panel.yaml"

SHORT_LABELS = {
    "tem1_horn": "TEM-1",
    "aqadk": "AdK",
    "blac_mtb": "BlaC",
    "blvrb": "BLVRB",
    "cypa": "CypA",
    "vhr": "VHR",
    "argkin": "ArgKin",
    "rnase_a": "RNase A",
    "chey": "CheY",
    "htra2_pdz": "HTRA2",
    "kras_switch2": "KRAS",
}


def load_panel() -> dict:
    with PANEL.open() as handle:
        return yaml.safe_load(handle)


def _rank_cpmg_crystal(spec: dict) -> ApoRanking:
    entry = cpmg_entry(spec["entry_id"])
    structure = prepare_structure(load_structure(spec["pdb_id"]), chain=spec["chain"])
    write_pdb(structure, PROCESSED / f"{spec['id']}_{spec['pdb_id']}.pdb")
    pockets = detect_pockets(structure, chain=spec["chain"], ligand_mode="exclude")
    protein = set(structure.residue_numbers(chain=spec["chain"]))
    resseqs, pdb_seq = chain_ca_sequence(structure, spec["chain"])
    index_to_resseq = map_cpmg_index_to_resseq(entry["sequence"], resseqs, pdb_seq)
    prior = {index_to_resseq[i] for i in entry["exchange_residues"] if i in index_to_resseq}
    mapped_frac = len(index_to_resseq) / max(len(entry["sequence"]), 1)
    ranking = rank_apo(
        pockets,
        protein_residues=protein,
        prior_residues=prior,
        prior_source="relaxdb_cpmg",
        case=spec["id"],
        tag="crystal",
        pdb_id=spec["pdb_id"],
    )
    ranking.extra.update(
        {
            "title": spec.get("title"),
            "entry_id": spec["entry_id"],
            "n_mapped": len(index_to_resseq),
            "mapped_frac": round(mapped_frac, 3),
            "n_prior_mapped": len(prior),
        }
    )
    return ranking


def rank_panel_protein(spec: dict) -> ApoRanking:
    if spec.get("existing_case"):
        prior = spec.get("prior", "relaxdb")
        ranking = rank_case(spec["id"], prior=prior)
        ranking.extra["title"] = spec.get("title")
        ranking.extra["entry_id"] = spec.get("entry_id")
        return ranking
    return _rank_cpmg_crystal(spec)


def plane_row_stats(ranking: ApoRanking) -> dict:
    scored = ranking.scored
    if not scored:
        return {
            "max_enrichment": 0.0,
            "dscore_at_max_enr": 0.0,
            "max_dscore": 0.0,
            "enrichment_at_max_dscore": 0.0,
            "n_enriched": 0,
        }
    top_enr = max(scored, key=lambda s: s.enrichment)
    top_geo = max(scored, key=lambda s: s.pocket.dscore)
    return {
        "max_enrichment": top_enr.enrichment,
        "dscore_at_max_enr": top_enr.pocket.dscore,
        "max_dscore": top_geo.pocket.dscore,
        "enrichment_at_max_dscore": top_geo.enrichment,
        "n_enriched": sum(1 for s in scored if s.enrichment > 1.0),
    }


def write_plane_report(rankings: list[ApoRanking], path: Path | None = None) -> Path:
    ensure_dirs()
    path = path or (REPORTS / "plane.md")
    lines = [
        "# Exchange × clearance",
        "",
        "Every apo cavity: lining μs–ms exchange enrichment vs clearance-maxima `dscore`.",
        "Dashed enrichment = 1 is protein background. YAML linings are not used.",
        "",
        "| Protein | PDB | Prior | Cavities | Mapped | Max enrichment | dscore at max enr | n(enr>1) | Enrichment at max dscore |",
        "|---------|-----|-------|----------|--------|----------------|-------------------|----------|--------------------------|",
    ]
    for ranking in rankings:
        stats = plane_row_stats(ranking)
        mapped = ranking.extra.get("mapped_frac")
        mapped_s = "" if mapped is None else f"{mapped:.3f}"
        label = SHORT_LABELS.get(ranking.case, ranking.case)
        lines.append(
            f"| {label} | `{ranking.pdb_id}` | {ranking.n_prior} | {len(ranking.scored)} | "
            f"{mapped_s} | {stats['max_enrichment']:.2f} | {stats['dscore_at_max_enr']:.3f} | "
            f"{stats['n_enriched']} | {stats['enrichment_at_max_dscore']:.2f} |"
        )
    lines.append("")
    path.write_text("\n".join(lines) + "\n")
    return path


def run_exchange_plane() -> list[ApoRanking]:
    ensure_dirs()
    panel = load_panel()
    rows: list[ApoRanking] = []
    for spec in list(panel.get("literature_extra") or []) + list(panel.get("proteins") or []):
        print(f"plane {spec['id']} {spec['pdb_id']} …", flush=True)
        rows.append(rank_panel_protein(spec))
    payload = {"proteins": [r.as_dict() for r in rows]}
    (PROCESSED / "cpmg_plane.json").write_text(json.dumps(payload, indent=2))
    return rows
