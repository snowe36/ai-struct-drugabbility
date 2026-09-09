from pocket_atlas.cases import load_case
from pocket_atlas.dynamics.labels import load_relaxdb_bundle, relaxdb_residues, resolve_prior


def test_relaxdb_kras_is_official_cpmg():
    residues = relaxdb_residues("kras_switch2")
    assert residues is not None
    assert len(residues) >= 50
    assert 12 in residues  # P-loop, absent from the literature switch subset
    assert 60 in residues


def test_relaxdb_does_not_map_blac_onto_tem1():
    assert relaxdb_residues("tem1_horn") is None
    bundle = load_relaxdb_bundle()
    blac = bundle["entries"]["BLAC_CPMG"]
    assert blac["uniprot"] == "P9WKD3"
    assert blac["sequence"].startswith("GDLADR")
    assert not blac["sequence"].startswith("HPETL")


def test_prior_literature_fallback_for_tem1():
    case = load_case("tem1_horn")
    residues, source = resolve_prior(case, prior="relaxdb")
    assert source == "nmr_literature"
    assert 166 in residues
    assert residues == set(case.nmr_residues)


def test_prior_relaxdb_for_kras():
    case = load_case("kras_switch2")
    residues, source = resolve_prior(case, prior="relaxdb")
    assert source == "relaxdb_cpmg"
    assert residues != set(case.nmr_residues)
    assert 12 in residues
