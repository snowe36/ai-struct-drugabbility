from pocket_atlas.cases import list_cases, load_case


def test_cases_load():
    names = list_cases()
    assert "tem1_horn" in names
    assert "kras_switch2" in names
    assert "vp35_iid" in names
    tem = load_case("tem1_horn")
    assert 70 in tem.residue_set("catalytic_site")
    assert tem.cryptic_residues
    assert tem.nmr_residues
    kras = load_case("kras_switch2")
    assert kras.residue_set("nucleotide_site")
