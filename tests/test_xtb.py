import pytest

from pocket_atlas.chem.xtb_strain import (
    XtbUnavailable,
    find_xtb,
    parse_xtb_energy,
    strain_kcal,
)


def test_strain_kcal_conversion():
    assert strain_kcal(0.01, 0.0) == pytest.approx(6.27509)


def test_parse_xtb_energy():
    text = "something\n          TOTAL ENERGY              -12.345678 Eh\n"
    assert parse_xtb_energy(text) == pytest.approx(-12.345678)


def test_find_xtb_fails_closed():
    from pocket_atlas.chem import xtb_strain as mod

    original = mod.shutil.which
    mod.shutil.which = lambda _name: None
    try:
        with pytest.raises(XtbUnavailable, match="xtb"):
            find_xtb()
    finally:
        mod.shutil.which = original
