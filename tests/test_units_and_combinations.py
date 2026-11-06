"""Unit conversions and the named load-combination presets."""

import pytest

from openstructure.combinations import apply_preset
from openstructure.model import Model
from openstructure.units import KSC_PER_MPA, kn_per_m_to_n_per_mm, mpa_to_ksc


def test_kn_per_m_matches_n_per_mm():
    assert kn_per_m_to_n_per_mm(12.5) == 12.5


def test_ksc_conversion_uses_standard_gravity():
    # 1 kgf/cm² = 0.0980665 MPa, so 1 MPa is about 10.197 ksc.
    assert mpa_to_ksc(1.0) == pytest.approx(KSC_PER_MPA)
    assert mpa_to_ksc(0.0980665) == pytest.approx(1.0)


def test_aci_gravity_factors_when_only_dead_and_live():
    model = Model()
    model.add_material("C", ec_mpa=25000, fc_mpa=28, fy_mpa=400)
    model.add_rectangular_section("S", "C", 300, 500)
    model.add_node("A", 0, 0)
    model.add_node("B", 4, 0)
    model.add_support("A", pinned=True)
    model.add_support("B", roller_x=True)
    model.add_member("AB", "A", "B", "S")
    model.add_load_case("DL", "dead")
    model.add_load_case("LL", "live")
    names = apply_preset(model, "ACI318-19")
    assert "ACI-1.4D" in names
    assert "ACI-1.2D+1.6L" in names
    assert model.combinations["ACI-1.4D"].factors == {"DL": 1.4}
    assert model.combinations["ACI-1.2D+1.6L"].factors == {"DL": 1.2, "LL": 1.6}
    assert all("W" not in name for name in names)


def test_aci_wind_is_reversible_and_cites_5_3_1():
    model = Model()
    model.add_material("C", ec_mpa=25000, fc_mpa=28, fy_mpa=400)
    model.add_rectangular_section("S", "C", 300, 500)
    model.add_node("A", 0, 0)
    model.add_node("B", 4, 0)
    model.add_support("A", fixed=True)
    model.add_support("B", fixed=True)
    model.add_member("AB", "A", "B", "S")
    model.add_load_case("DL", "dead")
    model.add_load_case("LL", "live")
    model.add_load_case("WX", "wind")
    names = apply_preset(model, "ACI318-19")
    assert any(name.endswith("[+WX]") for name in names)
    assert any(name.endswith("[-WX]") for name in names)
    positive = next(name for name in names if "0.9D+1.0W" in name and name.endswith("[+WX]"))
    negative = next(name for name in names if "0.9D+1.0W" in name and name.endswith("[-WX]"))
    assert model.combinations[positive].factors["WX"] == pytest.approx(1.0)
    assert model.combinations[negative].factors["WX"] == pytest.approx(-1.0)
    assert "5.3.1" in model.combinations[positive].reference


def test_eit_gravity_is_1_4_and_1_7_and_wind_is_flagged():
    model = Model()
    model.add_material("C", ec_mpa=25000, fc_mpa=28, fy_mpa=400)
    model.add_rectangular_section("S", "C", 300, 500)
    model.add_node("A", 0, 0)
    model.add_node("B", 4, 0)
    model.add_support("A", pinned=True)
    model.add_support("B", roller_x=True)
    model.add_member("AB", "A", "B", "S")
    model.add_load_case("DL", "dead")
    model.add_load_case("LL", "live")
    model.add_load_case("WX", "wind")
    apply_preset(model, "EIT1008")
    assert model.combinations["EIT-1.4D+1.7L"].factors == {"DL": 1.4, "LL": 1.7}
    wind = next(combo for name, combo in model.combinations.items() if "1.7W" in name and "[+WX]" in name)
    assert wind.factors["DL"] == pytest.approx(0.75 * 1.4)
    assert wind.factors["LL"] == pytest.approx(0.75 * 1.7)
    assert wind.factors["WX"] == pytest.approx(0.75 * 1.7)
    assert "ASSUMPTION" in wind.reference
