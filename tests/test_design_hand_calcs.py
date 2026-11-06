"""Design checks against hand calculations.

Sources for the factors are documented in docs/CODE_ASSUMPTIONS.md and in the
code-plugin assumptions. The arithmetic in each test is written out so it can
be repeated on a calculator. It does not call the library until the comparison.
"""

import math

import pytest

from openstructure.design import ACI318, EIT1008, check_column_uniaxial, design_beam_flexure, design_beam_shear
from openstructure.design.section import nominal_at_c, pure_compression
from openstructure.units import mpa_to_ksc, n_to_kn

# Bar areas and the section used for the column hand calculation.
B_MM = 400.0
H_MM = 400.0
COVER = 60.0
BAR_AREA = 500.0  # mm², a round number so the arithmetic stays readable
BARS = [
    (BAR_AREA, 140.0),
    (BAR_AREA, 140.0),
    (BAR_AREA, -140.0),
    (BAR_AREA, -140.0),
]
FC = 28.0
FY = 400.0


def test_aci_beta1_and_eit_beta1_differ_at_28_mpa():
    # ACI 318-19: β1 = 0.85 for f'c ≤ 28 MPa.
    assert ACI318().beta1(28.0) == 0.85
    assert ACI318().beta1(35.0) == pytest.approx(0.85 - 0.05 * (35.0 - 28.0) / 7.0)
    # EIT textbook rounding: reduction starts above 280 ksc (27.46 MPa).
    fc_ksc = mpa_to_ksc(28.0)
    expected = 0.85 - 0.05 * (fc_ksc - 280.0) / 70.0
    assert EIT1008().beta1(28.0) == pytest.approx(expected)
    assert EIT1008().beta1(28.0) < 0.85


def test_aci_singly_reinforced_beam_hand_calculation():
    """Mu = 200 kN·m, b = 300, d = 450, h = 500, f'c = 28, fy = 400.

    Rn = Mu / (φ b d²) with φ = 0.90.
    ρ = (0.85 f'c / fy) [1 − √(1 − 2 Rn / (0.85 f'c))].
    As,min = max(0.25√f'c, 1.4) bw d / fy, which does not govern.
    """
    mu_nmm = 200.0e6
    phi = 0.90
    b_mm, d_mm = 300.0, 450.0
    rn = mu_nmm / (phi * b_mm * d_mm**2)
    rho = (0.85 * FC / FY) * (1.0 - math.sqrt(1.0 - 2.0 * rn / (0.85 * FC)))
    as_hand = rho * b_mm * d_mm
    as_min = max(0.25 * math.sqrt(FC) / FY, 1.4 / FY) * b_mm * d_mm
    assert as_hand > as_min

    result = design_beam_flexure(200.0, b_mm, d_mm, 500.0, FC, FY, ACI318())
    assert result.status == "OK"
    assert result.doubly_reinforced is False
    assert result.as_mm2 == pytest.approx(as_hand, rel=1e-9)
    assert result.utilization == pytest.approx(1.0, rel=1e-3)
    assert "0.90" in result.report.to_markdown()


def test_large_moment_is_doubly_reinforced_for_aci_and_checked():
    # ACI c/d limit is 0.375. φMn at that limit is below 400 kN·m for this beam
    # (see the limit step in the report). EIT allows c/d = 0.75 * 0.6 = 0.45.
    aci = design_beam_flexure(400.0, 300.0, 450.0, 500.0, FC, FY, ACI318(), d_prime_mm=60.0)
    eit = design_beam_flexure(400.0, 300.0, 450.0, 500.0, FC, FY, EIT1008(), d_prime_mm=60.0)
    assert aci.doubly_reinforced is True
    assert aci.as_prime_mm2 > 0
    assert aci.status == "OK"
    assert aci.utilization == pytest.approx(1.0, rel=1e-3)
    # The EIT limit is higher, so this Mu may still be singly reinforced.
    assert eit.as_mm2 > 0
    assert eit.status == "OK"


def test_aci_shear_hand_calculation():
    """Vu = 200 kN, bw = 300, d = 450, f'c = 28, 2-leg 9 mm stirrups, fyt = 400.

    Vc = 0.17 √f'c bw d (ACI 318-19 §22.5.5.1 (a), Nu = 0).
    φ = 0.75. Vs,req = Vu/φ − Vc.
    """
    vc = 0.17 * math.sqrt(FC) * 300.0 * 450.0
    phi = 0.75
    vs_req = 200.0e3 / phi - vc
    av = 2.0 * math.pi * 9.0**2 / 4.0
    s_strength = av * 400.0 * 450.0 / vs_req
    result = design_beam_shear(200.0, 300.0, 450.0, FC, 400.0, ACI318(), stirrup_diameter_mm=9.0, legs=2)
    assert result.stirrups_required is True
    assert result.status == "OK"
    assert result.s_provided_mm <= s_strength + 1e-6
    assert result.s_provided_mm % 25 == 0
    # Provided spacing is the strength spacing rounded down, and it must satisfy Vu.
    vs_prov = av * 400.0 * 450.0 / result.s_provided_mm
    phi_vn = phi * (vc + vs_prov) / 1000.0
    assert result.phi_vn_kN == pytest.approx(phi_vn, rel=1e-9)
    assert result.utilization <= 1.0


def test_eit_shear_uses_0_53_and_phi_0_85():
    fc_ksc = mpa_to_ksc(FC)
    vc_kgf = 0.53 * math.sqrt(fc_ksc) * 30.0 * 45.0
    vc_n = vc_kgf * 9.80665
    code = EIT1008()
    assert code.vc_n(FC, 300.0, 450.0, True) == pytest.approx(vc_n, rel=1e-12)
    assert code.factors().phi_shear == 0.85
    result = design_beam_shear(80.0, 300.0, 450.0, FC, 240.0, code, stirrup_diameter_mm=9.0, legs=2)
    assert "0.85" in result.report.to_markdown()
    assert result.report.steps[1].result.startswith("Vc")


def test_column_po_and_phi_pn_max_hand_calculation():
    """Four 500 mm² bars in a 400 mm square tied column, f'c = 28, fy = 400.

    Po = 0.85 f'c (Ag − Ast) + fy Ast.
    ACI tied: φPn,max = 0.65 × 0.80 × Po.
    EIT tied: φPn,max = 0.70 × 0.80 × Po.
    """
    ag = 400.0 * 400.0
    ast = 4.0 * BAR_AREA
    po = 0.85 * FC * (ag - ast) + FY * ast
    point = pure_compression(400.0, 400.0, BARS, FC, FY)
    assert point.pn_n == pytest.approx(po)
    aci = check_column_uniaxial(0.0, 0.0, 400.0, 400.0, BARS, FC, FY, ACI318(), tied=True)
    eit = check_column_uniaxial(0.0, 0.0, 400.0, 400.0, BARS, FC, FY, EIT1008(), tied=True)
    assert aci.po_kN == pytest.approx(po / 1000.0)
    assert aci.phi_pn_max_kN == pytest.approx(0.65 * 0.80 * po / 1000.0)
    assert eit.phi_pn_max_kN == pytest.approx(0.70 * 0.80 * po / 1000.0)
    assert eit.phi_pn_max_kN > aci.phi_pn_max_kN


def test_column_utilization_scales_with_axial_demand():
    aci = ACI318()
    cap = check_column_uniaxial(0.0, 0.0, 400, 400, BARS, FC, FY, aci).phi_pn_max_kN
    half = check_column_uniaxial(0.5 * cap, 0.0, 400, 400, BARS, FC, FY, aci)
    over = check_column_uniaxial(1.2 * cap, 0.0, 400, 400, BARS, FC, FY, aci)
    assert half.status == "OK"
    assert half.utilization == pytest.approx(0.5, rel=0.02)
    assert over.status == "NG"
    assert over.utilization == pytest.approx(1.2, rel=0.02)


def test_balanced_point_matches_explicit_whitney_arithmetic():
    """Compression on the +y face, εt = εy = 0.002 at the bars y = -140 mm.

    d_t = 200 − (−140) = 340 mm.
    c = 0.003 / (0.003 + 0.002) × 340 = 204 mm.
    β1 = 0.85, a = 173.4 mm.
    Both compression bars (y = +140, 60 mm from the face) are inside the block
    and yielded. Both tension bars are at εt = εy, so fs = fy.
    """
    b_mm, h_mm = 400.0, 400.0
    c_mm = 0.003 / 0.005 * 340.0
    a_mm = 0.85 * c_mm
    cc = 0.85 * FC * b_mm * a_mm
    # Compression bars: fs' = min(fy, Es * 0.003 * (c - 60) / c), then subtract 0.85 f'c.
    eps_s = 0.003 * (c_mm - 60.0) / c_mm
    fs_comp = min(FY, eps_s * 200_000.0)
    force_comp_each = (fs_comp - 0.85 * FC) * BAR_AREA
    force_tens_each = -FY * BAR_AREA
    pn = cc + 2 * force_comp_each + 2 * force_tens_each
    # Moment about the centroid. Concrete lever = h/2 - a/2. Bar lever = y.
    mn = cc * (h_mm / 2.0 - a_mm / 2.0) + 2 * force_comp_each * 140.0 + 2 * force_tens_each * (-140.0)
    point = nominal_at_c(b_mm, h_mm, BARS, c_mm, FC, FY, 0.85, compression_sign=1)
    assert point.pn_n == pytest.approx(pn, rel=1e-9)
    assert point.mn_nmm == pytest.approx(mn, rel=1e-9)
    assert point.eps_t == pytest.approx(0.002, rel=1e-9)
    assert n_to_kn(pn) > 0
