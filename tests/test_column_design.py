"""Column section strength and moment magnifier against hand calculations (ACI 318-14)."""
import math

import pytest

from openstructure.design.column import magnifier
from openstructure.design.section import RectColumn, beta1, rect_bars


@pytest.fixture()
def sec():
    # 400 x 400, f'c 30, fy 400, 2 + 2 bars of 500 mm2 at 50 mm from the faces
    return RectColumn(D=400, B=400, fc=30, fy=400, bars=rect_bars(400, 400, 2, 2, 500.0, 50.0))


def test_axial_limits(sec):
    Po = 0.85 * 30 * (160000 - 2000) + 400 * 2000
    assert sec.Po == pytest.approx(Po)
    assert sec.phiPn_max == pytest.approx(0.8 * 0.65 * Po)
    assert sec.capacity_ratio(sec.phiPn_max / 1e3, 0, 0) == pytest.approx(1.0, rel=1e-6)
    assert sec.capacity_ratio(-0.9 * 400 * 2000 / 1e3, 0, 0) == pytest.approx(1.0, rel=1e-6)


def test_strain_compatibility_point(sec):
    c = 200.0
    a = beta1(30) * c
    Cc = 0.85 * 30 * a * 400
    Fs_top = (400 - 0.85 * 30) * 1000          # eps = 0.003 (1 - 50/200) = 0.00225 > ey
    Fs_bot = -400 * 1000                       # eps = 0.003 (1 - 350/200) = -0.00225
    P = Cc + Fs_top + Fs_bot
    M = Cc * (200 - a / 2) + Fs_top * 150 + Fs_bot * (-150)
    Pn, M3, M2, phi = sec.point(0.0, c)
    assert Pn == pytest.approx(P, rel=1e-9)
    assert M3 == pytest.approx(M, rel=1e-9)
    assert M2 == pytest.approx(0.0, abs=1e-3)
    assert phi == pytest.approx(0.65 + 0.25 * (0.00225 - 0.002) / 0.003)
    # same strength about the other axis for a doubly symmetric section
    Pn2, M3b, M2b, _ = sec.point(math.pi / 2, c)
    assert (Pn2, M2b) == pytest.approx((P, M))


def test_ray_ratio_uniaxial(sec):
    Pn, M3, M2, phi = sec.point(0.0, 200.0)
    assert sec.capacity_ratio(phi * Pn / 1e3, 0.0, phi * M3 / 1e6) == pytest.approx(1.0, rel=5e-3)
    assert sec.capacity_ratio(0.5 * phi * Pn / 1e3, 0.0, 0.5 * phi * M3 / 1e6) == pytest.approx(0.5, rel=5e-3)


def test_moment_magnifier():
    # 400 mm wide column, lu 3.05 m, Pu 6000 kN, negligible end moments, beta_dns 0.7, Ec 26500 MPa
    Ig = 800 * 400 ** 3 / 12
    m = magnifier(6000e3, 1e6, 2e6, 400, Ig, 26500, 3050, 0.7, sway=True)
    EI = 0.4 * 26500 * Ig / 1.7
    Pc = math.pi ** 2 * EI / 3050 ** 2
    assert m["kl_r"] == pytest.approx(3050 / 120)
    assert m["Pc"] == pytest.approx(Pc)
    assert m["Cm"] == 1.0                                   # M2 below M2,min
    assert m["delta"] == pytest.approx(1 / (1 - 6000e3 / (0.75 * Pc)))
    assert m["Mc"] == pytest.approx(m["delta"] * 6000e3 * (15 + 0.03 * 400))
    # nonsway, double curvature: limit 34 + 12 (M1/M2) = 40 -> not slender
    n = magnifier(6000e3, 300e6, -300e6, 400, Ig, 26500, 3050, 0.7, sway=False)
    assert n["limit"] == 40 and n["delta"] == 1.0
    # single curvature equal moments above M2,min: Cm = 1
    q = magnifier(6000e3, 300e6, 300e6, 400, Ig, 26500, 3050, 0.7, sway=True)
    assert q["Cm"] == pytest.approx(1.0)
