"""Equilibrium diagrams, independent of OpenSees.

The end actions below are the joint-on-member forces for a fixed beam with a
downward uniform load, in the sign convention checked against OpenSees
localForce (see tests/test_analysis.py). Section forces must recover the
classical beam results.
"""

import pytest

from openstructure.analysis import EndActions, LocalPointLoad, section_force_at

# w = 10 kN/m, L = 4 m. wL/2 = 20 kN, wL^2/12 = 13.333 kN·m, wL^2/24 = 6.667 kN·m.
FIXED = EndActions(
    fx_i_kN=0.0,
    fy_i_kN=20.0,
    m_i_kNm=40.0 / 3.0,
    fx_j_kN=0.0,
    fy_j_kN=20.0,
    m_j_kNm=-40.0 / 3.0,
)


def test_fixed_beam_udl_moments_and_shears():
    length = 4.0
    wy = -10.0
    at_0 = section_force_at(0.0, length, FIXED, 0.0, wy, [])
    at_mid = section_force_at(2.0, length, FIXED, 0.0, wy, [])
    at_end = section_force_at(4.0, length, FIXED, 0.0, wy, [])
    assert at_0.v_kN == pytest.approx(20.0)
    assert at_end.v_kN == pytest.approx(-20.0)
    assert at_0.m_kNm == pytest.approx(-10.0 * 16.0 / 12.0)
    assert at_end.m_kNm == pytest.approx(-10.0 * 16.0 / 12.0)
    assert at_mid.m_kNm == pytest.approx(10.0 * 16.0 / 24.0)
    assert at_mid.n_kN == pytest.approx(0.0)


def test_simple_beam_point_load_midspan():
    # Reactions 15 kN up, no end moments, point load -30 kN at midspan, L = 6 m.
    end = EndActions(0.0, 15.0, 0.0, 0.0, 15.0, 0.0)
    load = [LocalPointLoad(0.5, 0.0, -30.0)]
    before = section_force_at(3.0 - 1e-6, 6.0, end, 0.0, 0.0, load)
    after = section_force_at(3.0 + 1e-6, 6.0, end, 0.0, 0.0, load)
    mid = section_force_at(3.0, 6.0, end, 0.0, 0.0, load)
    assert before.v_kN == pytest.approx(15.0)
    assert after.v_kN == pytest.approx(-15.0)
    # Moment is continuous and equal to P*L/4 = 45 kN·m.
    assert mid.m_kNm == pytest.approx(45.0, rel=1e-6)
    assert before.m_kNm == pytest.approx(45.0, rel=1e-4)
