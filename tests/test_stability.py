"""Cantilever column checked against hand formulas (linear drift and P-Delta amplification)."""
import math
import textwrap

import pytest

from openstructure import Model
from openstructure.stability import story_drift, stability

E2K = '''
$ CONTROLS
  UNITS  "KN"  "M"  "C"
$ STORIES - IN SEQUENCE FROM TOP
  STORY "S1"  HEIGHT 4
  STORY "Base"  ELEV 0
$ MATERIAL PROPERTIES
  MATERIAL  "C"    TYPE "Concrete"    WEIGHTPERVOLUME 24
  MATERIAL  "C"    SYMTYPE "Isotropic"  E 25000000  U 0.2
$ FRAME SECTIONS
  FRAMESECTION  "COL"  MATERIAL "C"  SHAPE "Concrete Rectangular"  D 0.4 B 0.4
$ POINT COORDINATES
  POINT "1"  0 0
$ LINE CONNECTIVITIES
  LINE  "C1"  COLUMN  "1"  "1"  1
$ AREA CONNECTIVITIES
$ POINT ASSIGNS
  POINTASSIGN  "1"  "Base"  RESTRAINT "UX UY UZ RX RY RZ"
$ LINE ASSIGNS
  LINEASSIGN  "C1"  "S1"  SECTION "COL"  ANG  0
$ AREA ASSIGNS
$ LOAD PATTERNS
  LOADPATTERN "DL"  TYPE  "Dead"  SELFWEIGHT  1
'''
L, EI = 4.0, 25e6 * 0.4 ** 4 / 12


@pytest.fixture()
def model(tmp_path):
    p = tmp_path / "col.e2k"
    p.write_text(textwrap.dedent(E2K))
    return Model.from_e2k(p)


def test_linear_cantilever_drift(model):
    r = story_drift(model, "X", total_kN=10.0)[0]
    assert r["V_kN"] == pytest.approx(10.0)
    assert r["drift_cm_mm"] == pytest.approx(10.0 * L ** 3 / (3 * EI) * 1e3, rel=1e-6)


def test_pdelta_amplification(model):
    # self-weight lumped at the top node = half the column weight; scale it up to P = 500 kN
    w_top = 24 * 0.4 * 0.4 * L / 2
    P = 500.0
    res = stability(model, gravity_factors={"DL": P / w_top}, total_kN=10.0)
    r = res["X"][0]
    assert r["P_kN"] == pytest.approx(P)
    d0 = 10.0 * L ** 3 / (3 * EI)
    assert r["Q"] == pytest.approx(P * d0 / (10.0 * L), rel=1e-6)
    # PDelta transformation (one element) = linearised geometric stiffness P/L: amp = 1 / (1 - P L^2 / (3 EI))
    assert r["amp_pdelta"] == pytest.approx(1 / (1 - P * L ** 2 / (3 * EI)), rel=1e-4)
    assert r["amp_pdelta"] == pytest.approx(r["amp_1_over_1mQ"], rel=1e-4)


def test_base_reaction_includes_load_on_support(model):
    from openstructure.reactions import base_reactions
    res = base_reactions(model, {"U": {"DL": 1.4}})
    W = 24 * 0.4 * 0.4 * L
    assert res["DL"]["1"]["Fz"] == pytest.approx(W)          # full column weight, not only the top half
    assert res["U"]["1"]["Fz"] == pytest.approx(1.4 * W)
    assert (res["DL"]["1"]["x"], res["DL"]["1"]["y"]) == (0.0, 0.0)
