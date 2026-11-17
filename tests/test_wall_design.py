"""Wall section strength and wall section-cut forces against hand calculations."""
import math
import textwrap

import pytest

from openstructure import Model
from openstructure.design.section import beta1
from openstructure.design.wall import WallSection, check_walls

E2K = '''
$ CONTROLS
  UNITS  "KN"  "M"  "C"
$ STORIES - IN SEQUENCE FROM TOP
  STORY "S2"  HEIGHT 3
  STORY "S1"  HEIGHT 3
  STORY "Base"  ELEV 0
$ MATERIAL PROPERTIES
  MATERIAL  "C"    TYPE "Concrete"    WEIGHTPERVOLUME 24
  MATERIAL  "C"    SYMTYPE "Isotropic"  E 25000000  U 0.2
  MATERIAL  "C"    FC 30000
$ FRAME SECTIONS
$ WALL PROPERTIES
  SHELLPROP  "W20"  PROPTYPE  "Wall"  MATERIAL "C"  MODELINGTYPE "ShellThin"  WALLTHICKNESS 0.2
$ POINT COORDINATES
  POINT "1"  0 0
  POINT "2"  3 0
$ LINE CONNECTIVITIES
$ AREA CONNECTIVITIES
  AREA "W1"  PANEL  4  "1"  "2"  "2"  "1"  1  1  0  0
$ POINT ASSIGNS
$ LINE ASSIGNS
$ AREA ASSIGNS
  AREAASSIGN  "W1"  "S1"  SECTION "W20"  PIER  "P1"
  AREAASSIGN  "W1"  "S2"  SECTION "W20"  PIER  "P1"
$ LOAD PATTERNS
  LOADPATTERN "DL"  TYPE  "Dead"  SELFWEIGHT  1
'''


def test_planar_wall_section():
    s = WallSection([(0, 0, 3000, 0, 200)], fc=30, fy=400, rho_v=0.0, fiber=5.0)
    assert s.Ag == pytest.approx(600000)
    assert s.Po == pytest.approx(0.85 * 30 * 600000)
    c = 1000.0
    a = beta1(30) * c
    P, Mx, My, phi = s.point(0.0, c)                 # compression at the +x end
    assert P == pytest.approx(0.85 * 30 * 200 * a, rel=0.01)
    assert My == pytest.approx(-P * (1500 - a / 2), rel=0.01)
    assert Mx == pytest.approx(0.0, abs=1.0)
    vn = 0.17 * math.sqrt(30) + 0.0025 * 400
    assert s.shear_capacity("X")[0] == pytest.approx(0.75 * vn * 200 * 0.8 * 3000 / 1e3)
    assert s.shear_capacity("Y")[0] == 0.0


def test_wall_cut_forces(tmp_path):
    p = tmp_path / "wall.e2k"
    p.write_text(textwrap.dedent(E2K))
    m = Model.from_e2k(p)
    m.add_lateral_case("EX", "X", {"S2": 100.0})
    wf = m.wall_forces(["DL", "EX"])
    W = 24 * 0.2 * 3 * 6
    # self-weight lumped at the base nodes goes straight to the supports, so the base cut carries
    # the weight above the first row of elements (0.5 m mesh)
    assert wf["DL"][("P1", "S1")]["P"] == pytest.approx(W - 24 * 0.2 * 3 * 0.25, rel=1e-6)
    assert wf["DL"][("P1", "S2")]["P"] == pytest.approx(W / 2 - 24 * 0.2 * 3 * 0.25, rel=1e-6)
    b, t = wf["EX"][("P1", "S1")], wf["EX"][("P1", "S2")]
    assert abs(b["Vx"]) == pytest.approx(100.0, rel=1e-6)
    assert abs(b["My"]) == pytest.approx(600.0, rel=1e-6)      # 100 kN at 6 m
    assert abs(t["My"]) == pytest.approx(300.0, rel=1e-6)
    assert b["P"] == pytest.approx(0.0, abs=1e-6)
    rows, secs = check_walls(m, {"U": {"DL": 0.9, "EX": 1.0}}, rho_v=0.0025)
    assert len(rows) == 2 and rows[-1]["story"] == "S1"
    assert rows[-1]["DC_Vx"] == pytest.approx(100.0 / secs[("P1", "S1")].shear_capacity("X")[0])
