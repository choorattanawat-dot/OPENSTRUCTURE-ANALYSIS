"""Unit tests that need no project model. The full ETABS comparison is in test_reference.py."""
import textwrap

import pytest

from openstructure import e2k
from openstructure.etabs_tables import read_tables
from openstructure.units import factors, G_ACC

TEMPLATE = '''
$ CONTROLS
  UNITS  "{F}"  "{L}"  "C"
$ STORIES - IN SEQUENCE FROM TOP
  STORY "S2"  HEIGHT {h}
  STORY "S1"  HEIGHT {h}
  STORY "Base"  ELEV 0
$ MATERIAL PROPERTIES
  MATERIAL  "C"    TYPE "Concrete"    WEIGHTPERVOLUME {w}
  MATERIAL  "C"    SYMTYPE "Isotropic"  E {E}  U 0.2
$ FRAME SECTIONS
  FRAMESECTION  "COL"  MATERIAL "C"  SHAPE "Concrete Rectangular"  D {d} B {b}
$ SLAB PROPERTIES
  SHELLPROP  "S20"  PROPTYPE  "Slab"  MATERIAL "C"  MODELINGTYPE "ShellThin"  SLABTYPE "Slab"  SLABTHICKNESS {t}
$ POINT COORDINATES
  POINT "1"  0 0
  POINT "2"  {x} 0
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
$ SHELL UNIFORM LOAD SETS
  SHELLUNIFORMLOADSET "A"  LOADPAT "DL"  VALUE {q}
'''


def _parse(tmp_path, name, **kw):
    p = tmp_path / name
    p.write_text(textwrap.dedent(TEMPLATE.format(**kw)))
    return e2k.parse_e2k(p)


def test_unit_factors():
    assert factors(["N", "mm", "C"]) == (1.0, 1.0)
    assert factors(["KGF", "M", "C"]) == (G_ACC, 1000.0)
    assert factors(["KN", "M", "C"]) == (1000.0, 1000.0)
    with pytest.raises(ValueError):
        factors(["STONE", "M", "C"])


def test_same_model_in_three_unit_systems(tmp_path):
    """E = 25 GPa, gamma = 24 kN/m3, column 400x800, slab 200, q = 2 kPa, storey 3 m."""
    g = G_ACC
    a = _parse(tmp_path, "nmm.e2k", F="N", L="MM", h=3000, w=24e-6, E=25000, d=800, b=400, t=200, x=6000, q=0.002)
    b = _parse(tmp_path, "kgfm.e2k", F="KGF", L="M", h=3, w=24000 / g, E=25e9 / g, d=0.8, b=0.4, t=0.2, x=6, q=2000 / g)
    c = _parse(tmp_path, "knm.e2k", F="KN", L="M", h=3, w=24, E=25e6, d=0.8, b=0.4, t=0.2, x=6, q=2)
    for m in (a, b, c):
        assert m["elev"]["S2"] == pytest.approx(6000)
        assert m["mats"]["C"]["E"] == pytest.approx(25000)
        assert m["mats"]["C"]["WEIGHTPERVOLUME"] == pytest.approx(24e-6)
        assert m["fsec"]["COL"]["D"] == pytest.approx(800)
        assert m["sprop"]["S20"]["SLABTHICKNESS"] == pytest.approx(200)
        assert m["pts"]["2"][0] == pytest.approx(6000)
        assert m["lsets"]["A"]["DL"] == pytest.approx(0.002)
        assert m["patterns"]["DL"] == ("Dead", 1.0)


def test_table_reader_continuation(tmp_path):
    p = tmp_path / "t.txt"
    p.write_text('TABLE:  "BASE REACTIONS"\n   "Output Case"=DL   "Case Type"=LinStatic   FX=1 _\n        FZ=2.5\n'
                 'TABLE:  "OTHER"\n   A=1\n')
    t = read_tables(p)
    assert t["BASE REACTIONS"] == [{"Output Case": "DL", "Case Type": "LinStatic", "FX": 1.0, "FZ": 2.5}]
    assert t["OTHER"] == [{"A": 1.0}]
