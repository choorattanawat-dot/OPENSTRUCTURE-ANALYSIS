import math

from openstructure.design.pilecap import PileCap, design_pilecap, pile_layout


def test_layouts_keep_spacing():
    for n in range(2, 13):
        p = pile_layout(n, 1.0)
        assert len(p) == n
        dmin = min(math.dist(a, b) for i, a in enumerate(p) for b in p[i + 1:])
        assert dmin > 0.999
        assert abs(sum(x for x, _ in p)) < 1e-9 and abs(sum(y for _, y in p)) < 1e-9


def test_four_pile_hand_calc():
    cap = PileCap(pile_layout(4, 1.2), pile_size=0.4, h=0.9, col=(0.5, 0.5), fc=28, fy=420)
    assert cap.reactions(4000) == [1000.0] * 4
    R = cap.reactions(4000, Mx=240)                      # sum y^2 = 4 * 0.6^2 = 1.44 -> +-100 kN
    assert sorted(round(r) for r in R) == [900, 900, 1100, 1100]
    r = cap.check(4000)
    assert abs(r["x"]["Mu"] - 2 * 1000 * (0.6 - 0.25)) < 1e-6          # 700 kN-m at the column face
    # d = 0.75 m: piles 0.35 m from the face are inside the section at d -> no one-way shear
    assert r["x"]["Vu"] == 0.0
    # column critical section half-size (0.5 + 0.75) / 2 = 0.625: pile centre 0.025 inside -> share 0.4375
    assert abs(r["punch_col"]["Vu"] - 4 * 1000 * (0.5 - 0.025 / 0.4)) < 1e-6
    # corner pile: two sides of (0.4 + 0.75) / 2 + 0.4 edge
    assert abs(r["punch_pile"]["bo"] - 2 * (0.575 + 0.4)) < 1e-9


def test_design_finds_thickness():
    cap, res = design_pilecap(10, 0.35, (0.4, 0.8), (4800, 0, 0), [("U", 6200, 0, 0)], fc=24, fy=390, rotate=True)
    assert len(cap.piles) == 10 and cap.Ly > cap.Lx
    assert res["checks"]["U"]["DC_max"] <= 1.0
    thinner = PileCap(cap.piles, 0.35, cap.h - 0.05, (0.4, 0.8), fc=24, fy=390)
    assert thinner.check(6200)["DC_max"] > 1.0
