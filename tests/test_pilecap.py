import math

from openstructure.design.pilecap import PileCap, crsi_check, design_pilecap, pile_layout


def test_layouts_keep_spacing():
    for n in list(range(2, 14)) + [16]:
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


def test_layouts_match_crsi_table_6_1():
    """Pile group moments of inertia Ix, Iy / L^2 of CRSI Design Guide for Pile Caps, Table 6.1."""
    ref = {2: (0, 0.5), 3: (0.5, 0.5), 4: (1, 1), 5: (2, 2), 6: (1.5, 4), 7: (3, 3), 8: (4.5, 4.5), 9: (6, 6),
           10: (4.5, 9), 11: (6, 12), 12: (8, 15), 13: (7, 21), 16: (20, 20)}
    for n, (ix, iy) in ref.items():
        p = pile_layout(n, 1.0)
        assert abs(sum(y * y for _, y in p) - ix) < 1e-9 and abs(sum(x * x for x, _ in p) - iy) < 1e-9, n


def test_crsi_example_1():
    """CRSI Design Guide for Pile Caps, design example 1 (16 piles of 40 tons, 48 in. cap), kips -> kN."""
    ft, inch, kip = 0.3048, 0.0254, 4.4482216
    cap = PileCap(pile_layout(16, 3 * ft), pile_size=8 * inch, h=48 * inch, col=(22 * inch, 22 * inch), fc=3000 * 0.00689476,
                  fy=413.7, edge=15 * inch, d_offset=10 * inch, round_pile=True, gamma_c=0.150 * kip / ft ** 3)
    r = crsi_check(cap, 1.6 * 80 * kip, phi_v=0.85, tol=3 * inch)
    book = {"1 two-way d/2": (1433, 1698), "2/3 one-way at d, section normal to x": (494, 488),
            "4 two-way at column face": (1924, 3226), "5/6 one-way at column face, section normal to x": (971, 2440),
            "P1 pile punching": (128, 1022), "P2 two piles punching": (256, 1530), "P3 corner pile punching": (128, 468),
            "P4 corner pile one-way": (128, 271)}
    for k, (vu, pv) in book.items():
        assert abs(r[k]["Vu"] / kip - vu) < 0.01 * vu, (k, r[k]["Vu"] / kip)
        assert abs(r[k]["phiVc"] / kip - pv) < 0.01 * pv, (k, r[k]["phiVc"] / kip)
    f = r["flexure"]["x"]
    assert abs(f["Mu"] / kip / ft - 2704) < 0.01 * 2704
    assert abs(f["As_req"] / 6.4516 - 17.5) < 0.2                       # rho_min b d governs: 17.5 in2
