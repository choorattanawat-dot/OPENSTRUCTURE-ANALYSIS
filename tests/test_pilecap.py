import math

from openstructure.design.pilecap import PileCap, crsi_check, design_pilecap, pile_layout


def test_layouts_keep_spacing():
    for n in list(range(2, 25)) + [26, 28, 30]:
        p = pile_layout(n, 1.0)
        assert len(p) == n
        dmin = min(math.dist(a, b) for i, a in enumerate(p) for b in p[i + 1:])
        assert dmin > 0.999, n
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
    """Pile group inertias Ix, Iy / L^2 and largest edge pile force of CRSI Design Guide for Pile Caps, Tables 6.1
    and 6.2 (21 piles: Iy = 8 sqrt(3) + 21 = 34.86 as derived in its Appendix A; Table 6.1 prints 34.46)."""
    ref = {2: (0, 0.5), 3: (0.5, 0.5), 4: (1, 1), 5: (2, 2), 6: (1.5, 4), 7: (3, 3), 8: (4.5, 4.5), 9: (6, 6),
           10: (4.5, 9), 11: (6, 12), 12: (8, 15), 13: (7, 21), 14: (13.2, 14), 15: (16, 18), 16: (20, 20), 17: (16, 24),
           18: (19, 28.5), 19: (22, 33.86), 20: (25, 40), 21: (29.93, 34.86), 22: (31.5, 35), 23: (37.86, 40),
           24: (42.32, 45), 26: (40, 60.52), 28: (50, 67.25), 30: (60, 87.5)}
    for n, (ix, iy) in ref.items():
        p = pile_layout(n, 1.0)
        assert abs(sum(y * y for _, y in p) - ix) < 0.006 and abs(sum(x * x for x, _ in p) - iy) < 0.006, n
    edge = {4: (0.5, 0.5), 5: (0.35, 0.35), 6: (0.33, 0.25), 9: (0.17, 0.17), 10: (0.19, 0.17), 12: (0.13, 0.10),
            16: (0.075, 0.075), 20: (0.06, 0.05), 24: (0.044, 0.044), 30: (0.033, 0.029)}
    for n, (fx, fy) in edge.items():
        cap = PileCap(pile_layout(n, 1.0), 0.3, 1.0, (0.5, 0.5))
        assert abs(max(cap.reactions(0.0, Mx=1.0)) - fx) < 0.006, n
        assert abs(max(cap.reactions(0.0, My=1.0)) - fy) < 0.006, n


def _crsi(n, tons, dp, E, D, c, L=3, round_col=False):
    ft, inch, kip = 0.3048, 0.0254, 4.4482216
    cap = PileCap(pile_layout(n, L * ft), pile_size=dp * inch, h=D * inch, col=(c * inch, c * inch), fc=3000 * 0.00689476,
                  fy=413.7, edge=E * inch, d_offset=10 * inch, round_pile=True, round_col=round_col,
                  gamma_c=0.150 * kip / ft ** 3)
    return crsi_check(cap, 1.6 * 2 * tons * kip, phi_v=0.85, tol=3 * inch), kip, ft


def _compare(r, kip, book, tol=0.012):
    for k, (vu, pv) in book.items():
        assert abs(r[k]["Vu"] / kip - vu) <= tol * vu + 1, (k, r[k]["Vu"] / kip)
        assert abs(r[k]["phiVc"] / kip - pv) <= tol * pv, (k, r[k]["phiVc"] / kip)


def test_crsi_examples_2_to_5():
    """Design examples 2-5 of the CRSI Design Guide for Pile Caps (kips), every limit state the book evaluates."""
    X, Y = "5/6 one-way at column face, section normal to x", "5/6 one-way at column face, section normal to y"
    r, kip, ft = _crsi(5, 100, 10, 21, 43, 20)                               # example 2
    _compare(r, kip, {"1 two-way d/2": (1245, 1303), "4 two-way at column face": (1231, 1160), X: (620, 1111),
                      "P1 pile punching": (320, 830), "P2 two piles punching": (640, 1272),
                      "P3 corner pile punching": (320, 466), "P4 corner pile one-way": (320, 293)})
    assert abs(r["flexure"]["x"]["Mu"] / kip / ft - 1211) < 12 and abs(r["flexure"]["x"]["As_req"] / 6.4516 - 10.22) < 0.1
    r, kip, ft = _crsi(6, 100, 10, 21, 48, 22)                               # example 3: w differs in x and y
    _compare(r, kip, {"1 two-way d/2": (1245, 1698), X: (616, 627), Y: (938, 2018), "P1 pile punching": (320, 1069),
                      "P2 two piles punching": (640, 1578), "P3 corner pile punching": (320, 564),
                      "P4 corner pile one-way": (320, 338)})
    assert abs(r["4 two-way at column face"]["DC"] - 1.09) < 0.015           # average of 0.578 and 1.62
    assert abs(r["flexure"]["x"]["As_req"] / 6.4516 - 10.6) < 0.1 and abs(r["flexure"]["y"]["As_req"] / 6.4516 - 9.85) < 0.1
    r, kip, ft = _crsi(7, 40, 8, 15, 38, 17, round_col=True)                 # example 4: round column
    _compare(r, kip, {"1 two-way d/2": (725, 737), X: (363, 1038), Y: (236, 337), "P1 pile punching": (128, 589),
                      "P2 two piles punching": (256, 965)})
    assert "4 two-way at column face" not in r                               # radial w = 30.5 in. > d = 28 in.
    r, kip, ft = _crsi(5, 400, 20, 36, 69, 40, L=5)                          # example 5: high load piles
    _compare(r, kip, {"4 two-way at column face": (4900, 5055), X: (2472, 4300), "P1 pile punching": (1280, 2725),
                      "P2 two piles punching": (2560, 4043), "P3 corner pile punching": (1280, 1472),
                      "P4 corner pile one-way": (1280, 812)})
    assert r["1 two-way d/2"]["Vu"] == 0.0                                   # no pile outside the d/2 section
    assert abs(r["flexure"]["x"]["Mu"] / kip / ft - 7264) < 73 and abs(r["flexure"]["x"]["As_req"] / 6.4516 - 30.87) < 0.3


def test_crsi_example_6_lateral():
    """Example 6: 16 piles at 3 ft, 48 kips of pile capacity left for bending -> My = 1,920 k-ft (Table 6.2: 0.075 My / L)."""
    cap = PileCap(pile_layout(16, 3.0), 0.67, 4.0, (1.83, 1.83))
    assert abs(48 / max(cap.reactions(0.0, My=1.0)) - 1920) < 1


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
