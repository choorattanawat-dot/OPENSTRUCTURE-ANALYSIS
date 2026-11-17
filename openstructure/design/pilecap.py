"""Pile cap check per ACI 318-14 (13.4 deep foundations, 22.5 one-way shear, 22.6 two-way shear).

Rigid cap, piles carry axial load only:   R_i = P / n + Mx y_i / sum(y^2) - My x_i / sum(x^2)
Checks (sectional method; pile reactions act at pile centres, 13.4.2.2):
  1. service pile load (column load + cap self-weight + soil on the cap) against the safe pile load
  2. flexure at the column faces, both directions (13.2.7.1), minimum steel 7.6.1.1
  3. one-way shear at d from the column faces (13.4.2.5 partial count of piles near the section)
  4. two-way shear at d/2 around the column (22.6.5.2)
  5. two-way shear around the most loaded pile, with the perimeter cut by the cap edges
Units: kN, m, kN-m in and out; f'c, fy in MPa. Axes x, y are the cap axes with the column at the origin;
Mx, My are moments about x and y at the top of the cap (right-hand rule), Fz positive down on the cap.

Deep caps, where piles sit closer than d to the column so that the sectional checks see little shear, are
covered by crsi_check() below (CRSI Design Guide for Pile Caps 2015, ch. 5: checks at the column face).
Not covered: strut-and-tie (13.2.6.3), lateral pile loads, pile tension anchorage.
"""
import math
from dataclasses import dataclass, field

S3 = math.sqrt(3) / 2
# CRSI standard patterns (Design Guide for Pile Caps 2015, Fig. 4.4 and Appendix A) in units of the pile spacing L:
# {n: [(|y|, [|x|, ...]), ...]}; every row is mirrored in x and y. They reproduce Ix, Iy of Table 6.1 / Appendix A.
PATTERNS = {
    1: [(0, [0])], 2: [(0, [.5])], 4: [(.5, [.5])], 5: [(0, [0]), (math.sqrt(.5), [math.sqrt(.5)])],
    6: [(.5, [0, 1])], 7: [(0, [0, 1]), (S3, [.5])], 8: [(0, [.5]), (S3, [0, 1])], 9: [(0, [0, 1]), (1, [0, 1])],
    10: [(0, [.5, 1.5]), (S3, [0, 1])], 11: [(0, [0, 1]), (S3, [.5, 1.5])], 12: [(0, [.5, 1.5]), (1, [.5, 1.5])],
    13: [(0, [0, 2 * S3]), (.5, [S3]), (1, [0, 2 * S3])],
    14: [(.5, [.5, 1.5]), (.5 + S3, [0, 1])],
    15: [(0, [0]), (.5, [S3]), (1, [0, 2 * S3]), (1.5, [S3])],
    16: [(.5, [.5, 1.5]), (1.5, [.5, 1.5])],
    17: [(0, [0, 2 * S3]), (.5, [S3]), (1, [0, 2 * S3]), (1.5, [S3])],
    18: [(0, [S3]), (.5, [0, 2 * S3]), (1, [S3]), (1.5, [0, 2 * S3])],
    19: [(0, [0]), (.5, [S3, S3 + 1]), (1, [0]), (1.5, [S3, S3 + 1])],
    20: [(.5, [0, 1, 2]), (1.5, [0, 1, 2])],
    21: [(0, [0, 1]), (1, [0, 1]), (.5, [1 + S3]), (1.5, [1 + S3]), (1 + S3, [.5])],
    22: [(0, [.5, 1.5]), (S3, [0, 1, 2]), (2 * S3, [.5, 1.5])],
    23: [(0, [0, 1, 2]), (1, [0, 1, 2]), (1 + S3, [.5, 1.5])],
    24: [(0, [.5, 1.5]), (S3, [0, 1, 2]), (S3 + 1, [0, 1, 2])],
    26: [(0, [.5 + S3]), (1, [.5 + S3]), (2, [.5 + S3]), (.5, [.5, .5 + 2 * S3]), (1.5, [.5, .5 + 2 * S3])],
    28: [(0, [.5, .5 + 2 * S3]), (1, [.5, .5 + 2 * S3]), (2, [.5, .5 + 2 * S3]), (.5, [.5 + S3]), (1.5, [.5 + S3])],
    30: [(0, [.5, 1.5, 2.5]), (1, [.5, 1.5, 2.5]), (2, [.5, 1.5, 2.5])],
}


def pile_layout(n, spacing, rotate=False):
    """CRSI standard pattern for n piles (2-24, 26, 28, 30), centred on the origin, centre-to-centre distance
    >= spacing, long side A along x. rotate=True puts the long side along y."""
    if n == 3:                                  # triangle, centroid at the origin
        h = math.sqrt(3) / 2 * spacing
        pts = [(-spacing / 2, -h / 3), (spacing / 2, -h / 3), (0.0, 2 * h / 3)]
    else:
        pts = sorted({(sx * x * spacing, sy * y * spacing) for y, xs in PATTERNS[n] for x in xs
                      for sx in (1, -1) for sy in (1, -1)}, key=lambda p: (p[1], p[0]))
        pts = [(x + 0.0, y + 0.0) for x, y in pts]          # no negative zeros
    assert len(pts) == n, (n, len(pts))
    return [(y, x) for x, y in pts] if rotate else pts


@dataclass
class PileCap:
    piles: list                    # [(x, y)] m
    pile_size: float               # side of a square pile or diameter of a round one, m
    h: float                       # cap thickness, m
    col: tuple                     # column dimensions (along x, along y), m
    fc: float = 28.0
    fy: float = 390.0
    edge: float = None             # pile centre to cap edge, m (default: pile size, at least 0.3 m)
    d_offset: float = 0.15         # h - d: pile embedment + cover + bar layers, m
    round_pile: bool = False
    round_col: bool = False        # round column of diameter col[0] (CRSI uses one for 3- and 7-pile caps)
    gamma_c: float = 24.0
    Lx: float = field(init=False)
    Ly: float = field(init=False)

    def __post_init__(self):
        if self.edge is None:
            self.edge = max(self.pile_size, 0.30)
        self.Lx = 2 * (max(abs(x) for x, _ in self.piles) + self.edge)
        self.Ly = 2 * (max(abs(y) for _, y in self.piles) + self.edge)

    @property
    def d(self):
        return self.h - self.d_offset

    @property
    def weight(self):
        return self.Lx * self.Ly * self.h * self.gamma_c

    def reactions(self, P, Mx=0.0, My=0.0):
        """Pile loads (kN, compression +) for column load P down and moments Mx, My about x and y."""
        n = len(self.piles)
        sx = sum(x * x for x, _ in self.piles)
        sy = sum(y * y for _, y in self.piles)
        return [P / n + (Mx * y / sy if sy else 0.0) - (My * x / sx if sx else 0.0) for x, y in self.piles]

    # -- 13.4.2.5: share of a pile reaction that counts as shear on a section `dist` inside the pile centre
    def _share(self, dist):
        return min(max(dist / self.pile_size + 0.5, 0.0), 1.0)

    def check(self, Pu, Mxu=0.0, Myu=0.0):
        """Strength checks for one factored load set (column load only; cap weight is not added because
        it loads the piles without bending or shearing the cap about the column)."""
        R = self.reactions(Pu, Mxu, Myu)
        d, fc, rt = self.d, self.fc, math.sqrt(self.fc)
        out = {"Ru_max": max(R), "Ru_min": min(R), "d": d}
        # flexure and one-way shear, per direction
        for ax, (c, L, b) in {"x": (self.col[0], self.Lx, self.Ly), "y": (self.col[1], self.Ly, self.Lx)}.items():
            k = 0 if ax == "x" else 1
            Mu = Vu = 0.0
            for sgn in (1, -1):
                arm = [(sgn * p[k] - c / 2, r) for p, r in zip(self.piles, R)]
                Mu = max(Mu, sum(r * a for a, r in arm if a > 0))
                Vu = max(Vu, sum(r * self._share(a - d) for a, r in arm))
            Rn = Mu * 1e-3 / (0.9 * b * d * d)                       # MPa
            q = 1 - 2 * Rn / (0.85 * fc)
            rho = 0.85 * fc / self.fy * (1 - math.sqrt(q)) if q > 0 else float("nan")
            rho_min = 0.0020 if self.fy < 420 else 0.0018 * 420 / self.fy
            As, As_min = rho * b * d * 1e4, rho_min * b * self.h * 1e4   # cm2 over the full width b
            phiVc = 0.75 * 0.17 * rt * b * d * 1e3
            out[ax] = dict(Mu=Mu, width=b, As=As, As_min=As_min, As_req=max(As, As_min), Vu=Vu, phiVc=phiVc,
                           DC_shear=Vu / phiVc)
        # two-way shear around the column
        a, bb = (self.col[0] + d) / 2, (self.col[1] + d) / 2
        Vu2 = sum(r * self._share(max(abs(x) - a, abs(y) - bb)) for (x, y), r in zip(self.piles, R))
        bo = 4 * (a + bb)
        beta = max(self.col) / min(self.col)
        vc = min(0.33, 0.17 * (1 + 2 / beta), 0.083 * (40 * d / bo + 2)) * rt
        out["punch_col"] = dict(Vu=Vu2, bo=bo, phiVc=0.75 * vc * bo * d * 1e3, DC=Vu2 / (0.75 * vc * bo * d * 1e3))
        # two-way shear around a single pile (worst ratio), perimeter limited by the cap edges
        worst = None
        for (x, y), r in zip(self.piles, R):
            ex, ey = self.Lx / 2 - abs(x), self.Ly / 2 - abs(y)
            w = self.pile_size + d
            if self.round_pile:
                cand = [math.pi * w, math.pi * w / 2 + 2 * min(ex, ey), math.pi * w / 4 + ex + ey]
            else:
                cand = [4 * w, 2 * (w / 2 + ex) + w, 2 * (w / 2 + ey) + w, (w / 2 + ex) + (w / 2 + ey)]
            bo_p = min(cand)
            vcp = min(0.33, 0.083 * (40 * d / bo_p + 2)) * rt
            phiV = 0.75 * vcp * bo_p * d * 1e3
            if worst is None or r / phiV > worst["DC"]:
                worst = dict(pile=(x, y), Vu=r, bo=bo_p, phiVc=phiV, DC=r / phiV)
        out["punch_pile"] = worst
        out["DC_max"] = max(out["x"]["DC_shear"], out["y"]["DC_shear"], out["punch_col"]["DC"], worst["DC"])
        return out

    def service(self, P, Mx=0.0, My=0.0, soil=0.0):
        """Largest and smallest pile load under service load, including cap weight and `soil` (kN) on it."""
        R = self.reactions(P + self.weight + soil, Mx, My)
        return max(R), min(R)


def design_pilecap(n, pile_size, col, service, factored, spacing=None, h_min=0.6, h_step=0.05, rotate=False, **kw):
    """Smallest thickness (in h_step) of a standard n-pile cap that passes every shear check.
    service = (P, Mx, My), factored = [(name, Pu, Mxu, Myu), ...]. Returns (cap, result dict)."""
    pts = pile_layout(n, spacing or 3 * pile_size, rotate)
    h = h_min
    while True:
        cap = PileCap(pts, pile_size, h, col, **kw)
        res = {name: cap.check(Pu, Mx, My) for name, Pu, Mx, My in factored}
        if max(r["DC_max"] for r in res.values()) <= 1.0 or h > 5.0:
            break
        h = round(h + h_step, 3)
    gov = max(res, key=lambda k: res[k]["DC_max"])
    smax, smin = cap.service(*service)
    return cap, dict(checks=res, governing=gov, pile_service_max=smax, pile_service_min=smin)


def format_pilecap(cap, res):
    r = res["checks"][res["governing"]]
    o = [f"pile cap: {len(cap.piles)} piles {cap.pile_size * 1e3:.0f} mm, cap {cap.Lx:.2f} x {cap.Ly:.2f} x {cap.h:.2f} m "
         f"(d {cap.d:.2f} m), column {cap.col[0]:.2f} x {cap.col[1]:.2f} m, f'c {cap.fc:.1f} MPa, fy {cap.fy:.0f} MPa",
         f"cap weight {cap.weight:.0f} kN; service pile load max {res['pile_service_max']:.0f} kN, min {res['pile_service_min']:.0f} kN",
         f"governing factored case {res['governing']}: pile Ru max {r['Ru_max']:.0f} kN"]
    for ax in ("x", "y"):
        a = r[ax]
        o.append(f"  bars along {ax}: Mu {a['Mu']:.0f} kN-m, As {a['As']:.0f} cm2 (min {a['As_min']:.0f}) over {a['width']:.2f} m; "
                 f"one-way shear Vu {a['Vu']:.0f} / phiVc {a['phiVc']:.0f} kN = {a['DC_shear']:.2f}")
    for k, lab in (("punch_col", "punching at column"), ("punch_pile", "punching at pile")):
        p = r[k]
        o.append(f"  {lab}: Vu {p['Vu']:.0f} / phiVc {p['phiVc']:.0f} kN = {p['DC']:.2f} (bo {p['bo']:.2f} m)")
    return "\n".join(o)


# ---------------------------------------------------------------- CRSI Design Guide for Pile Caps (2015), chapter 5
PSI = 0.083036          # sqrt(f'c in psi) [psi]  =  PSI * sqrt(f'c in MPa) [MPa] / 0.00689476 ... i.e. k sqrt(psi) -> k * PSI sqrt(MPa)


def crsi_check(cap, R, phi_v=0.85, phi_f=0.90, tol=0.075, weight_factor=None, load_factor=1.6):
    """CRSI limit states for a PileCap and factored pile reactions R (kN; one value = the same in every pile).

    CRSI's own safety format is load factor 1.6 on (D + L) with phi 0.85 in shear and 0.90 in bending; with
    ACI 318-14 load combinations use phi_v = 0.75 and weight_factor = 1.2. The factored cap weight outside each
    critical section is deducted (weight_factor defaults to load_factor). tol = adverse pile placement
    tolerance (3 in.), added to every pile distance. Piles are treated as round with the same area.

    Limit states (CRSI numbering):  1 two-way at d/2 from the column,  2/3 one-way at d from the column,
    4 two-way at the column face (when w < d; mandatory when w < d/2),  5/6 one-way at the column face (only when w < d),
    P1 single pile punching, P2 two adjacent piles, P3 corner pile punching, P4 corner pile one-way.
    w = distance from the column face to the nearest pile line, including tol. Direction "x" means a section
    normal to x. For a rectangular column cx and cy are used per direction and (1 + d/c) of limit state 4 becomes
    bo / bs = 1 + 2 d / (cx + cy)  (CRSI itself assumes a square column).
    Returns {name: dict(Vu, phiVc, DC, ...)} plus "flexure" and "DC_max"."""
    n = len(cap.piles)
    R = [float(R)] * n if not hasattr(R, "__len__") else list(R)
    wf = load_factor if weight_factor is None else weight_factor
    d, rt = cap.d, math.sqrt(cap.fc)
    A, B, D = cap.Lx, cap.Ly, cap.h
    q = wf * cap.gamma_c * D                                # factored cap weight per unit area, kN/m2
    c = {"x": cap.col[0], "y": cap.col[1]}
    side = {"x": A, "y": B}
    dp = cap.pile_size if cap.round_pile else cap.pile_size * math.sqrt(4 / math.pi)
    out = {}

    def coord(p, ax):
        return abs(p[0] if ax == "x" else p[1]) + tol

    def w_of(ax):                                           # nearest pile line beyond the column face
        v = [coord(p, ax) - c[ax] / 2 for p in cap.piles if coord(p, ax) - tol > c[ax] / 2 - 1e-9]
        return min(v) if v else None

    # flexure: section at c/4 from the column centre, reinforcement by strain compatibility
    fl = {}
    for ax in ("x", "y"):
        k, width, L_ = (0 if ax == "x" else 1), side["y" if ax == "x" else "x"], side[ax]
        Mu = 0.0
        for sgn in (1, -1):
            m = sum(r * (sgn * p[k] + tol - c[ax] / 4) for p, r in zip(cap.piles, R) if sgn * p[k] > c[ax] / 4)
            Mu = max(Mu, m - q * width * (L_ / 2 - c[ax] / 4) ** 2 / 2)
        Rn = Mu * 1e-3 / (phi_f * width * d * d)
        t = 1 - 2 * Rn / (0.85 * cap.fc)
        As = 0.85 * cap.fc / cap.fy * (1 - math.sqrt(t)) * width * d * 1e4 if t > 0 else float("nan")
        rho = max(1.379 / cap.fy, 3 * PSI * rt / cap.fy)    # 200/fy and 3 sqrt(f'c)/fy in psi
        a1, a2, a3 = rho * width * d * 1e4, 4 / 3 * As, 0.0018 * width * D * 1e4
        req = As if As >= a1 else (a1 if a1 <= a2 else (a2 if a3 <= a2 else min(a3, a1)))
        fl[ax] = dict(Mu=Mu, width=width, As=As, As_req=req)
    beta = max(A, B) / min(A, B)
    short = "x" if A < B else "y"                           # bars spanning the short side
    if beta > 1.0 + 1e-9:                                   # uniform spacing instead of the ACI 13.3.3.3 band:
        f = fl[short]                                       # calculated As x 2 beta / (beta + 1), minimum steel unchanged
        f["As_req"] = max(f["As_req"], f["As"] * 2 * beta / (beta + 1))
    out["flexure"] = fl

    # 1: two-way at d/2 from the column
    a, b = (c["x"] + d) / 2, (c["y"] + d) / 2
    if cap.round_col:
        Vu = sum(r for p, r in zip(cap.piles, R) if math.hypot(*p) + tol > a) - q * (A * B - math.pi * a * a)
        bo, vc = 2 * math.pi * a, min(4.0, 40 * d / (2 * math.pi * a) + 2) * PSI * rt
    else:
        Vu = sum(r for p, r in zip(cap.piles, R) if coord(p, "x") > a or coord(p, "y") > b) - q * (A * B - 4 * a * b)
        bo = 4 * (a + b)
        bc = max(c.values()) / min(c.values())
        vc = min(4.0, 2 + 4 / bc, 40 * d / bo + 2) * PSI * rt
    out["1 two-way d/2"] = dict(Vu=max(Vu, 0.0), phiVc=phi_v * vc * bo * d * 1e3)
    # 2, 3: one-way at d from the column face
    for ax, tag in (("x", "2/3 one-way at d, section normal to x"), ("y", "2/3 one-way at d, section normal to y")):
        k, width, L_ = (0 if ax == "x" else 1), side["y" if ax == "x" else "x"], side[ax]
        xs = c[ax] / 2 + d
        Vu = max(sum(r for p, r in zip(cap.piles, R) if sgn * p[k] + tol > xs) for sgn in (1, -1))
        Vu -= q * width * max(L_ / 2 - xs, 0.0)
        out[tag] = dict(Vu=max(Vu, 0.0), phiVc=phi_v * 2 * PSI * rt * width * d * 1e3)
    # 4: two-way at the column face
    ws = {ax: w_of(ax) for ax in ("x", "y")}
    wv = sorted({round(w, 6) for w in ws.values() if w is not None})
    if cap.round_col:                                       # radial distance to the nearest pile outside the column
        rr = [math.hypot(*p) + tol - c["x"] / 2 for p in cap.piles if math.hypot(*p) > c["x"] / 2]
        wv = [min(rr)] if rr else []
    if wv and min(wv) < d:              # CRSI text: needed when w < d/2; its examples evaluate it whenever w < d
        bs = 2 * (c["x"] + c["y"])
        Vu = sum(r for p, r in zip(cap.piles, R) if coord(p, "x") - tol > c["x"] / 2 or coord(p, "y") - tol > c["y"] / 2)
        Vu -= q * (A * B - c["x"] * c["y"])
        caps = [phi_v * min(d / w * (1 + 4 * d / bs) * 2 * PSI * rt, 32 * PSI * rt) * bs * d * 1e3 for w in wv]
        dc4 = sum(Vu / v for v in caps) / len(caps)         # different w in x and y: average of the two ratios (example 3)
        out["4 two-way at column face"] = dict(Vu=Vu, phiVc=Vu / dc4, phiVc_each=caps, w=ws, required=min(wv) < d / 2)
    # 5, 6: one-way at the column face
    for ax in ("x", "y"):
        w = ws[ax]
        if w is None or w >= d:
            continue
        k, width, L_ = (0 if ax == "x" else 1), side["y" if ax == "x" else "x"], side[ax]
        best = None
        for sgn in (1, -1):
            sel = [(sgn * p[k] + tol - c[ax] / 2, r) for p, r in zip(cap.piles, R) if sgn * p[k] > c[ax] / 2]
            arm = L_ / 2 - c[ax] / 2
            V = sum(r for _, r in sel) - q * width * arm
            M = sum(r * x for x, r in sel) - q * width * arm ** 2 / 2
            if best is None or V > best[0]:
                best = (V, M)
        V, M = best
        mv = min(M / (V * d), 1.0)
        vc = min(d / w * (3.5 - 2.5 * mv) * (1.9 + 0.1 / mv) * PSI * rt, 10 * PSI * rt)
        out[f"5/6 one-way at column face, section normal to {ax}"] = dict(Vu=V, Mu=M, phiVc=phi_v * vc * width * d * 1e3, w=w)
    # piles
    Rmax = max(R)
    v4 = phi_v * 4 * PSI * rt * d * 1e3
    out["P1 pile punching"] = dict(Vu=Rmax, phiVc=v4 * math.pi * (dp + d))
    Lmin = min(math.dist(p, r) for i, p in enumerate(cap.piles) for r in cap.piles[i + 1:]) if n > 1 else None
    if Lmin:
        out["P2 two piles punching"] = dict(Vu=2 * Rmax, phiVc=v4 * (math.pi * (dp + d) + 2 * Lmin))
    out["P3 corner pile punching"] = dict(Vu=Rmax, phiVc=v4 * (math.pi * (dp + d) / 4 + 2 * cap.edge))
    bb = 2 * (cap.edge * math.sqrt(2) + dp / 2 + min(d, 0.330))
    out["P4 corner pile one-way"] = dict(Vu=Rmax, phiVc=phi_v * 2 * PSI * rt * bb * d * 1e3)
    for k, v in out.items():
        if k != "flexure":
            v["DC"] = v["Vu"] / v["phiVc"]
    out["DC_max"] = max(v["DC"] for k, v in out.items() if k != "flexure")
    return out


def format_crsi(res):
    o = [f"{'CRSI limit state':<52}{'Vu kN':>9}{'phiVc kN':>10}{'D/C':>7}"]
    for k, v in res.items():
        if k in ("flexure", "DC_max"):
            continue
        o.append(f"{k:<52}{v['Vu']:>9.0f}{v['phiVc']:>10.0f}{v['DC']:>7.2f}")
    for ax, f in res["flexure"].items():
        o.append(f"flexure, bars along {ax}: Mu {f['Mu']:.0f} kN-m, As {f['As']:.0f} cm2, required {f['As_req']:.0f} cm2 over {f['width']:.2f} m")
    return "\n".join(o)
