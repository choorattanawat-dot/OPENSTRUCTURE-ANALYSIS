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

Not covered: strut-and-tie for deep caps (13.2.6.3), where most piles sit closer than d to the column
the sectional one-way check sees almost no shear; lateral pile loads; pile tension anchorage.
"""
import math
from dataclasses import dataclass, field

ROWS = {1: [1], 2: [2], 3: [2, 1], 4: [2, 2], 5: [2, 1, 2], 6: [3, 3], 7: [2, 3, 2], 8: [3, 2, 3],
        9: [3, 3, 3], 10: [3, 4, 3], 11: [4, 3, 4], 12: [4, 4, 4]}


def pile_layout(n, spacing, rotate=False):
    """Standard pattern for n piles, centre-to-centre distance >= spacing, centred on the origin.
    Rows run along x; rows with different counts are staggered (row distance 0.866 spacing).
    rotate=True turns the pattern so that the rows run along y."""
    rows = ROWS[n]
    ys = [0.0]
    for a, b in zip(rows, rows[1:]):
        ys.append(ys[-1] + spacing * (1.0 if a == b else math.sqrt(3) / 2))
    pts = [((i - (m - 1) / 2) * spacing, y) for m, y in zip(rows, ys) for i in range(m)]
    cx = sum(p[0] for p in pts) / n
    cy = sum(p[1] for p in pts) / n
    pts = [(x - cx, y - cy) for x, y in pts]
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
