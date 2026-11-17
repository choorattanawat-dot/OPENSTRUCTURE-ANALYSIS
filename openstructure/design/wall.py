"""RC wall pier / lift core check per ACI 318-14 (ordinary structural walls, ch. 11 and 22).

Section: the horizontal cut through all wall panels of a pier at the bottom of a story (any plan
shape: box, C, L, planar), meshed into fibres along the wall centre lines. Vertical reinforcement
is smeared: As = rho_v * t per unit length, at the centre line.
  * axial + biaxial bending: strain compatibility (ecu 0.003, rectangular stress block 0.85 f'c over
    a = beta1 c), phi 0.65-0.90 from the extreme steel strain, phi Pn <= 0.80 phi Po;
    demand/capacity along the ray to (Pu, Mx, My) taken about the centroid of the cut
  * in-plane shear per direction (11.5.4): webs within 45 deg of the force direction,
    Vn = (0.17 sqrt(f'c) + rho_t fy) * sum(t * 0.8 lw) <= 0.83 sqrt(f'c) * sum(t * 0.8 lw), phi 0.75
Corner overlaps of the centre-line model are not deducted (area is slightly overestimated).
Not checked here: boundary elements, slenderness of the wall, lintel / coupling beams, torsion.
"""
import math

import numpy as np

from .section import ECU, ES, RectColumn, beta1


class WallSection:
    tied = True
    phi_c = 0.65

    def __init__(self, segs, fc, fy, rho_v=0.0025, rho_h=0.0025, xc=None, yc=None, fiber=50.0, name="wall"):
        """segs: [(x1, y1, x2, y2, t)] in mm (global plan coordinates)."""
        self.segs, self.fc, self.fy, self.rho_v, self.rho_h, self.name = segs, fc, fy, rho_v, rho_h, name
        pts, areas = [], []
        for x1, y1, x2, y2, t in segs:
            L = math.hypot(x2 - x1, y2 - y1)
            n = max(1, round(L / fiber))
            for i in range(n):
                f = (i + 0.5) / n
                pts.append((x1 + f * (x2 - x1), y1 + f * (y2 - y1)))
                areas.append(t * L / n)
        self.xy = np.array(pts)
        self.A = np.array(areas)
        self.Ag = float(self.A.sum())
        c = (self.xy * self.A[:, None]).sum(0) / self.Ag
        self.xc, self.yc = (c[0] if xc is None else xc), (c[1] if yc is None else yc)
        self.r = self.xy - np.array([self.xc, self.yc])
        self.As = self.A * rho_v
        self.Ast = float(self.As.sum())
        self.diag = float(np.ptp(self.xy[:, 0]) + np.ptp(self.xy[:, 1]))
        self._surf = None
        self._tri = None

    @property
    def Po(self):
        return 0.85 * self.fc * (self.Ag - self.Ast) + self.fy * self.Ast

    @property
    def phiPn_max(self):
        return 0.80 * self.phi_c * self.Po

    def point(self, theta, c):
        """Returns (Pn N, moment about the global X axis N-mm, moment about global Y N-mm, phi) for
        compression towards direction theta (from +x to +y). The slots match RectColumn.point:
        index 1 is reported as 'M3' = Mx, index 2 as 'M2' = My."""
        u = np.array([math.cos(theta), math.sin(theta)])
        d = self.r @ u
        depth = d.max() - d                         # from the extreme compression fibre
        a = beta1(self.fc) * c
        eps = ECU * (1 - depth / c)
        fs = np.clip(eps * ES, -self.fy, self.fy)
        inblock = depth <= a
        Fc = np.where(inblock, 0.85 * self.fc * (self.A - self.As), 0.0)
        F = Fc + fs * self.As
        P = F.sum()
        Mx = (F * self.r[:, 1]).sum()               # +Mx: compression on the +y side
        My = -(F * self.r[:, 0]).sum()              # +My: compression on the -x side
        et = max(0.0, -eps.min())
        ety = self.fy / ES
        phi = 0.65 if et <= ety else 0.9 if et >= 0.005 else 0.65 + 0.25 * (et - ety) / (0.005 - ety)
        return P, Mx, My, phi

    surface = RectColumn.surface
    _triangles = RectColumn._triangles
    _ray = RectColumn.capacity_ratio

    def capacity_ratio(self, Pu, Mx, My):
        """Demand/capacity for Pu (kN, compression +), Mx, My (kN-m, about the cut centroid)."""
        return self._ray(Pu, My, Mx)                # surface slots are (P, 'M2' = My, 'M3' = Mx)

    def shear_capacity(self, direction):
        """phi Vn (kN) for shear along 'X' or 'Y', and the web area used (mm2)."""
        Acv = 0.0
        for x1, y1, x2, y2, t in self.segs:
            dx, dy = abs(x2 - x1), abs(y2 - y1)
            if (dx >= dy) == (direction.upper() == "X"):
                Acv += t * 0.8 * math.hypot(dx, dy)
        vn = min(0.17 * math.sqrt(self.fc) + self.rho_h * self.fy, 0.83 * math.sqrt(self.fc))
        return 0.75 * vn * Acv / 1e3, Acv


def check_walls(model, combos=None, fc=None, fy=390.0, rho_v=0.0025, rho_h=0.0025):
    """Check every wall pier cut for every combination. fc MPa (default: wall material), fy MPa.
    rho_v, rho_h: vertical / horizontal reinforcement ratios (the e2k holds no wall reinforcement).
    Returns (rows sorted by story, sections)."""
    combos = combos or model.M["combos"]
    pats = sorted({p for c in combos.values() for p in c})
    wf = model.wall_forces(pats)
    cuts = model.wall_cuts()
    rows, secs = [], {}
    for key, c in cuts.items():
        mat = model.M["mats"][model.M["sprop"][c["sec"]]["mat"]]
        s = WallSection(c["segs"], fc or mat["FC"], fy, rho_v, rho_h, xc=c["xc"], yc=c["yc"], name=f"{key[0]} @ {key[1]}")
        secs[key] = s
        vx_cap, vy_cap = s.shear_capacity("X")[0], s.shear_capacity("Y")[0]
        best = None
        for cname, fac in combos.items():
            f = {q: sum(k * wf[p][key][q] for p, k in fac.items() if p in wf) for q in ("P", "Vx", "Vy", "Mx", "My", "T")}
            dc = s.capacity_ratio(f["P"], f["Mx"], f["My"])
            row = dict(pier=key[0], story=key[1], z_m=c["z"] / 1e3, combo=cname, Pu_kN=f["P"], Mx_kNm=f["Mx"], My_kNm=f["My"],
                       Vx_kN=f["Vx"], Vy_kN=f["Vy"], T_kNm=f["T"], DC_PMM=dc,
                       DC_Vx=abs(f["Vx"]) / vx_cap if vx_cap else float("inf"),
                       DC_Vy=abs(f["Vy"]) / vy_cap if vy_cap else float("inf"),
                       Pu_over_phiPnmax=f["P"] / (s.phiPn_max / 1e3))
            row["DC"] = max(row["DC_PMM"], row["DC_Vx"], row["DC_Vy"])
            if best is None or row["DC"] > best["DC"]:
                best = row
        rows.append(best)
    rows.sort(key=lambda r: (r["pier"], -r["z_m"]))
    return rows, secs


def format_walls(rows, secs):
    s = next(iter(secs.values()))
    o = [f"f'c {s.fc:.1f} MPa, fy {s.fy:.0f} MPa, rho_v {s.rho_v:.4f}, rho_h {s.rho_h:.4f} (assumed unless given), "
         f"Ag {s.Ag / 1e6:.3f} m2, phiPn,max {s.phiPn_max / 1e3:.0f} kN, phiVn X {s.shear_capacity('X')[0]:.0f} / Y {s.shear_capacity('Y')[0]:.0f} kN",
         f"\n{'pier':<6}{'story':<14}{'combo':<12}{'Pu kN':>8}{'Mx':>8}{'My':>8}{'Vx':>7}{'Vy':>7}{'PMM':>7}{'Vx':>6}{'Vy':>6}"]
    for r in rows:
        o.append(f"{r['pier']:<6}{r['story']:<14}{r['combo']:<12}{r['Pu_kN']:>8.0f}{r['Mx_kNm']:>8.0f}{r['My_kNm']:>8.0f}"
                 f"{r['Vx_kN']:>7.0f}{r['Vy_kN']:>7.0f}{r['DC_PMM']:>7.2f}{r['DC_Vx']:>6.2f}{r['DC_Vy']:>6.2f}"
                 f"{'  FAIL' if r['DC'] > 1 else ''}")
    return "\n".join(o)
