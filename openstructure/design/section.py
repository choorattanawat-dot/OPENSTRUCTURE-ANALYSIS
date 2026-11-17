"""Rectangular RC column section: factored P-M2-M3 interaction surface (ACI 318-14 ch. 22, 21.2).

Axes: y along ETABS local 2 (section depth D), z along ETABS local 3 (width B), origin at centroid.
M3 = moment about axis 3 (bends over the depth D), M2 = moment about axis 2 (bends over the width B).
Compression positive. Units N, mm, MPa; the surface is returned in kN and kN-m.
"""
import math
from dataclasses import dataclass, field

import numpy as np

ES = 200000.0
ECU = 0.003


def beta1(fc):
    return min(0.85, max(0.65, 0.85 - 0.05 * (fc - 28.0) / 7.0))


def rect_bars(D, B, n_face3, n_face2, bar_area, cover_to_center):
    """Bars on the perimeter. n_face3 = bars in each face parallel to axis 3 (spread over B),
    n_face2 = bars in each face parallel to axis 2 (spread over D); corner bars counted in both."""
    ys = np.linspace(-D / 2 + cover_to_center, D / 2 - cover_to_center, n_face2)
    zs = np.linspace(-B / 2 + cover_to_center, B / 2 - cover_to_center, n_face3)
    bars = [(y, z) for y in (ys[0], ys[-1]) for z in zs]
    bars += [(y, z) for y in ys[1:-1] for z in (zs[0], zs[-1])]
    return [(y, z, bar_area) for y, z in bars]


def _clip(poly, u, dmin):
    """Part of a convex polygon where p.u >= dmin (Sutherland-Hodgman, one half-plane)."""
    out = []
    n = len(poly)
    for i in range(n):
        a, b = poly[i], poly[(i + 1) % n]
        da, db = a[0] * u[0] + a[1] * u[1] - dmin, b[0] * u[0] + b[1] * u[1] - dmin
        if da >= 0:
            out.append(a)
        if (da >= 0) != (db >= 0):
            t = da / (da - db)
            out.append((a[0] + t * (b[0] - a[0]), a[1] + t * (b[1] - a[1])))
    return out


def _area_centroid(poly):
    if len(poly) < 3:
        return 0.0, 0.0, 0.0
    A = cx = cy = 0.0
    for i in range(len(poly)):
        x0, y0 = poly[i]
        x1, y1 = poly[(i + 1) % len(poly)]
        w = x0 * y1 - x1 * y0
        A += w
        cx += (x0 + x1) * w
        cy += (y0 + y1) * w
    A /= 2
    if abs(A) < 1e-9:
        return 0.0, 0.0, 0.0
    return abs(A), cx / (6 * A), cy / (6 * A)


@dataclass
class RectColumn:
    D: float                 # depth along local 2, mm
    B: float                 # width along local 3, mm
    fc: float                # MPa
    fy: float                # MPa
    bars: list               # [(y, z, area)]
    tied: bool = True
    _surf: object = field(default=None, repr=False)

    @property
    def Ag(self):
        return self.D * self.B

    @property
    def diag(self):
        return math.hypot(self.D, self.B)

    @property
    def Ast(self):
        return sum(b[2] for b in self.bars)

    @property
    def rho(self):
        return self.Ast / self.Ag

    @property
    def Po(self):
        """Nominal axial strength at zero eccentricity, N (22.4.2.2)."""
        return 0.85 * self.fc * (self.Ag - self.Ast) + self.fy * self.Ast

    @property
    def phi_c(self):
        return 0.65 if self.tied else 0.75

    @property
    def phiPn_max(self):
        """Design axial strength limit, N (22.4.2.1): 0.80 (tied) or 0.85 (spiral) of phi Po."""
        return (0.80 if self.tied else 0.85) * self.phi_c * self.Po

    def point(self, theta, c):
        """Nominal strength and phi for neutral-axis angle theta (direction of the compression face,
        measured from +y towards +z) and neutral-axis depth c. Returns (Pn N, M3 N-mm, M2 N-mm, phi)."""
        u = (math.cos(theta), math.sin(theta))
        corners = [(-self.D / 2, -self.B / 2), (self.D / 2, -self.B / 2), (self.D / 2, self.B / 2), (-self.D / 2, self.B / 2)]
        dmax = max(p[0] * u[0] + p[1] * u[1] for p in corners)
        dlow = min(p[0] * u[0] + p[1] * u[1] for p in corners)
        a = min(beta1(self.fc) * c, dmax - dlow)
        A, cy, cz = _area_centroid(_clip(corners, u, dmax - a))
        Cc = 0.85 * self.fc * A
        P, M3, M2 = Cc, Cc * cy, Cc * cz
        et = 0.0
        for y, z, As in self.bars:
            d = dmax - (y * u[0] + z * u[1])           # depth from the extreme compression fibre
            eps = ECU * (1 - d / c)
            fs = max(-self.fy, min(self.fy, eps * ES))
            if d <= a and fs > 0:
                fs -= 0.85 * self.fc                   # concrete displaced by the bar
            F = fs * As
            P += F
            M3 += F * y
            M2 += F * z
            et = max(et, -eps)
        ety = self.fy / ES
        lo = self.phi_c
        phi = lo if et <= ety else 0.9 if et >= 0.005 else lo + (0.9 - lo) * (et - ety) / (0.005 - ety)
        return P, M3, M2, phi

    def surface(self, n_theta=36, n_c=48):
        """Design surface points (phi Pn kN, phi M2 kN-m, phi M3 kN-m), shape (n_theta, n_c, 3),
        with phi Pn capped at phiPn_max. Cached."""
        if self._surf is None or self._surf.shape[:2] != (n_theta, n_c):
            diag = self.diag
            cs = np.concatenate([np.geomspace(20 * diag, 1.2 * diag, 8, endpoint=False),
                                 np.linspace(1.2 * diag, 0.02 * diag, n_c - 8)])
            S = np.zeros((n_theta, n_c, 3))
            for i in range(n_theta):
                th = 2 * math.pi * i / n_theta
                for j, c in enumerate(cs):
                    P, M3, M2, phi = self.point(th, c)
                    S[i, j] = (min(phi * P, self.phiPn_max) / 1e3, phi * M2 / 1e6, phi * M3 / 1e6)
            self._surf = S
        return self._surf

    def _triangles(self):
        S = self.surface()
        nt, nc, _ = S.shape
        top = np.array([self.phiPn_max / 1e3, 0.0, 0.0])
        bot = np.array([-0.9 * self.fy * self.Ast / 1e3, 0.0, 0.0])
        tri = []
        for i in range(nt):
            k = (i + 1) % nt
            tri.append((top, S[i, 0], S[k, 0]))
            tri.append((bot, S[i, -1], S[k, -1]))
            for j in range(nc - 1):
                tri.append((S[i, j], S[k, j], S[k, j + 1]))
                tri.append((S[i, j], S[k, j + 1], S[i, j + 1]))
        return np.array(tri)

    def capacity_ratio(self, Pu, M2u, M3u):
        """Demand/capacity along the ray from the origin to (Pu kN, M2u, M3u kN-m). > 1 fails."""
        if not hasattr(self, "_tri") or self._tri is None:
            self._tri = self._triangles()
        T = self._tri
        v = np.array([Pu, M2u, M3u], float)
        n = np.linalg.norm(v)
        if n < 1e-9:
            return 0.0
        dvec = v / n
        e1, e2 = T[:, 1] - T[:, 0], T[:, 2] - T[:, 0]
        h = np.cross(dvec, e2)
        a = np.einsum("ij,ij->i", e1, h)
        ok = np.abs(a) > 1e-12
        f = np.where(ok, 1.0 / np.where(ok, a, 1.0), 0.0)
        s = -T[:, 0]
        uu = f * np.einsum("ij,ij->i", s, h)
        q = np.cross(s, e1)
        vv = f * (q @ dvec)
        t = f * np.einsum("ij,ij->i", e2, q)
        hit = ok & (uu >= -1e-9) & (vv >= -1e-9) & (uu + vv <= 1 + 1e-9) & (t > 1e-9)
        if not hit.any():
            return float("inf")
        return float(n / t[hit].min())

    def uniaxial(self, axis=3, n_c=60):
        """Design P-M curve about one axis: list of (phi Pn kN, phi Mn kN-m), compression at +face."""
        th = 0.0 if axis == 3 else math.pi / 2
        h = self.D if axis == 3 else self.B
        out = []
        for c in np.concatenate([np.geomspace(20 * h, 1.2 * h, 8, endpoint=False), np.linspace(1.2 * h, 0.02 * h, n_c)]):
            P, M3, M2, phi = self.point(th, c)
            out.append((min(phi * P, self.phiPn_max) / 1e3, phi * (M3 if axis == 3 else M2) / 1e6))
        return out
