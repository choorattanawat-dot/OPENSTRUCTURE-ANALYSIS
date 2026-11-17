"""DPT 1301/1302-61 (มยผ. 1301/1302-61) seismic loads for sites outside the Bangkok basin.

Implemented: site coefficients Fa/Fv (Tables 1.4-2, 1.4-3), SDS/SD1 (eq 1.4-1..4), design spectra
(Figs 1.4-1..4), importance factor (Table 1.5-1), seismic design category (Tables 1.6-1, 1.6-2 with
the 0.8 Ts exception), approximate period (eq 3.3-1/2), equivalent lateral force (ch. 3),
design drift (eq 3.7-1, Table 2.11-1) and P-Delta coefficient theta (eq 3.8-1/2).
Not implemented yet: Bangkok basin spectra (Tables 1.4-4..7), accidental torsion, RSA (ch. 4).

Units: g for spectral values, kN, m, s.
"""
from dataclasses import dataclass

# Table 1.4-1: Ss, S1 (g) on rock. Add sites as they are read from the standard.
SITES = {"ขอนแก่น อ.เมือง": (0.053, 0.030), "khon kaen mueang": (0.053, 0.030)}

_FA_SS = (0.25, 0.5, 0.75, 1.0, 1.25)
FA = {"A": (0.8,) * 5, "B": (1.0,) * 5, "C": (1.2, 1.2, 1.1, 1.0, 1.0),
      "D": (1.6, 1.4, 1.2, 1.1, 1.0), "E": (2.5, 1.7, 1.2, 0.9, 0.9)}
_FV_S1 = (0.1, 0.2, 0.3, 0.4, 0.5)
FV = {"A": (0.8,) * 5, "B": (1.0,) * 5, "C": (1.7, 1.6, 1.5, 1.4, 1.3),
      "D": (2.4, 2.0, 1.8, 1.6, 1.5), "E": (3.5, 3.2, 2.8, 2.4, 2.4)}
IMPORTANCE = {"I": 1.0, "II": 1.0, "III": 1.25, "IV": 1.5}
SDC_ORDER = "กขคง"

# Table 2.3-1 (subset for RC buildings): R, Omega0, Cd
SYSTEMS = {
    "bearing_wall_ordinary_rc_wall": (4, 2.5, 4), "bearing_wall_special_rc_wall": (5, 2.5, 5),
    "bf_ordinary_rc_wall": (5, 2.5, 4.5), "bf_special_rc_wall": (6, 2.5, 5),
    "mf_ordinary_rc": (3, 3, 2.5), "mf_intermediate_rc": (5, 3, 4.5), "mf_special_rc": (8, 3, 5.5),
    "dual_imf_ordinary_rc_wall": (5.5, 2.5, 4.5), "dual_imf_special_rc_wall": (6.5, 2.5, 5),
    "dual_smf_special_rc_wall": (7, 2.5, 5.5), "interactive_ordinary_wall_frame": (4.5, 2.5, 4),
}
# Table 2.11-1, "other structures": allowable story drift ratio by occupancy category
DRIFT_LIMIT = {"I": 0.020, "II": 0.020, "III": 0.015, "IV": 0.010}


def _interp(x, xs, ys):
    if x <= xs[0]:
        return ys[0]
    if x >= xs[-1]:
        return ys[-1]
    for i in range(len(xs) - 1):
        if xs[i] <= x <= xs[i + 1]:
            return ys[i] + (ys[i + 1] - ys[i]) * (x - xs[i]) / (xs[i + 1] - xs[i])


@dataclass
class Site:
    """Design spectrum parameters for a site outside the Bangkok basin."""
    Ss: float
    S1: float
    soil: str = "D"       # section 1.4.2: class D when there is no soil data
    occupancy: str = "II"

    @classmethod
    def named(cls, name, **kw):
        Ss, S1 = SITES[name.lower() if name.lower() in SITES else name]
        return cls(Ss, S1, **kw)

    @property
    def Fa(self):
        return _interp(self.Ss, _FA_SS, FA[self.soil])

    @property
    def Fv(self):
        return _interp(self.S1, _FV_S1, FV[self.soil])

    @property
    def SDS(self):
        return 2 / 3 * self.Fa * self.Ss

    @property
    def SD1(self):
        return 2 / 3 * self.Fv * self.S1

    @property
    def Ts(self):
        return self.SD1 / self.SDS

    @property
    def I(self):
        return IMPORTANCE[self.occupancy]

    def sa_elf(self, T):
        """Figs 1.4-1 (SD1 <= SDS) and 1.4-2 (SD1 > SDS), 5 % damping."""
        SDS, SD1 = self.SDS, self.SD1
        if SD1 <= SDS:
            return SDS if T <= self.Ts else SD1 / T
        if T <= 0.2:
            return SDS
        if T <= 1.0:
            return SDS + (SD1 - SDS) * (T - 0.2) / 0.8
        return SD1 / T

    def sa_dynamic(self, T):
        """Figs 1.4-3 (SD1 <= SDS) and 1.4-4 (SD1 > SDS), 5 % damping."""
        SDS, SD1 = self.SDS, self.SD1
        if SD1 <= SDS:
            T0 = 0.2 * SD1 / SDS
            if T < T0:
                return SDS * (0.4 + 0.6 * T / T0)
            return SDS if T <= self.Ts else SD1 / T
        if T < 0.2:
            return SDS * (0.4 + 0.6 * T / 0.2)
        if T <= 1.0:
            return SDS + (SD1 - SDS) * (T - 0.2) / 0.8
        return SD1 / T

    def sdc(self, Ta=None):
        """Seismic design category. Outside the basin: the more severe of Tables 1.6-1 and 1.6-2;
        if Ta < 0.8 Ts, Table 1.6-1 alone may be used (section 1.6)."""
        col = {"I": 0, "II": 0, "III": 1, "IV": 2}[self.occupancy]
        cat = [["ก", "ก", "ก"], ["ข", "ข", "ค"], ["ค", "ค", "ง"], ["ง", "ง", "ง"]]

        def one(v, lim):
            return cat[sum(v >= x for x in lim)][col]

        a = one(self.SDS, (0.167, 0.33, 0.50))
        b = one(self.SD1, (0.067, 0.133, 0.20))
        if Ta is not None and Ta < 0.8 * self.Ts:
            return a
        return max(a, b, key=SDC_ORDER.index)


def period_approx(H, material="rc"):
    """Eq 3.3-1 / 3.3-2: T = 0.02 H (RC) or 0.03 H (steel), H in m."""
    return (0.02 if material == "rc" else 0.03) * H


def k_exponent(T):
    return 1.0 if T <= 0.5 else 2.0 if T >= 2.5 else 1 + (T - 0.5) / 2


def story_weights(model, dead=1.0, live=0.0):
    """Effective seismic weight per story (kN): dead-type patterns * dead + live-type * live
    (section 2.8.2: all dead load, 25 % of storage live load, partitions, equipment)."""
    from .stability import default_gravity_factors
    gf = default_gravity_factors(model, dead=dead, live=live)
    coords = model.B["coords"]
    out = {s: 0.0 for s in model.B["stories"]}
    for pat, f in gf.items():
        if not f:
            continue
        for n, w in model.loads.get(pat, {}).items():
            if n not in model.fixed:
                out[model.story_of(coords[n - 1][2])] += f * w
    base = model.B["stories"][0]
    out.pop(base, None)
    return {s: w for s, w in out.items() if w > 0}


def elf(model, site, system=None, R=None, Cd=None, T_model=None, material="rc", live=0.0, beta=1.0):
    """Equivalent lateral force procedure (ch. 3) on an openstructure Model, both directions.

    system: key of SYSTEMS, or give R and Cd directly.
    T_model: {"X": T, "Y": T} from modal analysis; the period used is min(T_model, 1.5 Ta).
    Returns base shear, story forces, design drift vs the limit and theta per story.
    """
    from .stability import story_drift
    if system:
        R, _, Cd = SYSTEMS[system]
    elev = model.B["elev"]
    z0 = min(elev.values())
    H = (max(elev.values()) - z0) / 1e3
    Ta = period_approx(H, material)
    w = story_weights(model, live=live)
    W = sum(w.values())
    Pgrav = story_weights(model, dead=1.0, live=1.0)        # total unfactored vertical load for theta
    out = {"H_m": H, "Ta_s": Ta, "W_kN": W, "R": R, "Cd": Cd, "I": site.I, "SDS": site.SDS, "SD1": site.SD1,
           "Ts": site.Ts, "SDC": site.sdc(Ta), "drift_limit": DRIFT_LIMIT[site.occupancy],
           "theta_max": min(0.5 / (beta * Cd), 0.25)}
    stories = [s for s in model.B["stories"] if s in w]
    for d in ("X", "Y"):
        T = min(T_model[d], 1.5 * Ta) if T_model else Ta
        Sa = site.sa_elf(T)
        Cs = max(Sa * site.I / R, 0.01)
        V = Cs * W
        k = k_exponent(T)
        wh = {s: w[s] * ((elev[s] - z0) / 1e3) ** k for s in stories}
        F = {s: V * v / sum(wh.values()) for s, v in wh.items()}
        rows = []
        for r in story_drift(model, d, forces=F):
            Px = sum(Pgrav.get(s, 0.0) for s in stories if elev[s] >= elev[r["story"]] - 1)
            dxe = abs(r["drift_cm_mm"]) / 1e3
            dx = Cd * dxe / site.I
            dx_edge = Cd * abs(r["drift_max_mm"]) / 1e3 / site.I
            theta = Px * dx / (r["V_kN"] * r["h_m"] * Cd) if r["V_kN"] else float("nan")
            rows.append(dict(story=r["story"], h_m=r["h_m"], F_kN=r["F_kN"], V_kN=r["V_kN"], Px_kN=Px,
                             drift_elastic_mm=dxe * 1e3, drift_design_mm=dx * 1e3,
                             drift_ratio=dx / r["h_m"], drift_ratio_edge=dx_edge / r["h_m"],
                             torsion_ratio=r["torsion_ratio"], theta=theta))
        out[d] = dict(T_used_s=T, Sa_g=Sa, Cs=Cs, V_kN=V, k=k, stories=rows)
    return out


def format_elf(res):
    o = [f"SDS {res['SDS']:.4f} g, SD1 {res['SD1']:.4f} g, Ts {res['Ts']:.3f} s, SDC {res['SDC']}, "
         f"H {res['H_m']:.1f} m, Ta {res['Ta_s']:.2f} s, W {res['W_kN']:.0f} kN, R {res['R']}, Cd {res['Cd']}, I {res['I']}"]
    for d in ("X", "Y"):
        r = res[d]
        o.append(f"\nDirection {d}: T used {r['T_used_s']:.2f} s, Sa {r['Sa_g']:.4f} g, Cs {r['Cs']:.4f}, "
                 f"V {r['V_kN']:.1f} kN, k {r['k']:.2f}")
        o.append(f"{'story':<14}{'F kN':>8}{'V kN':>9}{'dxe mm':>9}{'dx mm':>8}{'dx/h':>9}{'edge':>9}{'theta':>8}{'torsion':>9}")
        for s in r["stories"][::-1]:
            o.append(f"{s['story']:<14}{s['F_kN']:>8.1f}{s['V_kN']:>9.1f}{s['drift_elastic_mm']:>9.3f}{s['drift_design_mm']:>8.2f}"
                     f"{s['drift_ratio']:>9.5f}{s['drift_ratio_edge']:>9.5f}{s['theta']:>8.3f}{s['torsion_ratio']:>9.2f}")
    o.append(f"\nallowable drift ratio {res['drift_limit']:.3f}, theta_max {res['theta_max']:.3f} (theta <= 0.10: P-Delta may be ignored)")
    return "\n".join(o)
