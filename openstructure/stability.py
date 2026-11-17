"""Story drift, torsion ratio, stability index and P-Delta.

Definitions used here
  * story shear V_x      = sum of the applied lateral forces at and above the story
  * story drift (CM)     = mass-weighted mean drift of all vertical lines with nodes at both
                           levels (equals the centre-of-mass drift for a rigid floor, and stays
                           meaningful where two levels cover different parts of the plan)
  * story drift (max)    = largest drift of any vertical line that has nodes at both levels
  * torsion ratio        = d_max / d_avg with d_avg = (d_max + d_min) / 2
                           (DPT 1301/1302-61 section 2.4.2.1: > 1.2 type 1a, > 1.4 type 1b)
  * stability index      Q = sum(P) * drift_CM / (V * h)
        ACI 318-14 6.6.4.3: with factored gravity, Q > 0.05 means a sway story.
        DPT 1301/1302-61 eq 3.8-1: theta = Px * Delta / (Vx * hsx * Cd) with Delta = Cd * d_xe / I,
        i.e. the same expression with unfactored gravity and I = 1; theta <= 0.10 lets P-Delta be
        ignored and theta_max = 0.5 / (beta * Cd) <= 0.25.
  * P-Delta              = columns with the OpenSees PDelta transformation under constant gravity.
        Walls and slabs are shells without geometric stiffness, so only gravity carried by
        columns produces P-Delta here. The report gives the P-Delta amplification next to 1/(1-Q).
"""
from collections import defaultdict

import numpy as np
import openseespy.opensees as ops

DOF = {"X": 1, "Y": 2}


def default_gravity_factors(model, dead=1.2, live=1.0):
    """Factors by ETABS load pattern type: Dead and Super Dead -> dead, Live types -> live."""
    out = {}
    for name, (typ, _) in model.M["patterns"].items():
        t = typ.lower()
        if "dead" in t:
            out[name] = dead
        elif "live" in t:
            out[name] = live
    return out


def _levels(model):
    """[(story name, z mm)] bottom -> top, base first."""
    return sorted(((s, round(z)) for s, z in model.B["elev"].items()), key=lambda t: t[1])


def story_forces(model, total_kN=1000.0, k=1.0):
    """Lateral force per story proportional to story mass * height^k (DPT eq 3.4-2 shape)."""
    lv = dict(_levels(model))
    z0 = min(lv.values())
    sm = model.story_mass()
    w = {s: m * ((lv[s] - z0) / 1e3) ** k for s, m in sm.items() if lv[s] > z0}
    tot = sum(w.values())
    return {s: total_kN * v / tot for s, v in w.items()}


def gravity_above(model, factors):
    """Sum of factored gravity load applied above the bottom of each story, kN."""
    coords = model.B["coords"]
    lv = _levels(model)
    out = {}
    for (sb, zb), (st, zt) in zip(lv[:-1], lv[1:]):
        P = 0.0
        for pat, f in factors.items():
            for n, w in model.loads.get(pat, {}).items():
                if n not in model.fixed and coords[n - 1][2] > zb:
                    P += f * w
        out[st] = P
    return out


def _solve(nonlinear):
    ops.constraints("Plain")
    ops.numberer("RCM")
    ops.system("UmfPack")
    ops.test("NormDispIncr", 1e-8, 25)
    ops.algorithm("Newton" if nonlinear else "Linear")
    ops.integrator("LoadControl", 1.0)
    ops.analysis("Static")
    ok = ops.analyze(1)
    if ok != 0:
        raise RuntimeError("static analysis did not converge")


def story_drift(model, direction, forces=None, total_kN=1000.0, k=1.0, pdelta=False, gravity_factors=None):
    """Lateral static analysis in 'X' or 'Y'.

    forces: {story name: kN}; if None, total_kN distributed ~ mass * height^k.
    pdelta: True applies gravity (gravity_factors, default 1.2D + 1.0L) first, holds it constant
            and uses the PDelta transformation for columns.
    Returns a list of story rows, bottom -> top.
    """
    dof = DOF[direction.upper()]
    coords = model.B["coords"]
    lv = _levels(model)
    forces = forces or story_forces(model, total_kN, k)
    model._ensure(pdelta)
    ops.wipeAnalysis()
    ops.reset()
    ops.setTime(0.0)
    nn = len(coords)

    # nodes with mass at each level, used both to spread the story force and to average displacement
    by_level = defaultdict(list)
    for n, m in model.mass.items():
        if n not in model.fixed and m > 0:
            by_level[coords[n - 1][2]].append((n, m))

    u0 = np.zeros(nn + 1)
    if pdelta:
        gf = gravity_factors or default_gravity_factors(model)
        ops.timeSeries("Linear", 901)
        ops.pattern("Plain", 901, 901)
        tot = defaultdict(float)
        for pat, f in gf.items():
            for n, w in model.loads.get(pat, {}).items():
                if n not in model.fixed:
                    tot[n] += f * w
        for n, w in tot.items():
            ops.load(n, 0, 0, -w, 0, 0, 0)
        _solve(True)
        ops.loadConst("-time", 0.0)
        u0 = np.array([0.0] + [ops.nodeDisp(n, dof) for n in range(1, nn + 1)])

    ops.timeSeries("Linear", 902)
    ops.pattern("Plain", 902, 902)
    zname = {z: s for s, z in lv}
    for z, lst in by_level.items():
        F = forces.get(zname.get(z), 0.0)
        if not F:
            continue
        msum = sum(m for _, m in lst)
        for n, m in lst:
            f6 = [0.0] * 6
            f6[dof - 1] = F * m / msum
            ops.load(n, *f6)
    _solve(pdelta)
    u = np.array([0.0] + [ops.nodeDisp(n, dof) for n in range(1, nn + 1)]) - u0

    pos = {tuple(c): i for i, c in enumerate(coords, 1)}
    cm = {}
    for s, z in lv:
        lst = by_level.get(z, [])
        cm[s] = sum(u[n] * m for n, m in lst) / sum(m for _, m in lst) if lst else 0.0
    rows = []
    for (sb, zb), (st, zt) in zip(lv[:-1], lv[1:]):
        h = (zt - zb) / 1e3
        V = sum(F for s, F in forces.items() if dict(lv)[s] > zb)
        lines = [(nh, pos[(x, y, zb)]) for (x, y, z), nh in pos.items() if z == zt and (x, y, zb) in pos]
        d = np.array([u[a] - u[b] for a, b in lines])
        w = np.array([model.mass.get(a, 0.0) if a not in model.fixed else 0.0 for a, b in lines])
        if w.sum() <= 0:
            w = np.ones(len(lines))
        dcm = float((d * w).sum() / w.sum())          # mass-weighted mean drift of the story
        sign = 1.0 if dcm >= 0 else -1.0
        dmax, dmin = (sign * d).max(), (sign * d).min()
        rows.append(dict(story=st, h_m=h, F_kN=forces.get(st, 0.0), V_kN=V, disp_cm_mm=cm[st] * 1e3,
                         drift_cm_mm=dcm * 1e3, drift_max_mm=dmax * 1e3, drift_min_mm=dmin * 1e3,
                         drift_ratio_cm=abs(dcm) / h, drift_ratio_max=abs(dmax) / h,
                         torsion_ratio=2 * dmax / (dmax + dmin) if dmax + dmin > 0 else float("nan")))
    ops.remove("loadPattern", 902)
    ops.remove("timeSeries", 902)
    if pdelta:
        ops.remove("loadPattern", 901)
        ops.remove("timeSeries", 901)
    ops.wipeAnalysis()
    ops.reset()
    ops.setTime(0.0)
    return rows


def stability(model, gravity_factors=None, forces=None, total_kN=1000.0, k=1.0, pdelta=True):
    """Stability index per story in X and Y, plus the P-Delta drift amplification.

    forces may be {"X": {story: kN}, "Y": {...}} (e.g. wind or ELF story forces). Q is a ratio,
    so the magnitude of the lateral load does not matter, only its distribution over the height.
    """
    gf = gravity_factors or default_gravity_factors(model)
    P = gravity_above(model, gf)
    out = {"gravity_factors": gf}
    for d in ("X", "Y"):
        f = forces.get(d) if forces else None
        lin = story_drift(model, d, forces=f, total_kN=total_kN, k=k)
        pd = story_drift(model, d, forces=f, total_kN=total_kN, k=k, pdelta=True, gravity_factors=gf) if pdelta else None
        rows = []
        for i, r in enumerate(lin):
            Q = P[r["story"]] * abs(r["drift_cm_mm"]) / 1e3 / (r["V_kN"] * r["h_m"]) if r["V_kN"] else float("nan")
            row = dict(story=r["story"], h_m=r["h_m"], P_kN=P[r["story"]], V_kN=r["V_kN"],
                       drift_cm_mm=r["drift_cm_mm"], drift_max_mm=r["drift_max_mm"],
                       torsion_ratio=r["torsion_ratio"], Q=Q, sway=Q > 0.05, amp_1_over_1mQ=1 / (1 - Q))
            if pd:
                row["drift_cm_pdelta_mm"] = pd[i]["drift_cm_mm"]
                row["amp_pdelta"] = pd[i]["drift_cm_mm"] / r["drift_cm_mm"] if r["drift_cm_mm"] else float("nan")
            rows.append(row)
        out[d] = rows
    return out


def format_stability(res):
    out = ["gravity: " + " + ".join(f"{f:g} {p}" for p, f in res["gravity_factors"].items())]
    for d in ("X", "Y"):
        out.append(f"\nDirection {d}")
        out.append(f"{'story':<14}{'h m':>6}{'P kN':>10}{'V kN':>9}{'drift mm':>10}{'Q':>8}{'1/(1-Q)':>9}{'P-Delta':>9}{'torsion':>9}  class")
        for r in res[d][::-1]:
            out.append(f"{r['story']:<14}{r['h_m']:>6.2f}{r['P_kN']:>10.0f}{r['V_kN']:>9.1f}{r['drift_cm_mm']:>10.3f}"
                       f"{r['Q']:>8.3f}{r['amp_1_over_1mQ']:>9.3f}{r.get('amp_pdelta', float('nan')):>9.3f}"
                       f"{r['torsion_ratio']:>9.2f}  {'sway' if r['sway'] else 'nonsway'}")
    return "\n".join(out)
