"""Base reactions per support point, for load patterns and combinations.

Supports are reported by ETABS point label. Wall-base mesh nodes that are not ETABS points are added
to the nearest labelled support of the same wall line group (ETABS reports wall reactions only at the
panel corner points). Loads lumped directly on restrained nodes (e.g. half of the bottom column
segment's self-weight) never enter the OpenSees solution, so they are added back to Fz here.
Units kN, kN-m; Fz positive = compression on the foundation (upward reaction); negative = uplift.
"""
import math
from collections import defaultdict


def support_groups(model):
    """{label: [support nodes]} and {label: (x m, y m)}."""
    coords = model.B["coords"]
    labels = model.B["support_labels"]
    lab_xy = {lab: coords[n - 1][:2] for n, lab in labels.items()}
    groups = defaultdict(list)
    for n in model.fixed:
        if n in labels:
            groups[labels[n]].append(n)
        else:
            x, y, _ = coords[n - 1]
            near = min(lab_xy, key=lambda L: math.hypot(lab_xy[L][0] - x, lab_xy[L][1] - y))
            groups[near].append(n)
    return dict(groups), {lab: (xy[0] / 1e3, xy[1] / 1e3) for lab, xy in lab_xy.items()}


def base_reactions(model, combos=None, patterns=None):
    """Returns {case name: {label: dict(x, y, Fx, Fy, Fz, Mx, My, Mz)}} for every load pattern
    (or the listed ones) and every combination {name: {pattern: factor}}."""
    groups, xy = support_groups(model)
    pats = set(patterns or [p for p in model.patterns if model.loads.get(p) or model.lateral.get(p)])
    for c in (combos or {}).values():
        pats |= set(c)
    pats = [p for p in model.patterns if p in pats and (model.loads.get(p) or model.lateral.get(p))]
    R = model.reactions(pats)
    per = {}
    for p in pats:
        per[p] = {}
        for lab, nodes in groups.items():
            v = [sum(R[p][n][i] for n in nodes) for i in range(6)]
            v[2] += sum(model.loads.get(p, {}).get(n, 0.0) for n in nodes)      # loads applied on the support itself
            per[p][lab] = v
    out = {p: per[p] for p in (patterns or pats) if p in per}
    for name, fac in (combos or {}).items():
        out[name] = {lab: [sum(f * per[p][lab][i] for p, f in fac.items() if p in per) for i in range(6)] for lab in groups}
    keys = ("Fx", "Fy", "Fz", "Mx", "My", "Mz")
    return {case: {lab: dict(x=xy[lab][0], y=xy[lab][1], **dict(zip(keys, v))) for lab, v in d.items()} for case, d in out.items()}


def format_reactions(res, cases=None, component="Fz"):
    cases = cases or list(res)
    labs = sorted(next(iter(res.values())), key=lambda s: (len(s), s))
    o = [f"{'point':<7}{'x m':>7}{'y m':>7}" + "".join(f"{c:>12}" for c in cases)]
    for lab in labs:
        r0 = res[cases[0]][lab]
        o.append(f"{lab:<7}{r0['x']:>7.2f}{r0['y']:>7.2f}" + "".join(f"{res[c][lab][component]:>12.1f}" for c in cases))
    o.append(f"{'total':<21}" + "".join(f"{sum(res[c][l][component] for l in labs):>12.1f}" for c in cases))
    return "\n".join(o)
