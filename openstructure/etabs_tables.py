"""Reader for ETABS tables exported as text (File > Export > ETABS tables to text/.e2k-style).

Format:  TABLE:  "NAME"  followed by rows of  Key=Value  pairs; a row continues on the
next line when it ends with " _". Keys/values may be quoted.
"""
import re
from collections import defaultdict

from .units import factors

_PAIR = re.compile(r'("[^"]*"|[^\s=]+)=("[^"]*"|\S+)')


def _clean(s):
    return s[1:-1] if s.startswith('"') and s.endswith('"') else s


def _val(s):
    s = _clean(s)
    try:
        return float(s)
    except ValueError:
        return s


def read_tables(path, wanted=None):
    """Return {table name: [row dict, ...]}. `wanted` limits which tables are kept."""
    tables = defaultdict(list)
    name, buf = None, ""
    with open(path, encoding="utf-8", errors="replace") as fh:
        for raw in fh:
            ln = raw.rstrip("\n")
            if ln.startswith("TABLE:"):
                name = _clean(ln.split(":", 1)[1].strip())
                buf = ""
                continue
            if name is None or (wanted and name not in wanted):
                continue
            s = ln.strip()
            if not s:
                continue
            if s.endswith(" _"):
                buf += s[:-2] + " "
                continue
            row = buf + s
            buf = ""
            tables[name].append({_clean(k): _val(v) for k, v in _PAIR.findall(row)})
    return dict(tables)


def etabs_reference(path):
    """Extract what we verify against: base reactions FZ per static case (kN),
    modal periods and participating mass ratios. Units read from PROGRAM CONTROL."""
    T = read_tables(path, wanted={"PROGRAM CONTROL", "BASE REACTIONS", "MODAL PERIODS AND FREQUENCIES",
                                  "MODAL PARTICIPATING MASS RATIOS"})
    units = [u.strip().upper() for u in T["PROGRAM CONTROL"][0]["CurrUnits"].split(",")]
    fF, fL = factors(units)
    to_kN = fF / 1e3
    out = {"units": units, "base_FZ_kN": {}, "periods_s": [], "UX": [], "UY": [], "RZ": []}
    for r in T.get("BASE REACTIONS", []):
        if r.get("Case Type") == "LinStatic" and "Step Number" not in r:
            out["base_FZ_kN"][str(r["Output Case"])] = r["FZ"] * to_kN
    for r in T.get("MODAL PARTICIPATING MASS RATIOS", []):
        out["periods_s"].append(r["Period"])
        for k in ("UX", "UY", "RZ"):
            out[k].append(r[k])
    return out


def story_forces(path, case, bottom_story=None):
    """Story forces of a static case from the ETABS "STORY FORCES" table, in kN and m.
    Returns {step number (1 if the case has no steps): {"X": {story: Fx}, "Y": {story: Fy},
    "y_of_X": {story: y of the Fx line of action}, "x_of_Y": {story: x of the Fy line of action}}},
    obtained by differencing story shear and torsion from the top down to bottom_story (inclusive).
    Use bottom_story for wind cases whose exposure stops above the base."""
    T = read_tables(path, wanted={"PROGRAM CONTROL", "STORY FORCES"})
    units = [u.strip().upper() for u in T["PROGRAM CONTROL"][0]["CurrUnits"].split(",")]
    fF, fL = factors(units)
    k, km = fF / 1e3, fF * fL / 1e6
    shear = defaultdict(dict)
    order = []
    for r in T["STORY FORCES"]:
        if r.get("Output Case") != case or r.get("Location") != "Bottom":
            continue
        st = str(r["Story"])
        if st not in order:
            order.append(st)
        shear[int(r.get("Step Number", 1))][st] = (r["VX"] * k, r["VY"] * k, r["T"] * km)
    out = {}
    for step, sh in shear.items():
        res = {"X": {}, "Y": {}, "y_of_X": {}, "x_of_Y": {}}
        px = py = pt = 0.0
        vmax = max(max(abs(v[0]), abs(v[1])) for v in sh.values())
        for st in order:                       # the table lists stories from the top down
            vx, vy, t = sh[st]
            dx, dy, dt = vx - px, vy - py, t - pt
            px, py, pt = vx, vy, t
            # story shear is minus the applied force; a case acts in one direction per step
            if abs(dx) > 1e-3 * vmax and abs(dx) >= abs(dy):
                res["X"][st] = -dx
                res["y_of_X"][st] = -dt / dx
            elif abs(dy) > 1e-3 * vmax:
                res["Y"][st] = -dy
                res["x_of_Y"][st] = dt / dy
            if st == bottom_story:
                break
        out[step] = res
    return out
