"""Parser for the ETABS .e2k text export.

All numbers are normalised to N and mm on reading, whatever the file's own
UNITS line says (KGF-M, KN-M, N-MM, TONF-M, LB-IN, ...). Only the user's own
text export is read; nothing of ETABS itself is reverse-engineered.
"""
import re
from collections import defaultdict, OrderedDict

from .units import factors


def _q(s):
    return re.findall(r'"([^"]*)"', s)


def _num(line, key, default=None):
    m = re.search(rf'\b{key}\s+([-\d.Ee+]+)', line)
    return float(m.group(1)) if m else default


def read_sections(path):
    txt = open(path, encoding="utf-8", errors="replace").read()
    secs = OrderedDict()
    sec = None
    for ln in txt.splitlines():
        if ln.startswith("$ "):
            sec = ln[2:].strip()
            secs[sec] = []
        elif sec and ln.strip():
            secs[sec].append(ln.strip())
    return secs


def parse_e2k(path):
    """Return a dict describing the model, in N and mm."""
    secs = read_sections(path)
    M = {"secs": secs, "path": str(path)}
    units = _q(secs["CONTROLS"][0])
    fF, fL = factors(units)
    M["units_in_file"] = units
    stress = fF / fL ** 2          # file F/L^2 -> N/mm^2
    density = fF / fL ** 3         # file F/L^3 -> N/mm^3

    # ---- stories (bottom -> top), elevation in mm
    st = []
    for l in secs["STORIES - IN SEQUENCE FROM TOP"]:
        n = _q(l)[0]
        h = _num(l, "HEIGHT", 0.0) * fL
        e = _num(l, "ELEV", None)
        st.append((n, h, None if e is None else e * fL))
    st = st[::-1]
    elev = OrderedDict()
    z = 0.0
    for i, (n, h, e) in enumerate(st):
        z = (e if e is not None else 0.0) if i == 0 else z + h
        elev[n] = z
    M["elev"] = elev
    M["story_list"] = list(elev)

    # ---- points: X, Y, DZ (below the story the object is assigned to)
    pts = {}
    for l in secs["POINT COORDINATES"]:
        t = l.split()
        pid = _q(l)[0]
        dz = float(t[4]) if len(t) > 4 else 0.0
        pts[pid] = (float(t[2]) * fL, float(t[3]) * fL, dz * fL)
    M["pts"] = pts

    # ---- materials: E, FC, FY in N/mm^2, WEIGHTPERVOLUME in N/mm^3
    mats = defaultdict(dict)
    for l in secs["MATERIAL PROPERTIES"]:
        n = _q(l)[0]
        for k, f in (("E", stress), ("U", 1.0), ("FC", stress), ("FY", stress), ("WEIGHTPERVOLUME", density)):
            v = _num(l, k)
            if v is not None and k not in mats[n]:
                mats[n][k] = v * f
        m = re.search(r'TYPE "([^"]+)"', l)
        if m:
            mats[n]["TYPE"] = m.group(1)
    M["mats"] = dict(mats)

    # ---- frame sections (dimensions in mm)
    fs = defaultdict(dict)
    for l in secs["FRAME SECTIONS"]:
        n = _q(l)[0]
        d = fs[n]
        m = re.search(r'MATERIAL "([^"]+)"', l)
        if m:
            d["mat"] = m.group(1)
        m = re.search(r'SHAPE "([^"]+)"', l)
        if m:
            d["shape"] = m.group(1)
        for k in ("D", "B", "TF", "TW"):
            v = _num(l, k)
            if v is not None:
                d[k] = v * fL
        for k in ("AMOD", "A2MOD", "A3MOD", "JMOD", "I2MOD", "I3MOD", "MMOD", "WMOD"):
            v = _num(l, k)
            if v is not None:
                d[k] = v
    M["fsec"] = dict(fs)

    # ---- shell properties (thickness in mm)
    sp = defaultdict(dict)
    for sname in ("SLAB PROPERTIES", "WALL PROPERTIES", "DECK PROPERTIES"):
        for l in secs.get(sname, []):
            n = _q(l)[0]
            d = sp[n]
            m = re.search(r'MATERIAL "([^"]+)"', l)
            if m and "mat" not in d:
                d["mat"] = m.group(1)
            for k in ("SLABTHICKNESS", "WALLTHICKNESS"):
                v = _num(l, k)
                if v is not None:
                    d[k] = v * fL
            for k in ("F11MOD", "F22MOD", "F12MOD", "M11MOD", "M22MOD", "M12MOD", "MMOD", "WMOD"):
                v = _num(l, k)
                if v is not None:
                    d[k] = v
            m = re.search(r'MODELINGTYPE "([^"]+)"', l)
            if m:
                d["mtype"] = m.group(1)
    M["sprop"] = dict(sp)

    # ---- plan objects
    lines = {}
    for l in secs["LINE CONNECTIVITIES"]:
        t = _q(l)
        w = l.split()
        lines[t[0]] = (w[2], t[1], t[2], int(w[-1]))
    M["lines"] = lines
    areas = {}
    for l in secs["AREA CONNECTIVITIES"]:
        t = _q(l)
        w = l.split()
        n = int(w[3])
        ids = t[1:1 + n]
        flags = [int(x) for x in w[4 + n:4 + 2 * n]]
        areas[t[0]] = (w[2], ids, flags)
    M["areas"] = areas

    la = []
    for l in secs["LINE ASSIGNS"]:
        t = _q(l)
        d = {"name": t[0], "story": t[1], "sec": re.search(r'SECTION "([^"]+)"', l).group(1),
             "ang": _num(l, "ANG", 0.0), "cardinal": int(_num(l, "CARDINALPT", 10))}
        for k in ("PROPMODA", "PROPMODA2", "PROPMODA3", "PROPMODT", "PROPMODI22", "PROPMODI33", "PROPMODM", "PROPMODW"):
            v = _num(l, k)
            if v is not None:
                d[k] = v
        la.append(d)
    M["lassign"] = la
    aa = []
    for l in secs["AREA ASSIGNS"]:
        t = _q(l)
        m = re.search(r'SECTION "([^"]+)"', l)
        if not m:
            continue
        d = {"name": t[0], "story": t[1], "sec": m.group(1)}
        for k in ("PROPMODF11", "PROPMODF22", "PROPMODF12", "PROPMODM11", "PROPMODM22", "PROPMODM12", "PROPMODM", "PROPMODW"):
            v = _num(l, k)
            if v is not None:
                d[k] = v
        aa.append(d)
    M["aassign"] = aa

    rest = {}
    for l in secs["POINT ASSIGNS"]:
        if "RESTRAINT" in l:
            t = _q(l)
            rest[(t[0], t[1])] = re.search(r'RESTRAINT "([^"]+)"', l).group(1).split()
    M["rest"] = rest

    # ---- load patterns: name -> (type, self-weight multiplier)
    pats = OrderedDict()
    for l in secs.get("LOAD PATTERNS", []):
        if l.startswith("LOADPATTERN"):
            t = _q(l)
            m = re.search(r'TYPE\s+"([^"]+)"', l)
            pats[t[0]] = (m.group(1) if m else "Other", _num(l, "SELFWEIGHT", 0.0))
    if not pats:
        pats["Dead"] = ("Dead", 1.0)
    M["patterns"] = pats

    # ---- area loads (N/mm^2)
    lsets = defaultdict(dict)
    for l in secs.get("SHELL UNIFORM LOAD SETS", []):
        t = _q(l)
        lsets[t[0]][t[1]] = _num(l, "VALUE") * stress
    M["lsets"] = dict(lsets)
    aload = []
    for l in secs.get("SHELL OBJECT LOADS", []):
        t = _q(l)
        if "UNIFLOADSET" in l:
            aload.append((t[0], t[1], "SET", t[3]))
        else:
            m = re.search(r'LOADPAT "([^"]+)"\s+VALUE\s+([-\d.Ee+]+)', l)
            if m:
                aload.append((t[0], t[1], m.group(1), float(m.group(2)) * stress))
    M["aload"] = aload

    # ---- features read but not modelled yet (reported, never silently dropped)
    unsupported = []
    if secs.get("POINT OBJECT LOADS"):
        unsupported.append(f"{len(secs['POINT OBJECT LOADS'])} point loads")
    if secs.get("FRAME OBJECT LOADS"):
        unsupported.append(f"{len(secs['FRAME OBJECT LOADS'])} frame loads")
    if secs.get("TENDON SECTIONS"):
        unsupported.append("tendons / PT balancing loads")
    M["unsupported"] = unsupported
    return M
