"""Mesh builder: parsed .e2k (N, mm) -> nodes, shells and frame segments.

Fixed conventions (v1.0 of the in-house converter, extended in openstructure 0.1):
  * one global structured plan grid through every key coordinate of assigned objects
    (merge tol 25 mm), subdivided to <= SLAB_H (ETABS FLOORMESHMAXSIZE), <= WALL_H in wall ranges
  * walls meshed vertically <= WALL_H and always at every story elevation they cross
  * area/line connectivity flag k = vertex k stories below the assigned story
  * frames split at every node lying on them
Not modelled (reported): rigid end zones, cardinal-point offsets, panel zones, links,
tendons/PT balancing loads, wind, frame shear deformation.
"""
import math
from collections import defaultdict

from .units import G_ACC

# ---------------------------------------------------------------- geometry helpers
def _pip(x, y, poly):
    inside = False; n = len(poly)
    for i in range(n):
        x1, y1 = poly[i]; x2, y2 = poly[(i + 1) % n]
        if (y1 > y) != (y2 > y):
            xi = x1 + (y - y1) * (x2 - x1) / (y2 - y1)
            if xi > x: inside = not inside
    return inside

def _subdivide(coords, hmax, fine=None, hfine=None):
    """coords: sorted key coordinates (mm). Split intervals longer than hmax
    (or hfine inside the [a,b] ranges listed in `fine`)."""
    out = [coords[0]]
    for a, b in zip(coords[:-1], coords[1:]):
        h = hmax
        if fine and any(lo - 1 <= a and b <= hi + 1 for lo, hi in fine): h = hfine
        n = max(1, math.ceil((b - a) / h - 1e-9))
        out += [a + (b - a) * k / n for k in range(1, n + 1)]
    return out

def _merge(vals, tol):
    vals = sorted(vals); out = [vals[0]]
    for v in vals[1:]:
        if v - out[-1] > tol: out.append(v)
    return out

# ---------------------------------------------------------------- model builder
def build(M, SLAB_H=1250.0, WALL_H=500.0, MERGE_TOL=25.0):
    pts, elev, stories = M["pts"], M["elev"], M["story_list"]
    sidx = {s: i for i, s in enumerate(stories)}
    def below(s, k=1):
        return stories[max(0, sidx[s] - k)]
    asg_line = defaultdict(list)
    for d in M["lassign"]: asg_line[d["name"]].append(d)
    asg_area = defaultdict(list)
    for d in M["aassign"]: asg_area[d["name"]].append(d)
    warn = []

    def pz(pid, story, flag):
        s = below(story, flag) if flag else story
        return elev[s] - pts[pid][2]

    # ---- key plan coordinates from objects that are actually assigned
    kx, ky = set(), set(); wall_ranges_x, wall_ranges_y = [], []
    for nm, (typ, i, j, nst) in M["lines"].items():
        if nm in asg_line:
            for p in (i, j): kx.add(pts[p][0]); ky.add(pts[p][1])
    for nm, (typ, ids, fl) in M["areas"].items():
        if nm in asg_area:
            for p in ids: kx.add(pts[p][0]); ky.add(pts[p][1])
            if typ == "PANEL" and any(d["sec"] in M["sprop"] and "WALLTHICKNESS" in M["sprop"][d["sec"]] for d in asg_area[nm]):
                xs = [pts[p][0] for p in ids]; ys = [pts[p][1] for p in ids]
                wall_ranges_x.append((min(xs), max(xs))); wall_ranges_y.append((min(ys), max(ys)))
    gx = _subdivide(_merge(kx, MERGE_TOL), SLAB_H, [r for r in wall_ranges_x if r[1] > r[0]], WALL_H)
    gy = _subdivide(_merge(ky, MERGE_TOL), SLAB_H, [r for r in wall_ranges_y if r[1] > r[0]], WALL_H)
    snapx = lambda v: min(gx, key=lambda g: abs(g - v)); snapy = lambda v: min(gy, key=lambda g: abs(g - v))

    nodes = {}                       # key (mm ints) -> tag
    coords = []
    def node(x, y, z):
        k = (round(x), round(y), round(z))
        if k not in nodes:
            nodes[k] = len(coords) + 1; coords.append(k)
        return nodes[k]

    shells, frames = [], []          # shells: dict(n=[4], sec, obj, story, kind)
    def lines_between(a, b, grid):
        lo, hi = min(a, b), max(a, b)
        return [g for g in grid if lo - 1 <= g <= hi + 1]

    # ---- slabs (FLOOR objects at story level, no DZ)
    slab_polys = defaultdict(list)   # story -> list of (poly, sec, obj)
    stair_quads = []                 # sloped / offset floors and panels
    wall_panels = []
    for nm, (typ, ids, fl) in M["areas"].items():
        for d in asg_area.get(nm, []):
            s, sec = d["story"], d["sec"]
            zs = [pz(p, s, f) for p, f in zip(ids, fl)]
            P = [(snapx(pts[p][0]), snapy(pts[p][1])) for p in ids]
            if typ == "FLOOR" and max(zs) - min(zs) < 1 and abs(zs[0] - elev[s]) < 1:
                slab_polys[s].append((P, sec, nm, d))
            elif typ == "PANEL" and len(set(P)) == 2:
                wall_panels.append((P, zs, sec, nm, s, d))
            else:
                stair_quads.append((P, zs, sec, nm, s, d))
    for s, polys in slab_polys.items():
        z = elev[s]
        for i in range(len(gx) - 1):
            for j in range(len(gy) - 1):
                xc, yc = (gx[i] + gx[i + 1]) / 2, (gy[j] + gy[j + 1]) / 2
                hit = [(sec, nm, d) for P, sec, nm, d in polys if _pip(xc, yc, P)]
                if hit:
                    # overlapping floor objects: the first one gives section/weight/mass (ETABS does not
                    # double the slab), but area loads assigned to every overlapping object are applied
                    sec, nm, d = hit[0]
                    n = [node(gx[i], gy[j], z), node(gx[i + 1], gy[j], z), node(gx[i + 1], gy[j + 1], z), node(gx[i], gy[j + 1], z)]
                    shells.append(dict(n=n, sec=sec, obj=nm, objs=[h[1] for h in hit], story=s, kind="slab", d=d,
                                       area=(gx[i + 1] - gx[i]) * (gy[j + 1] - gy[j]) / 1e6))
    # ---- walls (vertical panels)
    for P, zs, sec, nm, s, d in wall_panels:
        (x1, y1), (x2, y2) = sorted(set(P))
        zb, zt = min(zs), max(zs)
        if abs(x1 - x2) < 1:   hl = [(x1, y) for y in lines_between(y1, y2, gy)]
        elif abs(y1 - y2) < 1: hl = [(x, y1) for x in lines_between(x1, x2, gx)]
        else: warn.append(f"wall {nm} not axis-aligned - skipped"); continue
        zk = sorted({zb, zt} | {z for z in elev.values() if zb + 1 < z < zt - 1})
        zl = [zk[0]]
        for za_, zb_2 in zip(zk[:-1], zk[1:]):
            nv = max(2 if len(zk) == 2 else 1, math.ceil((zb_2 - za_) / WALL_H - 1e-9))
            zl += [za_ + (zb_2 - za_) * k / nv for k in range(1, nv + 1)]
        for a, b in zip(hl[:-1], hl[1:]):
            L = math.hypot(b[0] - a[0], b[1] - a[1])
            for z0, z1 in zip(zl[:-1], zl[1:]):
                n = [node(*a, z0), node(*b, z0), node(*b, z1), node(*a, z1)]
                shells.append(dict(n=n, sec=sec, obj=nm, story=s, kind="wall", d=d, area=L * (z1 - z0) / 1e6))
    # ---- stairs / offset or sloped areas: plan-rectangle quads with bilinear z
    for P, zs, sec, nm, s, d in stair_quads:
        xs = sorted(set(p[0] for p in P)); ys = sorted(set(p[1] for p in P))
        if len(P) != 4 or len(xs) != 2 or len(ys) != 2:
            if len(xs) == 1 or len(ys) == 1:            # vertical panel with sloped edges -> treat as wall strip
                warn.append(f"area {nm}@{s}: vertical non-rectangular panel - skipped"); continue
            warn.append(f"area {nm}@{s}: not a plan rectangle - skipped"); continue
        zc = {(p[0], p[1]): z for p, z in zip(P, zs)}
        (xa, xb), (ya, yb) = xs, ys
        z00, z10, z11, z01 = zc[(xa, ya)], zc[(xb, ya)], zc[(xb, yb)], zc[(xa, yb)]
        def zb_(x, y):
            u, v = (x - xa) / (xb - xa), (y - ya) / (yb - ya)
            return z00 * (1 - u) * (1 - v) + z10 * u * (1 - v) + z11 * u * v + z01 * (1 - u) * v
        lx, ly = lines_between(xa, xb, gx), lines_between(ya, yb, gy)
        for i in range(len(lx) - 1):
            for j in range(len(ly) - 1):
                c = [(lx[i], ly[j]), (lx[i + 1], ly[j]), (lx[i + 1], ly[j + 1]), (lx[i], ly[j + 1])]
                n = [node(x, y, zb_(x, y)) for x, y in c]
                P3 = [(x, y, zb_(x, y)) for x, y in c]
                ux = [P3[1][k] - P3[0][k] for k in range(3)]; uy = [P3[3][k] - P3[0][k] for k in range(3)]
                cr = [ux[1] * uy[2] - ux[2] * uy[1], ux[2] * uy[0] - ux[0] * uy[2], ux[0] * uy[1] - ux[1] * uy[0]]
                shells.append(dict(n=n, sec=sec, obj=nm, story=s, kind="stair", d=d, area=math.sqrt(sum(v * v for v in cr)) / 1e6))
    # ---- frames
    for nm, (typ, pi, pj, nst) in M["lines"].items():
        for d in asg_line.get(nm, []):
            s = d["story"]
            if typ == "COLUMN":
                x, y = snapx(pts[pi][0]), snapy(pts[pi][1])
                zb, zt = elev[below(s, max(nst, 1))] - pts[pi][2], elev[s] - pts[pj][2]
                frames.append(dict(kind="col", obj=nm, story=s, d=d, ends=((x, y, zb), (x, y, zt))))
            else:
                a = (snapx(pts[pi][0]), snapy(pts[pi][1]), elev[s] - pts[pi][2])
                b = (snapx(pts[pj][0]), snapy(pts[pj][1]), elev[s] - pts[pj][2])
                frames.append(dict(kind="beam", obj=nm, story=s, d=d, ends=(a, b)))
    # split frames at every node lying on them (connect to slab/wall/stair nodes)
    by_xy = defaultdict(list)
    for k in nodes: by_xy[(k[0], k[1])].append(k[2])
    segs = []
    for f in frames:
        (xa, ya, za), (xb, yb, zb) = f["ends"]
        if f["kind"] == "col":
            zz = sorted(set([round(za), round(zb)] + [z for z in by_xy.get((round(xa), round(ya)), []) if round(za) < z < round(zb)]))
            for z0, z1 in zip(zz[:-1], zz[1:]):
                segs.append(dict(f, n=(node(xa, ya, z0), node(xa, ya, z1)), L=(z1 - z0) / 1e3))
        else:
            L = math.dist((xa, ya, za), (xb, yb, zb)); ts = {0.0, 1.0}
            for k in nodes:
                v = (k[0] - xa, k[1] - ya, k[2] - za); dv = (xb - xa, yb - ya, zb - za)
                t = sum(p * q for p, q in zip(v, dv)) / (L * L)
                if 1e-6 < t < 1 - 1e-6:
                    pp = [xa + t * dv[0], ya + t * dv[1], za + t * dv[2]]
                    if math.dist(pp, k) < 2: ts.add(t)
            ts = sorted(ts)
            for t0, t1 in zip(ts[:-1], ts[1:]):
                p0 = [xa + t0 * (xb - xa), ya + t0 * (yb - ya), za + t0 * (zb - za)]
                p1 = [xa + t1 * (xb - xa), ya + t1 * (yb - ya), za + t1 * (zb - za)]
                segs.append(dict(f, n=(node(*p0), node(*p1)), L=(t1 - t0) * L / 1e3))
    # restraints (base points + every node at a restrained level/position)
    fixed = {}
    base_z = elev[stories[0]]
    for (pid, s), dofs in M["rest"].items():
        if s == stories[0]:
            x, y = snapx(pts[pid][0]), snapy(pts[pid][1])
            k = (round(x), round(y), round(elev[s] - pts[pid][2]))
            if k in nodes: fixed[nodes[k]] = dofs
    # wall base mesh nodes (ETABS ADDRESTRAINT on walls / auto-mesh restraint at base)
    wall_base = set()
    for sh in shells:
        if sh["kind"] == "wall":
            for t in sh["n"]:
                if coords[t - 1][2] == round(base_z): wall_base.add(t)
    for t in wall_base: fixed.setdefault(t, ["UX", "UY", "UZ"])
    return dict(coords=coords, shells=shells, frames=segs, fixed=fixed, gx=gx, gy=gy, warn=warn,
                elev=dict(elev), stories=stories)

# ---------------------------------------------------------------- section properties
def rect_J(b, d):
    a, c = max(b, d), min(b, d); r = c / a
    return a * c ** 3 * (1 / 3 - 0.21 * r * (1 - r ** 4 / 12))

def frame_props(M, d):
    s = M["fsec"][d["sec"]]; mat = M["mats"][s["mat"]]
    E = mat["E"] * 1e3; nu = mat.get("U", 0.2); G = E / (2 * (1 + nu)); gam = mat["WEIGHTPERVOLUME"] * 1e6
    if s["shape"] != "Concrete Rectangular":
        raise NotImplementedError(s["shape"])
    t3, t2 = s["D"] / 1e3, s["B"] / 1e3          # depth along local-2, width along local-3
    A = t2 * t3; I33 = t2 * t3 ** 3 / 12; I22 = t3 * t2 ** 3 / 12; J = rect_J(t2, t3)
    fA = s.get("AMOD", 1) * d.get("PROPMODA", 1)
    f22 = s.get("I2MOD", 1) * d.get("PROPMODI22", 1); f33 = s.get("I3MOD", 1) * d.get("PROPMODI33", 1)
    fJ = s.get("JMOD", 1) * d.get("PROPMODT", 1)
    fM = s.get("MMOD", 1) * d.get("PROPMODM", 1); fW = s.get("WMOD", 1) * d.get("PROPMODW", 1)
    return dict(E=E, G=G, A=A * fA, I22=I22 * f22, I33=I33 * f33, J=J * fJ, t2=t2, t3=t3,
                wpl=gam * A * fW, mpl=gam * A * fM / G_ACC)

def shell_props(M, d):
    p = M["sprop"][d["sec"]]; mat = M["mats"][p["mat"]]
    t = (p.get("SLABTHICKNESS") or p.get("WALLTHICKNESS")) / 1e3
    fF = p.get("F11MOD", 1) * d.get("PROPMODF11", 1)       # membrane (F11=F22=F12 assumed)
    fM_ = p.get("M11MOD", 1) * d.get("PROPMODM11", 1)      # bending
    fM12 = p.get("M12MOD", 1) * d.get("PROPMODM12", 1)     # twisting
    mm = p.get("MMOD", 1) * d.get("PROPMODM", 1); wm = p.get("WMOD", 1) * d.get("PROPMODW", 1)
    gam = mat["WEIGHTPERVOLUME"] * 1e6
    return dict(E=mat["E"] * 1e3, nu=mat.get("U", 0.2), t=t, fF=fF, fM=fM_, fM12=fM12, wpa=gam * t * wm, mpa=gam * t * mm / G_ACC,
                mtype=p.get("mtype", "ShellThin"))
