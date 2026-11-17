"""High-level API: ETABS .e2k -> OpenSeesPy 3D linear model -> analyses.

    from openstructure import Model
    m = Model.from_e2k("building.e2k")
    m.gravity()        # base reactions per load pattern
    m.modal(12)        # periods + modal participating mass ratios
    m.lateral_share()  # wall vs column base shear, story drift, torsion ratio

Units inside the OpenSees model: kN, m, s, t.
OpenSeesPy has one global model per Python process, so only one Model is "live" at a
time; calling an analysis on a Model rebuilds it in OpenSees if another one was built since.
"""
import math
import os
from collections import defaultdict

os.environ.setdefault("MKL_THREADING_LAYER", "SEQUENTIAL")  # MKL crash on some Windows conda builds

import numpy as np
import openseespy.opensees as ops

from . import e2k as _e2k
from . import mesh as _mesh

_LIVE = None   # id of the Model currently held by OpenSees
DOF = {"UX": 0, "UY": 1, "UZ": 2, "RX": 3, "RY": 4, "RZ": 5}


class Model:
    N_LAYERS = 20
    twist_modifier = True   # honour ETABS M12MOD different from M11MOD on slabs

    def __init__(self, parsed, built, name="model"):
        self.M, self.B, self.name = parsed, built, name
        self.coords = np.array(built["coords"], float) / 1e3          # m
        self.story_z = sorted(round(z) for z in built["elev"].values())  # mm
        self.fixed = set(built["fixed"])
        self.tcl = []
        self._frame_forces, self._reactions = {}, {}
        self._prepare_loads_and_mass()

    # ------------------------------------------------------------------ construction
    @classmethod
    def from_e2k(cls, path, slab_mesh=1.25, wall_mesh=0.5):
        """slab_mesh / wall_mesh: max mesh size in m (ETABS FLOORMESHMAXSIZE default 1.25 m)."""
        M = _e2k.parse_e2k(path)
        B = _mesh.build(M, SLAB_H=slab_mesh * 1e3, WALL_H=wall_mesh * 1e3)
        return cls(M, B, name=os.path.splitext(os.path.basename(str(path)))[0])

    @property
    def patterns(self):
        return list(self.M["patterns"])

    def _prepare_loads_and_mass(self):
        M, B = self.M, self.B
        sw_pats = {p: f for p, (typ, f) in M["patterns"].items() if f}
        loads = {p: defaultdict(float) for p in M["patterns"]}
        mass_raw = defaultdict(float)
        wsum = defaultdict(float)
        for sh in B["shells"]:
            sp = _mesh.shell_props(M, sh["d"])
            w, m = sp["wpa"] * sh["area"], sp["mpa"] * sh["area"]
            for p, f in sw_pats.items():
                for n in sh["n"]:
                    loads[p][n] += f * w / 4
            for n in sh["n"]:
                mass_raw[n] += m / 4
            wsum["self-weight " + sh["kind"]] += w
        for fr in B["frames"]:
            pr = _mesh.frame_props(M, fr["d"])
            w, m = pr["wpl"] * fr["L"], pr["mpl"] * fr["L"]
            for p, f in sw_pats.items():
                for n in fr["n"]:
                    loads[p][n] += f * w / 2
            for n in fr["n"]:
                mass_raw[n] += m / 2
            wsum["self-weight " + fr["kind"]] += w
        # uniform area loads: N/mm^2 * 1000 = kPa
        shells_by_obj = defaultdict(list)
        for sh in B["shells"]:
            for o in sh.get("objs", [sh["obj"]]):
                shells_by_obj[(o, sh["story"])].append(sh)
        for obj, story, typ, val in M["aload"]:
            items = M["lsets"][val].items() if typ == "SET" else [(typ, val)]
            for pat, q in items:
                q *= 1e3
                for sh in shells_by_obj[(obj, story)]:
                    for n in sh["n"]:
                        loads.setdefault(pat, defaultdict(float))[n] += q * sh["area"] / 4
                    wsum[f"{pat} area load"] += q * sh["area"]
        # mass source = element self-mass, lateral only, lumped to story levels on the same vertical
        xy_levels = defaultdict(dict)
        for i, k in enumerate(B["coords"], 1):
            xy_levels[(k[0], k[1])][k[2]] = i
        mass = defaultdict(float)
        for n, m in mass_raw.items():
            x, y, z = B["coords"][n - 1]
            if z in self.story_z:
                mass[n] += m
                continue
            lv = xy_levels[(x, y)]
            lo = max([s for s in self.story_z if s < z and s in lv], default=None)
            hi = min([s for s in self.story_z if s > z and s in lv], default=None)
            if lo is None or hi is None:
                mass[n] += m
                continue
            a = (z - lo) / (hi - lo)
            mass[lv[lo]] += m * (1 - a)
            mass[lv[hi]] += m * a
        self.loads, self.mass, self.weight_breakdown = loads, mass, dict(wsum)

    def _cmd(self, name, *a):
        self.tcl.append(name + " " + " ".join(str(x) for x in a))
        getattr(ops, name)(*a)

    def build_ops(self, pdelta=False):
        """Create the model in OpenSees (nodes, fixities, elements, masses).
        pdelta=True uses the PDelta geometric transformation for columns."""
        global _LIVE
        M, B, cmd = self.M, self.B, self._cmd
        self.tcl = []
        ops.wipe()
        cmd("model", "basic", "-ndm", 3, "-ndf", 6)
        for i, (x, y, z) in enumerate(self.coords, 1):
            cmd("node", i, float(x), float(y), float(z))
        for n, dofs in B["fixed"].items():
            f = [0] * 6
            for d in dofs:
                f[DOF[d]] = 1
            cmd("fix", n, *f)
        sec_tags, etag = {}, 0
        for sh in B["shells"]:
            sp = _mesh.shell_props(M, sh["d"])
            key = (sh["sec"], sp["fF"], sp["fM"], sp["fM12"])
            if key not in sec_tags:
                sec_tags[key] = tag = len(sec_tags) + 1
                if self.twist_modifier and abs(sp["fM12"] - sp["fM"]) > 1e-9 and sh["kind"] == "slab":
                    # ETABS M12 (twisting) modifier differs from M11/M22: orthotropic layered plate.
                    # Bending: D11 = D22 = fM*D, twist = fM12*D(1-nu)/2. Membrane uses E*fM too
                    # (slab membrane is diaphragm action; reported as an approximation).
                    E, nu, t = sp["E"] * sp["fM"], sp["nu"], sp["t"]
                    G = sp["E"] / (2 * (1 + nu))
                    cmd("nDMaterial", "ElasticOrthotropic", 1000 + tag, E, E, E, nu, nu, nu, G * sp["fM12"], G, G)
                    cmd("nDMaterial", "PlateFiber", 2000 + tag, 1000 + tag)
                    nl = self.N_LAYERS
                    lay = []
                    for _ in range(nl):
                        lay += [2000 + tag, t / nl]
                    cmd("section", "LayeredShell", tag, nl, *lay)
                else:
                    cmd("section", "ElasticMembranePlateSection", tag, sp["E"] * sp["fF"], sp["nu"], sp["t"], 0.0, sp["fM"] / sp["fF"])
            etag += 1
            cmd("element", "ShellMITC4", etag, *sh["n"], sec_tags[key])
            sh["tag"] = etag
        transf = {}
        for fr in B["frames"]:
            pr = _mesh.frame_props(M, fr["d"])
            i, j = fr["n"]
            if fr["kind"] == "col":
                a = math.radians(fr["d"]["ang"])
                v = (round(math.cos(a), 12), round(math.sin(a), 12), 0.0)   # ETABS local-2
            else:
                if abs(fr["d"]["ang"]) > 1e-6:
                    raise NotImplementedError(f"rotated beam {fr['obj']}")
                v = (0.0, 0.0, 1.0)
            v = v + (fr["kind"] == "col",)
            if v not in transf:
                transf[v] = len(transf) + 1
                cmd("geomTransf", "PDelta" if pdelta and fr["kind"] == "col" else "Linear", transf[v], *v[:3])
            etag += 1
            # vecxz = ETABS local-2 -> Iz = I22, Iy = I33
            cmd("element", "elasticBeamColumn", etag, i, j, pr["A"], pr["E"], pr["G"], pr["J"], pr["I33"], pr["I22"], transf[v])
            fr["tag"] = etag
        for n, m in self.mass.items():
            if n not in self.fixed and m > 0:
                cmd("mass", n, m, m, 0.0, 0.0, 0.0, 0.0)
        _LIVE = (id(self), bool(pdelta))

    def _ensure(self, pdelta=False):
        if _LIVE != (id(self), bool(pdelta)):
            self.build_ops(pdelta)

    # ------------------------------------------------------------------ results helpers
    @property
    def total_mass(self):
        """Lateral mass on free nodes, t."""
        return sum(m for n, m in self.mass.items() if n not in self.fixed)

    def story_of(self, z_mm):
        return min(self.B["elev"], key=lambda k: abs(self.B["elev"][k] - z_mm))

    def story_mass(self):
        out = defaultdict(float)
        for n, m in self.mass.items():
            if n not in self.fixed:
                out[self.story_of(self.B["coords"][n - 1][2])] += m
        return {s: out[s] for s in self.B["stories"] if out[s] > 0}

    def info(self):
        return {"nodes": len(self.coords), "shells": len(self.B["shells"]), "frame_segments": len(self.B["frames"]),
                "supports": len(self.fixed), "units_in_file": self.M["units_in_file"],
                "patterns": self.M["patterns"], "warnings": self.B["warn"], "not_modelled": self.M["unsupported"]}

    @staticmethod
    def _static_setup():
        ops.constraints("Plain")
        ops.numberer("RCM")
        ops.system("UmfPack")
        ops.test("NormDispIncr", 1e-10, 10)
        ops.algorithm("Linear")
        ops.integrator("LoadControl", 1.0)
        ops.analysis("Static")

    # ------------------------------------------------------------------ analyses
    def gravity(self, patterns=None):
        """Linear static per load pattern. Returns {pattern: {applied_kN, base_Rz_kN, min_uz_mm}}."""
        self._ensure()
        out = {}
        for k, pat in enumerate(patterns or self.patterns, 1):
            nl = self.loads.get(pat)
            if not nl:
                continue
            ops.timeSeries("Linear", k)
            ops.pattern("Plain", k, k)
            for n, w in nl.items():
                if n not in self.fixed:
                    ops.load(n, 0, 0, -w, 0, 0, 0)
            self._static_setup()
            ok = ops.analyze(1)
            ops.reactions()
            self._frame_forces[pat] = self._read_frame_forces()
            self._reactions[pat] = {n: ops.nodeReaction(n) for n in self.fixed}
            out[pat] = dict(ok=ok, applied_kN=sum(nl.values()),
                            applied_free_kN=sum(w for n, w in nl.items() if n not in self.fixed),
                            base_Rz_kN=sum(ops.nodeReaction(n, 3) for n in self.fixed),
                            min_uz_mm=1e3 * min(ops.nodeDisp(n, 3) for n in range(1, len(self.coords) + 1)))
            ops.remove("loadPattern", k)
            ops.remove("timeSeries", k)
            ops.wipeAnalysis()
            ops.setTime(0.0)
            ops.reset()
        return out

    def _read_frame_forces(self):
        """End forces of every frame object (column per story, beam) in ETABS-like terms.
        P compression positive (kN); M2, M3 = bending moment about ETABS local 2 and 3 (kN-m) at the
        bottom (i) and top (j) end, same sign convention at both ends (same sign = single curvature);
        V2, V3 shear magnitudes."""
        groups = defaultdict(list)
        for fr in self.B["frames"]:
            groups[(fr["obj"], fr["story"], fr["kind"])].append(fr)
        coords = self.B["coords"]
        out = {}
        for (obj, story, kind), segs in groups.items():
            key = (lambda f: coords[f["n"][0] - 1][2]) if kind == "col" else (lambda f: coords[f["n"][0] - 1][:2])
            segs = sorted(segs, key=key)
            fi = ops.eleResponse(segs[0]["tag"], "localForce")
            fj = ops.eleResponse(segs[-1]["tag"], "localForce")
            # OpenSees local z = ETABS local 2, local y = ETABS local 3 (see build_ops)
            out[(obj, story)] = dict(kind=kind, sec=segs[0]["d"]["sec"], L=sum(f["L"] for f in segs),
                                     P=fi[0], P_top=-fj[6], V2=abs(fi[2]), V3=abs(fi[1]), T=fi[3],
                                     M3_i=fi[4], M3_j=-fj[10], M2_i=fi[5], M2_j=-fj[11])
        return out

    def frame_forces(self, patterns=None):
        """{pattern: {(object, story): forces}} from linear static analysis of each load pattern."""
        pats = [p for p in (patterns or self.patterns) if self.loads.get(p)]
        if any(p not in self._frame_forces for p in pats):
            self.gravity(pats)
        return {p: self._frame_forces[p] for p in pats}

    def reactions(self, patterns=None):
        """{pattern: {support node: [Fx, Fy, Fz, Mx, My, Mz]}} (kN, kN-m)."""
        pats = [p for p in (patterns or self.patterns) if self.loads.get(p)]
        if any(p not in self._reactions for p in pats):
            self.gravity(pats)
        return {p: self._reactions[p] for p in pats}

    def modal(self, n_modes=12):
        """Eigen analysis. Returns periods (s) and modal participating mass ratios (UX, UY, RZ)."""
        self._ensure()
        ops.wipeAnalysis()
        ops.constraints("Plain")
        ops.numberer("RCM")
        ops.system("UmfPack")
        lam = ops.eigen(n_modes)
        T = [2 * math.pi / math.sqrt(l) for l in lam]
        mp = ops.modalProperties("-return")
        self.modes = np.array([[ops.nodeEigenvector(n, k + 1)[:3] for n in range(1, len(self.coords) + 1)]
                               for k in range(n_modes)])
        return {"periods_s": T,
                "UX": list(mp["partiMassRatiosMX"]), "UY": list(mp["partiMassRatiosMY"]),
                "RZ": list(mp["partiMassRatiosRMZ"]), "total_mass_t": mp["totalMass"][:3]}

    def lateral_share(self, total_kN=1000.0):
        """Static pattern F ~ m*z in X and Y: base shear wall vs columns, story column shear,
        story drift and torsion ratio 2*dmax/(dmax+dmin)."""
        self._ensure()
        coords = self.B["coords"]
        z0 = self.story_z[0]
        base_wall = {t for sh in self.B["shells"] if sh["kind"] == "wall" for t in sh["n"] if coords[t - 1][2] == z0}
        pos = {tuple(k): i for i, k in enumerate(coords, 1)}
        wts = {n: m * coords[n - 1][2] / 1e3 for n, m in self.mass.items() if n not in self.fixed and m > 0}
        S = sum(wts.values())
        res = {}
        for kk, (lab, dof) in enumerate((("X", 1), ("Y", 2)), start=20):
            ops.timeSeries("Linear", kk)
            ops.pattern("Plain", kk, kk)
            for n, v in wts.items():
                f6 = [0.0] * 6
                f6[dof - 1] = total_kN * v / S
                ops.load(n, *f6)
            self._static_setup()
            ok = ops.analyze(1)
            ops.reactions()
            Rw = -sum(ops.nodeReaction(n, dof) for n in self.fixed if n in base_wall)
            Rc = -sum(ops.nodeReaction(n, dof) for n in self.fixed if n not in base_wall)
            rows = []
            for si in range(1, len(self.story_z)):
                zb, zt = self.story_z[si - 1], self.story_z[si]
                Vtot = sum(total_kN * v / S for n, v in wts.items() if coords[n - 1][2] > zb)
                Vcol = sum(-ops.eleResponse(f["tag"], "globalForce")[dof - 1] for f in self.B["frames"]
                           if f["kind"] == "col" and coords[f["n"][0] - 1][2] == zb)
                d = [ops.nodeDisp(nh, dof) - ops.nodeDisp(pos[(x, y, zb)], dof)
                     for (x, y, z), nh in pos.items() if z == zt and (x, y, zb) in pos]
                if not d:
                    continue
                d = np.abs(np.array(d))
                dmax, dmin = d.max(), d.min()
                rows.append(dict(story=self.story_of(zt), V_total_kN=Vtot, V_columns_kN=Vcol,
                                 col_share_pct=100 * Vcol / Vtot if Vtot else float("nan"),
                                 drift_max_mm=dmax * 1e3, drift_min_mm=dmin * 1e3,
                                 torsion_ratio=2 * dmax / (dmax + dmin) if dmax + dmin > 0 else float("nan")))
            res[lab] = dict(ok=ok, base_walls_kN=Rw, base_columns_kN=Rc,
                            wall_share_pct=100 * Rw / (Rw + Rc), stories=rows)
            ops.remove("loadPattern", kk)
            ops.remove("timeSeries", kk)
            ops.wipeAnalysis()
            ops.setTime(0.0)
            ops.reset()
        return res

    def story_drift(self, direction, **kw):
        """Story drift, torsion ratio and optional P-Delta. See openstructure.stability.story_drift."""
        from .stability import story_drift
        return story_drift(self, direction, **kw)

    def stability(self, **kw):
        """Stability index Q / theta per story, with a P-Delta check. See openstructure.stability.stability."""
        from .stability import stability
        return stability(self, **kw)

    def write_tcl(self, path):
        """Write every OpenSees command used to build the model (for audit / OpenSees Tcl)."""
        open(path, "w").write("\n".join(self.tcl))
