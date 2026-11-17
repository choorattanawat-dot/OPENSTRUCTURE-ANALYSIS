"""RC column check per ACI 318-14: P-M2-M3 strength with nonsway moment magnification.

For every column object and story, and every load combination:
  1. Pu, M2, M3 at both ends by linear superposition of the load pattern results
  2. slenderness per axis (6.2.5, 6.6.4.5):  k = 1,  r = 0.3 h,  lu = story height - lu_deduct
       slender if k lu / r > 22 (sway story) or > 34 + 12 (M1/M2) <= 40 (nonsway; M1/M2 negative
       in single curvature)
       (EI)eff = 0.4 Ec Ig / (1 + beta_dns)          (6.6.4.4.4a)
       Pc = pi^2 (EI)eff / (k lu)^2,  Cm = 0.6 - 0.4 M1/M2,  M2,min = Pu (15 + 0.03 h) mm
       delta_ns = Cm / (1 - Pu / 0.75 Pc) >= 1
  3. demand/capacity along the ray to (Pu, delta2 M2, delta3 M3) on the design interaction surface,
     with phi Pn capped at 0.80 phi Po (tied)
Sway magnification (delta_s) is not applied here: use P-Delta analysis results for lateral combos.
"""
import math

from .section import RectColumn, rect_bars


def column_section(model, sec_name, fc=None, fy=None, n_face3=None, n_face2=None, tie_dia=None):
    """Build a RectColumn from the e2k frame section + CONCRETE SECTIONS data.

    PATTERN "R-a-b" is read as a bars in each face parallel to local 3 and b bars in each face
    parallel to local 2 (the order of the ETABS dialog). Override with n_face3 / n_face2 if the
    drawing differs. Bar centre = clear cover to ties + tie diameter + bar diameter / 2.
    """
    M = model.M
    fs = M["fsec"][sec_name]
    rc = M["rcsec"].get(sec_name, {})
    if fs.get("shape") != "Concrete Rectangular":
        raise NotImplementedError(f"{sec_name}: {fs.get('shape')}")
    mat = M["mats"][fs["mat"]]
    fc = fc or mat["FC"]
    fy = fy or M["mats"][rc["longbarmaterial"]]["FY"]
    if n_face3 is None or n_face2 is None:
        m = rc.get("pattern", "R-3-3").split("-")
        n_face3, n_face2 = n_face3 or int(m[1]), n_face2 or int(m[2])
    Ab = rc["longbararea"]
    db = math.sqrt(4 * Ab / math.pi)
    tie = tie_dia or math.sqrt(4 * rc.get("confinebararea", 63.6) / math.pi)
    cc = rc.get("cover", 40.0) + tie + db / 2
    sec = RectColumn(D=fs["D"], B=fs["B"], fc=fc, fy=fy, bars=rect_bars(fs["D"], fs["B"], n_face3, n_face2, Ab, cc),
                     tied=rc.get("transreinf", "TIES").upper() != "SPIRAL")
    sec.Ec = mat["E"] if fc == mat["FC"] else 4700 * math.sqrt(fc)
    sec.name, sec.n_face3, sec.n_face2, sec.bar_dia = sec_name, n_face3, n_face2, db
    return sec


def magnifier(Pu, Ma, Mb, h, Ig, Ec, lu, beta_dns, k=1.0, sway=True):
    """Nonsway magnification about one axis. Pu N, end moments N-mm (same sign = single curvature),
    h = section dimension in the bending direction (mm), lu mm. Returns dict."""
    M2 = max(abs(Ma), abs(Mb))
    M1 = min(abs(Ma), abs(Mb))
    ratio = 0.0 if M2 < 1e-9 else -M1 / M2 * (1 if Ma * Mb >= 0 else -1)   # ACI sign: negative = single curvature
    r = 0.3 * h
    slender_limit = 22.0 if sway else min(34 + 12 * ratio, 40.0)      # 6.2.5 (a) sway, (b) nonsway
    kl_r = k * lu / r
    M2min = Pu * (15 + 0.03 * h)
    out = dict(kl_r=kl_r, limit=slender_limit, M2=M2, M2min=M2min, delta=1.0, Cm=1.0, Pc=float("inf"), stable=True)
    Mdes = max(M2, M2min) if Pu > 0 else M2
    if Pu > 0 and kl_r > slender_limit:
        EI = 0.4 * Ec * Ig / (1 + beta_dns)
        Pc = math.pi ** 2 * EI / (k * lu) ** 2
        Cm = 1.0 if M2 < M2min else 0.6 - 0.4 * ratio
        if Pu >= 0.75 * Pc:
            out.update(Pc=Pc, Cm=Cm, delta=float("inf"), stable=False)
        else:
            out.update(Pc=Pc, Cm=Cm, delta=max(1.0, Cm / (1 - Pu / (0.75 * Pc))))
    out["Mc"] = Mdes * out["delta"]
    return out


def check_columns(model, combos=None, fc=None, fy=None, lu_deduct=0.0, sway=True, sections=None, **sec_kw):
    """Check every column. combos: {name: {pattern: factor}} (default: combinations in the e2k).
    fc, fy in MPa override the model materials. lu_deduct (m) is subtracted from the story height
    (e.g. slab thickness). sway: True/False for the whole building, or {story: bool} (e.g. from
    Model.stability(): Q > 0.05); sway stories use the k lu / r <= 22 limit of 6.2.5(a) and always get
    the along-the-length check of 6.6.4.6.4, nonsway stories use 34 + 12 M1/M2 <= 40. Returns a list of rows sorted by descending D/C (governing combo per column)."""
    combos = combos or model.M["combos"]
    if not combos:
        raise ValueError("no load combinations in the model: pass combos={name: {pattern: factor}}")
    pats = sorted({p for c in combos.values() for p in c})
    ff = model.frame_forces([p for p in pats if model.loads.get(p)])
    dead = {p for p, (typ, _) in model.M["patterns"].items() if "dead" in typ.lower()}
    secs = dict(sections or {})
    rows = []
    keys = [k for k, v in next(iter(ff.values())).items() if v["kind"] == "col"]
    for key in keys:
        f0 = next(iter(ff.values()))[key]
        name = f0["sec"]
        if name not in secs:
            secs[name] = column_section(model, name, fc=fc, fy=fy, **sec_kw)
        s = secs[name]
        lu = max(f0["L"] - lu_deduct, 0.1) * 1e3
        sw = sway.get(key[1], True) if isinstance(sway, dict) else bool(sway)
        best = None
        for cname, fac in combos.items():
            def comb(q):
                return sum(f * ff[p][key][q] for p, f in fac.items() if p in ff)
            Pu = comb("P") * 1e3
            Pd = sum(f * ff[p][key]["P"] for p, f in fac.items() if p in ff and p in dead) * 1e3
            bd = min(max(Pd / Pu, 0.0), 1.0) if Pu > 0 else 0.0
            m3 = magnifier(Pu, comb("M3_i") * 1e6, comb("M3_j") * 1e6, s.D, s.B * s.D ** 3 / 12, s.Ec, lu, bd, sway=sw)
            m2 = magnifier(Pu, comb("M2_i") * 1e6, comb("M2_j") * 1e6, s.B, s.D * s.B ** 3 / 12, s.Ec, lu, bd, sway=sw)
            stable = m2["stable"] and m3["stable"]
            dc = s.capacity_ratio(Pu / 1e3, m2["Mc"] / 1e6, m3["Mc"] / 1e6) if stable else float("inf")
            row = dict(column=key[0], story=key[1], section=name, combo=cname, lu_m=lu / 1e3, Pu_kN=Pu / 1e3,
                       M2_kNm=m2["M2"] / 1e6, M3_kNm=m3["M2"] / 1e6, kl_r2=m2["kl_r"], kl_r3=m3["kl_r"],
                       delta2=m2["delta"], delta3=m3["delta"], Mc2_kNm=m2["Mc"] / 1e6, Mc3_kNm=m3["Mc"] / 1e6,
                       beta_dns=bd, Pu_over_phiPnmax=Pu / s.phiPn_max, DC=dc,
                       DC_no_slender=s.capacity_ratio(Pu / 1e3, max(m2["M2"], m2["M2min"]) / 1e6,
                                                      max(m3["M2"], m3["M2min"]) / 1e6))
            if best is None or row["DC"] > best["DC"]:
                best = row
        rows.append(best)
    rows.sort(key=lambda r: -r["DC"])
    return rows, secs


def format_columns(rows, secs, top=None):
    o = []
    for s in secs.values():
        o.append(f"{s.name}: {s.B:.0f} x {s.D:.0f} mm, f'c {s.fc:.1f} MPa, fy {s.fy:.0f} MPa, {len(s.bars)} bars d{s.bar_dia:.0f} "
                 f"(rho {100 * s.rho:.2f} %; {s.n_face3} per face along axis 3, {s.n_face2} per face along axis 2), "
                 f"phiPn,max {s.phiPn_max / 1e3:.0f} kN")
    o.append(f"\n{'column':<7}{'story':<14}{'combo':<9}{'Pu kN':>8}{'M2':>7}{'M3':>7}{'kl/r2':>7}{'d2':>6}{'d3':>6}"
             f"{'Mc2':>7}{'Mc3':>7}{'Pu/Pmax':>9}{'D/C':>7}")
    for r in rows[:top]:
        o.append(f"{r['column']:<7}{r['story']:<14}{r['combo']:<9}{r['Pu_kN']:>8.0f}{r['M2_kNm']:>7.0f}{r['M3_kNm']:>7.0f}"
                 f"{r['kl_r2']:>7.1f}{r['delta2']:>6.2f}{r['delta3']:>6.2f}{r['Mc2_kNm']:>7.0f}{r['Mc3_kNm']:>7.0f}"
                 f"{r['Pu_over_phiPnmax']:>9.2f}{r['DC']:>7.2f}{'  FAIL' if r['DC'] > 1 else ''}")
    return "\n".join(o)
