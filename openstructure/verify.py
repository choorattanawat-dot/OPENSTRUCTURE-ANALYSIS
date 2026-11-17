"""Compare an openstructure Model with ETABS result tables: ETABS | OPENSEES | % DIFF."""
from .etabs_tables import etabs_reference


def _row(label, e, o, unit=""):
    d = 100 * (o - e) / e if e else float("nan")
    return dict(item=label, etabs=e, opensees=o, diff_pct=d, unit=unit)


def verify(model, tables_path, n_modes=12):
    ref = etabs_reference(tables_path)
    rows = []
    g = model.gravity()
    for pat, r in g.items():
        if pat in ref["base_FZ_kN"]:
            rows.append(_row(f"Base FZ {pat}", ref["base_FZ_kN"][pat], r["applied_kN"], "kN"))
    mod = model.modal(n_modes)
    for i, (te, to) in enumerate(zip(ref["periods_s"], mod["periods_s"]), 1):
        rows.append(_row(f"T{i}", te, to, "s"))
    part = []
    for i in range(min(n_modes, len(ref["periods_s"]))):
        part.append(dict(mode=i + 1, etabs=(ref["UX"][i], ref["UY"][i], ref["RZ"][i]),
                         opensees=(mod["UX"][i] / 100, mod["UY"][i] / 100, mod["RZ"][i] / 100)))
    return dict(rows=rows, participation=part, gravity=g, modal=mod, reference=ref)


def format_report(v):
    out = [f"{'ITEM':<22}{'ETABS':>14}{'OPENSEES':>14}{'DIFF %':>9}"]
    for r in v["rows"]:
        out.append(f"{r['item']:<22}{r['etabs']:>14.4f}{r['opensees']:>14.4f}{r['diff_pct']:>9.2f}  {r['unit']}")
    out.append("")
    out.append("Modal participating mass ratio (UX, UY, RZ)   ETABS | OPENSEES")
    for p in v["participation"]:
        e, o = p["etabs"], p["opensees"]
        out.append(f"  mode {p['mode']:>2}:  {e[0]:.3f} {e[1]:.3f} {e[2]:.3f}  |  {o[0]:.3f} {o[1]:.3f} {o[2]:.3f}")
    return "\n".join(out)
