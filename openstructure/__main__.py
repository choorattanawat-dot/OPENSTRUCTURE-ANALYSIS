"""Command line:  python -m openstructure MODEL.e2k [ETABS_TABLES.e2k] [--modes 12] [--out results.json]"""
import argparse
import json

from . import Model, __version__
from .verify import verify, format_report


def main():
    ap = argparse.ArgumentParser(prog="openstructure", description=f"openstructure {__version__}")
    ap.add_argument("e2k", help="ETABS model text export (.e2k)")
    ap.add_argument("tables", nargs="?", help="ETABS result tables exported as text (optional, for verification)")
    ap.add_argument("--modes", type=int, default=12)
    ap.add_argument("--stability", action="store_true",
                    help="story drift, stability index Q and P-Delta amplification (1.2D + 1.0L, F ~ m*z)")
    ap.add_argument("--out", help="write results to this JSON file")
    a = ap.parse_args()
    m = Model.from_e2k(a.e2k)
    info = m.info()
    print(f"nodes {info['nodes']}, shells {info['shells']}, frame segments {info['frame_segments']}, "
          f"supports {info['supports']}, file units {'-'.join(info['units_in_file'][:2])}")
    for w in info["warnings"]:
        print("[!] " + w)
    if info["not_modelled"]:
        print("not modelled: " + "; ".join(info["not_modelled"]))
    print(f"total lateral mass {m.total_mass:.1f} t")
    if a.tables:
        v = verify(m, a.tables, a.modes)
        print(format_report(v))
        res = {"rows": v["rows"], "participation": v["participation"]}
    else:
        g, mod = m.gravity(), m.modal(a.modes)
        for p, r in g.items():
            print(f"{p:<8} applied {r['applied_kN']:.1f} kN, base Rz {r['base_Rz_kN']:.1f} kN")
        for i, t in enumerate(mod["periods_s"], 1):
            print(f"T{i} = {t:.4f} s   UX {mod['UX'][i-1]:.1f} %  UY {mod['UY'][i-1]:.1f} %  RZ {mod['RZ'][i-1]:.1f} %")
        res = {"gravity": g, "modal": mod}
    if a.stability:
        from .stability import format_stability
        st = m.stability()
        print(format_stability(st))
        res["stability"] = st
    if a.out:
        res.update(info=info, mass_total_t=m.total_mass, story_mass_t=m.story_mass())
        json.dump(res, open(a.out, "w"), indent=1, default=float)


if __name__ == "__main__":
    main()
