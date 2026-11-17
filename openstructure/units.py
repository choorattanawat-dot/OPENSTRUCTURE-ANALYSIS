"""Unit factors for ETABS text exports.

Internally the parser normalises every ETABS file to N and mm (the convention the
mesh builder was written and regression-tested in). The OpenSees model itself is
built in kN, m, s, t (mass in t = kN*s^2/m).
"""
G_ACC = 9.80665  # m/s^2

# force unit -> N
FORCE_TO_N = {
    "N": 1.0, "KN": 1e3, "MN": 1e6,
    "KGF": G_ACC, "TONF": 1e3 * G_ACC,
    "LB": 4.4482216152605, "KIP": 4448.2216152605,
}
# length unit -> mm
LENGTH_TO_MM = {"MM": 1.0, "CM": 10.0, "M": 1000.0, "IN": 25.4, "FT": 304.8}


def factors(units):
    """units = ["KGF", "M", "C"] as written on the CONTROLS/UNITS line.
    Returns (force factor to N, length factor to mm)."""
    f, l = units[0].upper(), units[1].upper()
    if f not in FORCE_TO_N:
        raise ValueError(f"Unsupported ETABS force unit {units[0]!r}")
    if l not in LENGTH_TO_MM:
        raise ValueError(f"Unsupported ETABS length unit {units[1]!r}")
    return FORCE_TO_N[f], LENGTH_TO_MM[l]
