"""Uniaxial P–M check for a short rectangular column.

Axial force is positive in compression. Moment is positive when it compresses
the positive-y face (y measured from the geometric centroid). The check builds
the nominal interaction curve, cuts it at Pn,max, applies φ, and reports the
radial utilization of (Pu, Mu).

This is a section check. Slenderness and P-delta are not included.
"""

from __future__ import annotations

import math
from dataclasses import dataclass

from openstructure.design.codes import EPS_CU, ES_MPA, DesignCode
from openstructure.design.report import CalcReport, fmt
from openstructure.design.section import (
    NominalResult,
    bars_are_symmetric,
    c_grid,
    layer_ratio,
    nominal_at_c,
    pure_compression,
    pure_tension,
)
from openstructure.errors import OpenStructureError
from openstructure.units import kn_to_n, knm_to_nmm, n_to_kn, nmm_to_knm

Bar = tuple[float, float]


@dataclass
class ColumnResult:
    status: str
    report: CalcReport
    utilization: float
    phi_pn_max_kN: float
    po_kN: float
    inside: bool


def perimeter_bars(
    b_mm: float,
    h_mm: float,
    n_along_b: int,
    n_along_h: int,
    diameter_mm: float,
    cover_to_center_mm: float,
) -> list[Bar]:
    """Bars around a rectangular perimeter.

    ``h_mm`` is the bending depth (y). ``n_along_h`` counts bars through that
    depth, including the corner bars. Intermediate rows have two side bars.
    Extreme rows have ``n_along_b`` bars. Only the y coordinate affects uniaxial
    bending about the axis perpendicular to h.
    """
    if n_along_b < 2 or n_along_h < 2:
        raise OpenStructureError("Perimeter layout needs at least 2 bars along each side.")
    if cover_to_center_mm <= 0 or cover_to_center_mm >= h_mm / 2:
        raise OpenStructureError("cover_to_center_mm must lie inside the section.")
    area = math.pi * diameter_mm**2 / 4.0
    y_top = h_mm / 2.0 - cover_to_center_mm
    y_bot = -y_top
    bars: list[Bar] = []
    for _ in range(n_along_b):
        bars.append((area, y_top))
        bars.append((area, y_bot))
    if n_along_h > 2:
        for index in range(1, n_along_h - 1):
            y_mm = y_bot + (y_top - y_bot) * index / (n_along_h - 1)
            bars.append((area, y_mm))
            bars.append((area, y_mm))
    return bars


def check_column_uniaxial(
    pu_kN: float,
    mu_kNm: float,
    b_mm: float,
    h_mm: float,
    bars: list[Bar],
    fc_mpa: float,
    fy_mpa: float,
    code: DesignCode,
    tied: bool = True,
) -> ColumnResult:
    if min(b_mm, h_mm, fc_mpa, fy_mpa) <= 0:
        raise OpenStructureError("Column check: b, h, f'c, and fy must be positive.")
    if not bars:
        raise OpenStructureError("Column check: provide at least one bar.")
    factors = code.factors()
    beta1 = code.beta1(fc_mpa)
    symmetric = bars_are_symmetric(bars)
    ratio = layer_ratio(h_mm, bars)
    po = pure_compression(b_mm, h_mm, bars, fc_mpa, fy_mpa)
    pn_cap = code.pn_max_factor(tied) * po.pn_n
    phi_c = code.phi_compression(tied)
    phi_pn_max = phi_c * pn_cap
    ast = sum(area for area, _y in bars)
    ag = b_mm * h_mm
    rho_g = ast / ag
    pb_n = _balanced_pn(b_mm, h_mm, bars, fc_mpa, fy_mpa, beta1)

    report = CalcReport(
        title="RC column, axial force and uniaxial moment",
        status="OK",
        code=factors.code_name,
        assumptions=list(factors.assumptions),
    )
    if not symmetric:
        report.assumptions.insert(
            0,
            "Reinforcement is not symmetric. Moment is about the geometric centroid; "
            "v0.1 does not move the origin to the plastic centroid.",
        )
    report.add(
        "Inputs",
        "Pu positive in compression. Mu positive when it compresses the +y face.",
        (
            f"Pu = {fmt(pu_kN, 2)} kN, Mu = {fmt(mu_kNm, 2)} kN·m, "
            f"b = {fmt(b_mm, 1)} mm, h = {fmt(h_mm, 1)} mm, "
            f"f'c = {fmt(fc_mpa, 2)} MPa, fy = {fmt(fy_mpa, 1)} MPa, "
            f"{'tied' if tied else 'spiral'}"
        ),
        (
            f"Ast = {fmt(ast, 1)} mm², Ag = {fmt(ag, 1)} mm², ρg = {fmt(rho_g, 4)}, "
            f"β1 = {fmt(beta1, 4)}, bars symmetric = {symmetric}"
        ),
    )
    report.add(
        "Concentric nominal strength and the axial cap",
        "Po = 0.85 f'c (Ag − Ast) + fy Ast;  Pn,max = 0.80 Po (tied) or 0.85 Po (spiral);  "
        "φPn,max uses the compression-controlled φ",
        (
            f"0.85 f'c (Ag − Ast) = {fmt(n_to_kn(0.85 * fc_mpa * (ag - ast)), 2)} kN; "
            f"fy Ast = {fmt(n_to_kn(fy_mpa * ast), 2)} kN"
        ),
        (
            f"Po = {fmt(n_to_kn(po.pn_n), 2)} kN; Pn,max = {fmt(n_to_kn(pn_cap), 2)} kN; "
            f"φ = {fmt(phi_c, 2)}; φPn,max = {fmt(n_to_kn(phi_pn_max), 2)} kN"
        ),
        reference="ACI 318-19 §22.4.2. The same Po and 0.80/0.85 factors are in the วสท. 1008 notes.",
    )

    diagram = _design_diagram(
        b_mm, h_mm, bars, fc_mpa, fy_mpa, beta1, code, tied, pn_cap, symmetric, ratio, pb_n
    )
    pu_n = kn_to_n(pu_kN)
    mu_nmm = knm_to_nmm(mu_kNm)
    scale = _ray_scale(pu_n, mu_nmm, diagram)
    if scale is None:
        utilization = float("inf")
        report.add(
            "Demand ray does not meet the design diagram",
            "Scale s such that s·(Pu, Mu) lies on the φPn–φMn diagram",
            f"Pu = {fmt(pu_kN, 2)} kN, Mu = {fmt(mu_kNm, 2)} kN·m",
            "No intersection was found.",
        )
    else:
        utilization = 0.0 if math.isinf(scale) else 1.0 / scale
        cap_text = ""
        if pu_kN > 0:
            axial_util = pu_n / phi_pn_max
            if axial_util > utilization:
                utilization = axial_util
                cap_text = f" Axial cap utilization Pu/φPn,max = {fmt(axial_util, 3)} governs."
        report.add(
            "Radial utilization",
            "s is the factor on (Pu, Mu) that lands on the design diagram. Utilization = 1/s.",
            f"Pu = {fmt(pu_kN, 2)} kN, Mu = {fmt(mu_kNm, 2)} kN·m, {len(diagram)} diagram points",
            (
                f"s = {fmt(scale, 3)}; utilization = {fmt(utilization, 3)}; "
                f"capacity along this ray ≈ {fmt(pu_kN * scale, 2)} kN, "
                f"{fmt(mu_kNm * scale, 2)} kN·m.{cap_text}"
            ),
            note="A value above 1 means (Pu, Mu) is outside the design diagram.",
        )

    rho_ok = 0.01 - 1.0e-9 <= rho_g <= 0.08 + 1.0e-9
    report.add(
        "Longitudinal reinforcement limits",
        "0.01 Ag ≤ Ast ≤ 0.08 Ag",
        f"Ast = {fmt(ast, 1)} mm², Ag = {fmt(ag, 1)} mm²",
        f"ρg = {fmt(rho_g, 4)} → {'inside 0.01 to 0.08' if rho_ok else 'outside 0.01 to 0.08'}",
        reference="ACI 318-19 §10.6.1.1. วสท. 1008 teaching notes state the same limits.",
    )
    _add_diagram_table(report, diagram)

    inside = utilization <= 1.001 and rho_ok
    status = "OK" if inside else "NG"
    if scale is None:
        status = "NG"
    report.status = status
    report.summary = [
        f"Po = {fmt(n_to_kn(po.pn_n), 2)} kN; φPn,max = {fmt(n_to_kn(phi_pn_max), 2)} kN",
        f"Pu = {fmt(pu_kN, 2)} kN, Mu = {fmt(mu_kNm, 2)} kN·m",
        f"Utilization = {fmt(utilization, 3)}",
        f"ρg = {fmt(rho_g, 4)} ({'OK' if rho_ok else 'outside 0.01–0.08'})",
        "Short-column section check only. Slenderness is not included.",
    ]
    return ColumnResult(
        status=status,
        report=report,
        utilization=utilization,
        phi_pn_max_kN=n_to_kn(phi_pn_max),
        po_kN=n_to_kn(po.pn_n),
        inside=inside,
    )


def _balanced_pn(
    b_mm: float,
    h_mm: float,
    bars: list[Bar],
    fc_mpa: float,
    fy_mpa: float,
    beta1: float,
) -> float:
    """Smaller balanced axial force of the two faces. Used only for the EIT φ ramp."""
    eps_y = fy_mpa / ES_MPA
    values: list[float] = []
    for sign in (1, -1):
        d_t = max(h_mm / 2.0 - sign * y for _area, y in bars)
        if d_t <= 0:
            continue
        c_b = EPS_CU / (EPS_CU + eps_y) * d_t
        point = nominal_at_c(b_mm, h_mm, bars, c_b, fc_mpa, fy_mpa, beta1, sign)
        values.append(point.pn_n)
    return min(values) if values else 0.0


def _design_diagram(
    b_mm: float,
    h_mm: float,
    bars: list[Bar],
    fc_mpa: float,
    fy_mpa: float,
    beta1: float,
    code: DesignCode,
    tied: bool,
    pn_cap: float,
    symmetric: bool,
    ratio: float,
    pb_n: float,
) -> list[tuple[float, float]]:
    nominal_points: list[NominalResult] = [pure_tension(bars, fy_mpa)]
    for sign in (1, -1):
        for c_mm in c_grid(h_mm):
            point = nominal_at_c(b_mm, h_mm, bars, c_mm, fc_mpa, fy_mpa, beta1, sign)
            if point.pn_n <= pn_cap:
                nominal_points.append(point)
        capped = _point_at_axial_cap(
            b_mm, h_mm, bars, fc_mpa, fy_mpa, beta1, sign, pn_cap
        )
        if capped is not None:
            nominal_points.append(capped)

    design: list[tuple[float, float]] = []
    for point in nominal_points:
        phi = code.phi_column(
            eps_t=point.eps_t,
            fy_mpa=fy_mpa,
            tied=tied,
            pn_n=max(point.pn_n, 0.0) if point.pn_n > 0 else point.pn_n,
            fc_mpa=fc_mpa,
            ag_mm2=b_mm * h_mm,
            pb_n=pb_n,
            symmetric=symmetric,
            layer_ratio=ratio,
        )
        design.append((phi * point.pn_n, phi * point.mn_nmm))
    design.sort(key=lambda item: math.atan2(item[1], item[0]))
    return design


def _point_at_axial_cap(
    b_mm: float,
    h_mm: float,
    bars: list[Bar],
    fc_mpa: float,
    fy_mpa: float,
    beta1: float,
    sign: int,
    pn_cap: float,
) -> NominalResult | None:
    low = h_mm / 400.0
    high = h_mm * 40.0
    low_point = nominal_at_c(b_mm, h_mm, bars, low, fc_mpa, fy_mpa, beta1, sign)
    high_point = nominal_at_c(b_mm, h_mm, bars, high, fc_mpa, fy_mpa, beta1, sign)
    if low_point.pn_n > pn_cap or high_point.pn_n < pn_cap:
        return None
    for _ in range(50):
        mid = 0.5 * (low + high)
        point = nominal_at_c(b_mm, h_mm, bars, mid, fc_mpa, fy_mpa, beta1, sign)
        if point.pn_n > pn_cap:
            high = mid
        else:
            low = mid
    return nominal_at_c(b_mm, h_mm, bars, 0.5 * (low + high), fc_mpa, fy_mpa, beta1, sign)


def _ray_scale(
    pu_n: float,
    mu_nmm: float,
    diagram: list[tuple[float, float]],
) -> float | None:
    """Factor s > 0 such that s·(Pu, Mu) lies on the closed design diagram."""
    if abs(pu_n) < 1.0e-6 and abs(mu_nmm) < 1.0e-3:
        return float("inf")  # zero demand: utilization 0
    best: float | None = None
    count = len(diagram)
    for index in range(count):
        p1, m1 = diagram[index]
        p2, m2 = diagram[(index + 1) % count]
        scale = _segment_scale(pu_n, mu_nmm, p1, m1, p2, m2)
        if scale is not None and scale > 1.0e-8 and (best is None or scale < best):
            best = scale
    return best


def _segment_scale(
    pu: float,
    mu: float,
    p1: float,
    m1: float,
    p2: float,
    m2: float,
) -> float | None:
    dp = p2 - p1
    dm = m2 - m1
    det = -pu * dm + dp * mu
    if abs(det) < 1.0e-8:
        return None
    scale = (-p1 * dm + dp * m1) / det
    position = (pu * m1 - p1 * mu) / det
    if scale > 0 and -1.0e-8 <= position <= 1.0 + 1.0e-8:
        return scale
    return None


def _add_diagram_table(report: CalcReport, diagram: list[tuple[float, float]]) -> None:
    if not diagram:
        return
    # Show a short sample spread around the diagram, in kN and kN·m.
    step = max(1, len(diagram) // 8)
    rows = []
    for index in range(0, len(diagram), step):
        phi_p, phi_m = diagram[index]
        rows.append(f"φPn = {fmt(n_to_kn(phi_p), 1)} kN, φMn = {fmt(nmm_to_knm(phi_m), 2)} kN·m")
    report.add(
        "Design-strength diagram, sample points",
        "φPn and φMn after the Pn,max cut. Positive moment compresses the +y face.",
        f"{len(diagram)} points were used for the ray test. Eight of them:",
        "; ".join(rows[:8]),
    )
