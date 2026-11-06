"""Rectangular beam flexure: singly and doubly reinforced.

Proportioning is algebraic so the report can show every substitution. The
resulting reinforcement is then checked with the shared Whitney-block solver.
Moment is factored Mu. A negative value is treated as hogging of the same
rectangular section and is reported as such.
"""

from __future__ import annotations

import math
from dataclasses import dataclass

from openstructure.design.codes import EPS_CU, ES_MPA, DesignCode
from openstructure.design.report import CalcReport, fmt
from openstructure.design.section import nominal_at_c, solve_neutral_axis
from openstructure.errors import OpenStructureError
from openstructure.units import knm_to_nmm, n_to_kn, nmm_to_knm


@dataclass
class FlexureResult:
    status: str
    report: CalcReport
    as_mm2: float
    as_prime_mm2: float
    phi_mn_kNm: float
    utilization: float
    doubly_reinforced: bool


def design_beam_flexure(
    mu_kNm: float,
    b_mm: float,
    d_mm: float,
    h_mm: float,
    fc_mpa: float,
    fy_mpa: float,
    code: DesignCode,
    d_prime_mm: float | None = None,
) -> FlexureResult:
    if min(b_mm, d_mm, h_mm, fc_mpa, fy_mpa) <= 0:
        raise OpenStructureError("Beam flexure: b, d, h, f'c, and fy must be positive.")
    if d_mm >= h_mm:
        raise OpenStructureError("Beam flexure: d must be less than h.")
    sign_note = ""
    if mu_kNm < 0:
        sign_note = "Mu was negative (hogging). The same rectangular section is designed for |Mu|."
        mu_kNm = abs(mu_kNm)

    factors = code.factors()
    beta1 = code.beta1(fc_mpa)
    c_over_d = code.c_over_d_limit(fy_mpa)
    phi = factors.phi_flexure
    mu_nmm = knm_to_nmm(mu_kNm)
    report = CalcReport(
        title="RC beam flexure",
        status="OK",
        code=factors.code_name,
        assumptions=list(factors.assumptions),
    )
    if sign_note:
        report.assumptions.insert(0, sign_note)

    report.add(
        "Inputs",
        "Mu, b, d, h, f'c, fy, φ, β1",
        (
            f"Mu = {fmt(mu_kNm)} kN·m, b = {fmt(b_mm, 1)} mm, d = {fmt(d_mm, 1)} mm, "
            f"h = {fmt(h_mm, 1)} mm, f'c = {fmt(fc_mpa, 2)} MPa, fy = {fmt(fy_mpa, 1)} MPa"
        ),
        f"φ flexure = {fmt(phi, 3)}, β1 = {fmt(beta1, 4)}, c/d limit = {fmt(c_over_d, 4)}",
        reference=factors.code_name,
    )

    a_max = beta1 * c_over_d * d_mm
    cc_max = 0.85 * fc_mpa * b_mm * a_max
    as_limit = cc_max / fy_mpa
    mn_limit = as_limit * fy_mpa * (d_mm - a_max / 2.0)
    phi_mn_limit = phi * mn_limit
    report.add(
        "Singly-reinforced limit at the code's maximum reinforcement",
        "a_max = β1 (c/d)_limit d;  Cc = 0.85 f'c b a;  Mn = As fy (d − a/2);  As = Cc / fy",
        (
            f"a_max = {fmt(beta1, 4)} × {fmt(c_over_d, 4)} × {fmt(d_mm, 1)} = {fmt(a_max, 2)} mm; "
            f"Cc = 0.85 × {fmt(fc_mpa, 2)} × {fmt(b_mm, 1)} × {fmt(a_max, 2)} = {fmt(cc_max, 1)} N"
        ),
        (
            f"As at the limit = {fmt(as_limit, 1)} mm²; "
            f"φMn,limit = {fmt(nmm_to_knm(phi_mn_limit), 2)} kN·m"
        ),
        reference="Whitney block. ACI uses εt = 0.005. EIT uses 0.75 ρb (assumption).",
    )

    doubly = mu_nmm > phi_mn_limit * (1.0 + 1.0e-9)
    as_prime = 0.0
    if not doubly:
        as_mm2, a_mm, c_mm, eps_t, rho = _singly(
            report, mu_nmm, phi, b_mm, d_mm, fc_mpa, fy_mpa
        )
    else:
        if d_prime_mm is None:
            raise OpenStructureError(
                "Mu exceeds the singly-reinforced limit. Provide d_prime_mm to design compression steel."
            )
        if d_prime_mm <= 0 or d_prime_mm >= d_mm:
            raise OpenStructureError("d' must be positive and less than d.")
        as_mm2, as_prime, a_mm, c_mm, eps_t = _doubly(
            report,
            mu_nmm,
            phi,
            b_mm,
            d_mm,
            d_prime_mm,
            fc_mpa,
            fy_mpa,
            beta1,
            c_over_d,
            a_max,
            cc_max,
            as_limit,
            mn_limit,
        )
        rho = as_mm2 / (b_mm * d_mm)

    as_min = code.as_min_beam_mm2(fc_mpa, fy_mpa, b_mm, d_mm)
    report.add(
        "Minimum flexural reinforcement",
        "As,min from the code plugin (see assumptions for the expression)",
        f"b = {fmt(b_mm, 1)} mm, d = {fmt(d_mm, 1)} mm, f'c = {fmt(fc_mpa, 2)} MPa, fy = {fmt(fy_mpa, 1)} MPa",
        f"As,min = {fmt(as_min, 1)} mm²; As required by strength = {fmt(as_mm2, 1)} mm²",
    )
    governed_by_min = as_mm2 < as_min
    if governed_by_min:
        as_mm2 = as_min
        report.add(
            "Minimum reinforcement governs",
            "As = As,min",
            f"As,min = {fmt(as_min, 1)} mm² > strength requirement",
            f"Provide As = {fmt(as_mm2, 1)} mm². φMn will exceed Mu.",
        )

    y_tension = h_mm / 2.0 - d_mm
    bars = [(as_mm2, y_tension)]
    if as_prime > 0 and d_prime_mm is not None:
        bars.append((as_prime, h_mm / 2.0 - d_prime_mm))
    checked = solve_neutral_axis(b_mm, h_mm, bars, fc_mpa, fy_mpa, beta1, compression_sign=1)
    phi_check = code.phi_flexure_for_strain(checked.eps_t, fy_mpa)
    if code.code_id == "ACI318-19":
        phi_check = code.phi_flexure_for_strain(checked.eps_t, fy_mpa)
        strain_ok = checked.eps_t >= 0.004 - 1.0e-6
        strain_note = "ACI 318-19 §9.3.3.1 requires εt ≥ 0.004 for a flexural member."
    else:
        eps_y = fy_mpa / ES_MPA
        strain_ok = checked.eps_t >= eps_y - 1.0e-6
        phi_check = factors.phi_flexure if strain_ok else factors.phi_flexure
        strain_note = (
            "EIT beam check (assumption): tension steel must yield (εt ≥ εy), "
            "which is the under-reinforced limit behind ρmax = 0.75 ρb. "
            "φ is not increased above the flexure value."
        )
    phi_mn = phi_check * checked.mn_nmm
    utilization = mu_nmm / phi_mn if phi_mn > 0 else float("inf")
    report.add(
        "Independent section check of the provided reinforcement",
        "Neutral axis solved so that Pn = 0 (bisection). Mn from the Whitney block.",
        (
            f"As = {fmt(as_mm2, 1)} mm² at y = {fmt(y_tension, 1)} mm from the centroid; "
            f"As' = {fmt(as_prime, 1)} mm²"
        ),
        (
            f"c = {fmt(checked.c_mm, 2)} mm, a = {fmt(checked.a_mm, 2)} mm, "
            f"εt = {fmt(checked.eps_t, 5)}, φ = {fmt(phi_check, 3)}, "
            f"φMn = {fmt(nmm_to_knm(phi_mn), 2)} kN·m, "
            f"Mu / φMn = {fmt(utilization, 3)}"
        ),
        note=strain_note,
    )

    status = "OK"
    if not strain_ok:
        status = "NG"
    elif utilization > 1.001:
        status = "NG"
    report.status = status
    report.summary = [
        f"As = {fmt(as_mm2, 1)} mm²" + (f", As' = {fmt(as_prime, 1)} mm²" if as_prime else ", singly reinforced"),
        f"φMn = {fmt(nmm_to_knm(phi_mn), 2)} kN·m for Mu = {fmt(mu_kNm, 2)} kN·m",
        f"Utilization Mu/φMn = {fmt(utilization, 3)}",
        f"ρ = As/(bd) = {fmt(as_mm2 / (b_mm * d_mm), 5)}",
        "Rectangular beam. No flange, no compression steel unless Mu requires it. "
        "Development, crack control, and deflection are not checked.",
    ]
    if not strain_ok:
        report.summary.insert(0, "NG: tensile strain is below the code limit for a beam.")
    del a_mm, c_mm, eps_t, rho  # the solver result is the one reported
    return FlexureResult(
        status=status,
        report=report,
        as_mm2=as_mm2,
        as_prime_mm2=as_prime,
        phi_mn_kNm=nmm_to_knm(phi_mn),
        utilization=utilization,
        doubly_reinforced=as_prime > 0,
    )


def _singly(
    report: CalcReport,
    mu_nmm: float,
    phi: float,
    b_mm: float,
    d_mm: float,
    fc_mpa: float,
    fy_mpa: float,
) -> tuple[float, float, float, float, float]:
    rn = mu_nmm / (phi * b_mm * d_mm**2)
    disc = 1.0 - 2.0 * rn / (0.85 * fc_mpa)
    if disc <= 0:
        raise OpenStructureError("Singly-reinforced discriminant is not positive. Increase the section.")
    rho = (0.85 * fc_mpa / fy_mpa) * (1.0 - math.sqrt(disc))
    as_mm2 = rho * b_mm * d_mm
    a_mm = as_mm2 * fy_mpa / (0.85 * fc_mpa * b_mm)
    # β1 is not required for a; c needs it. Recover β1 from a later in the solver.
    # εt needs c = a/β1. The caller recomputes with the solver; show a here.
    report.add(
        "Singly-reinforced steel for Mu",
        "Rn = Mu / (φ b d²);  ρ = (0.85 f'c / fy) [1 − √(1 − 2 Rn / (0.85 f'c))];  As = ρ b d",
        (
            f"Rn = {fmt(mu_nmm, 1)} / ({fmt(phi, 2)} × {fmt(b_mm, 1)} × {fmt(d_mm, 1)}²) "
            f"= {fmt(rn, 4)} MPa"
        ),
        (
            f"ρ = {fmt(rho, 5)}; As = {fmt(as_mm2, 1)} mm²; "
            f"a = As fy / (0.85 f'c b) = {fmt(a_mm, 2)} mm"
        ),
    )
    # c and εt are finalized by the section solver. Placeholders keep the return shape.
    return as_mm2, a_mm, 0.0, 0.0, rho


def _doubly(
    report: CalcReport,
    mu_nmm: float,
    phi: float,
    b_mm: float,
    d_mm: float,
    d_prime_mm: float,
    fc_mpa: float,
    fy_mpa: float,
    beta1: float,
    c_over_d: float,
    a_max: float,
    cc_max: float,
    as1: float,
    mn1: float,
) -> tuple[float, float, float, float, float]:
    c_mm = c_over_d * d_mm
    eps_s = EPS_CU * (c_mm - d_prime_mm) / c_mm
    fs_prime = min(fy_mpa, max(0.0, eps_s * ES_MPA))
    report.add(
        "Compression reinforcement is required",
        "Set c to the code limit. fs' = Es × 0.003 (1 − d'/c), not above fy.",
        (
            f"c = {fmt(c_over_d, 4)} × {fmt(d_mm, 1)} = {fmt(c_mm, 2)} mm; "
            f"d' = {fmt(d_prime_mm, 1)} mm; εs' = {fmt(eps_s, 5)}"
        ),
        f"fs' = {fmt(fs_prime, 2)} MPa; 0.85 f'c = {fmt(0.85 * fc_mpa, 2)} MPa",
    )
    net = fs_prime - 0.85 * fc_mpa
    if d_prime_mm >= a_max or net <= 0:
        raise OpenStructureError(
            "Compression steel at this d' does not add compression capacity "
            f"(d' = {d_prime_mm:.1f} mm, a_max = {a_max:.1f} mm, fs' − 0.85 f'c = {net:.2f} MPa). "
            "Increase h or reduce d'."
        )
    mu2 = mu_nmm - phi * mn1
    mn2 = mu2 / phi
    as2 = mn2 / (fy_mpa * (d_mm - d_prime_mm))
    as_prime = as2 * fy_mpa / net
    as_total = as1 + as2
    eps_t = EPS_CU * (d_mm - c_mm) / c_mm
    # Equilibrium shown with forces in N.
    tension = as_total * fy_mpa
    compression = cc_max + as_prime * net
    report.add(
        "Doubly-reinforced areas",
        "Mn2 = Mu/φ − Mn1;  As2 = Mn2 / [fy (d − d')];  As' = As2 fy / (fs' − 0.85 f'c);  As = As1 + As2",
        (
            f"Mn1 = {fmt(nmm_to_knm(mn1), 2)} kN·m; Mu2/φ = {fmt(nmm_to_knm(mn2), 2)} kN·m; "
            f"lever arm d − d' = {fmt(d_mm - d_prime_mm, 1)} mm"
        ),
        (
            f"As2 = {fmt(as2, 1)} mm²; As' = {fmt(as_prime, 1)} mm²; "
            f"As = {fmt(as_total, 1)} mm²; εt = {fmt(eps_t, 5)}"
        ),
        note=(
            f"Force check: tension As fy = {fmt(n_to_kn(tension), 2)} kN; "
            f"concrete plus net compression steel = {fmt(n_to_kn(compression), 2)} kN."
        ),
    )
    if abs(tension - compression) / tension > 0.02:
        raise OpenStructureError("Doubly-reinforced force balance is off by more than 2%.")
    del b_mm, beta1
    return as_total, as_prime, a_max, c_mm, eps_t


def beam_bars_for_capacity(
    h_mm: float,
    d_mm: float,
    as_mm2: float,
    as_prime_mm2: float = 0.0,
    d_prime_mm: float | None = None,
) -> list[tuple[float, float]]:
    """Bar groups for a sagging beam, +y toward the compression (top) face."""
    bars = [(as_mm2, h_mm / 2.0 - d_mm)]
    if as_prime_mm2 > 0:
        if d_prime_mm is None:
            raise OpenStructureError("d_prime_mm is required when As' is given.")
        bars.append((as_prime_mm2, h_mm / 2.0 - d_prime_mm))
    return bars


def beam_capacity_phi_mn_nmm(
    b_mm: float,
    h_mm: float,
    bars: list[tuple[float, float]],
    fc_mpa: float,
    fy_mpa: float,
    code: DesignCode,
) -> tuple[float, float, float]:
    """Return (φMn in N·mm, φ, εt) for pure flexure, compression on +y."""
    point = solve_neutral_axis(
        b_mm, h_mm, bars, fc_mpa, fy_mpa, code.beta1(fc_mpa), compression_sign=1
    )
    if code.code_id == "EIT1008":
        phi = code.factors().phi_flexure
    else:
        phi = code.phi_flexure_for_strain(point.eps_t, fy_mpa)
    return phi * point.mn_nmm, phi, point.eps_t


# nominal_at_c is re-exported for tests that want a fixed c rather than Pn = 0.
__all__ = [
    "FlexureResult",
    "beam_bars_for_capacity",
    "beam_capacity_phi_mn_nmm",
    "design_beam_flexure",
    "nominal_at_c",
]
