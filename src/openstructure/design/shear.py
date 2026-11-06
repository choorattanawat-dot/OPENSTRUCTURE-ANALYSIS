"""Beam shear reinforcement (vertical stirrups).

Vu is the factored shear at the station the caller chooses. v0.1 does not
move the critical section to a distance d from the support face; the pipeline
can do that and pass the resulting Vu. Axial force is taken as zero.
"""

from __future__ import annotations

import math
from dataclasses import dataclass

from openstructure.design.codes import DesignCode
from openstructure.design.report import CalcReport, fmt
from openstructure.errors import OpenStructureError
from openstructure.units import kn_to_n, n_to_kn


@dataclass
class ShearResult:
    status: str
    report: CalcReport
    stirrups_required: bool
    av_mm2: float
    s_provided_mm: float
    phi_vn_kN: float
    utilization: float


def design_beam_shear(
    vu_kN: float,
    b_mm: float,
    d_mm: float,
    fc_mpa: float,
    fyt_mpa: float,
    code: DesignCode,
    stirrup_diameter_mm: float = 9.0,
    legs: int = 2,
    spacing_step_mm: float = 25.0,
) -> ShearResult:
    if min(b_mm, d_mm, fc_mpa, fyt_mpa, stirrup_diameter_mm) <= 0:
        raise OpenStructureError("Beam shear: b, d, f'c, fyt, and the bar diameter must be positive.")
    if legs < 1:
        raise OpenStructureError("A stirrup needs at least one leg.")
    if vu_kN < 0:
        vu_kN = abs(vu_kN)
    vu_n = kn_to_n(vu_kN)
    factors = code.factors()
    fyt_use = code.fyt_for_design(fyt_mpa)
    av = legs * math.pi * stirrup_diameter_mm**2 / 4.0
    # First estimate Vc as if minimum reinforcement will be provided. Refined below.
    vc = code.vc_n(fc_mpa, b_mm, d_mm, av_meets_min=True)
    phi = factors.phi_shear
    phi_vc = phi * vc
    trigger = code.av_min_trigger_n(fc_mpa, b_mm, d_mm)
    av_per_s = code.av_min_per_mm(fc_mpa, b_mm, fyt_mpa)

    report = CalcReport(
        title="RC beam shear",
        status="OK",
        code=factors.code_name,
        assumptions=list(factors.assumptions),
    )
    report.add(
        "Inputs",
        "Vu ≤ φ (Vc + Vs);  Vs = Av fyt d / s",
        (
            f"Vu = {fmt(vu_kN, 2)} kN, b = {fmt(b_mm, 1)} mm, d = {fmt(d_mm, 1)} mm, "
            f"f'c = {fmt(fc_mpa, 2)} MPa, fyt specified = {fmt(fyt_mpa, 1)} MPa, "
            f"fyt used = {fmt(fyt_use, 1)} MPa"
        ),
        (
            f"φ = {fmt(phi, 2)}; stirrup {legs} legs of DB{fmt(stirrup_diameter_mm, 0)} "
            f"(Av = {fmt(av, 1)} mm²)"
        ),
    )
    report.add(
        "Concrete shear strength",
        "Vc from the code plugin (Nu = 0). See assumptions for the coefficient.",
        f"b = {fmt(b_mm, 1)} mm, d = {fmt(d_mm, 1)} mm, f'c = {fmt(fc_mpa, 2)} MPa",
        f"Vc = {fmt(n_to_kn(vc), 2)} kN; φVc = {fmt(n_to_kn(phi_vc), 2)} kN",
    )

    needs_strength = vu_n > phi_vc
    needs_min = vu_n > trigger
    if not needs_strength and not needs_min:
        report.add(
            "Stirrups are not required",
            "Vu ≤ φVc and Vu is not above the Av,min trigger",
            f"Vu = {fmt(vu_kN, 2)} kN, φVc = {fmt(n_to_kn(phi_vc), 2)} kN, trigger = {fmt(n_to_kn(trigger), 2)} kN",
            "No shear reinforcement required for strength or minimum steel.",
            note="A practical beam still needs hanger bars and support-region stirrups. v0.1 does not add those.",
        )
        report.summary = [
            "No stirrups required by the strength or Av,min rules at this station.",
            f"φVc = {fmt(n_to_kn(phi_vc), 2)} kN ≥ Vu = {fmt(vu_kN, 2)} kN",
        ]
        return ShearResult(report.status, report, False, av, 0.0, n_to_kn(phi_vc), vu_n / phi_vc if phi_vc else 0.0)

    vs_req = max(0.0, vu_n / phi - vc) if needs_strength else 0.0
    vs_max = code.vs_max_n(fc_mpa, b_mm, d_mm, av_meets_min=True)
    report.add(
        "Required steel shear strength",
        "Vs,req = Vu/φ − Vc, but not less than 0 when only Av,min is required",
        f"Vu/φ = {fmt(n_to_kn(vu_n / phi), 2)} kN; Vc = {fmt(n_to_kn(vc), 2)} kN",
        f"Vs,req = {fmt(n_to_kn(vs_req), 2)} kN; Vs,max = {fmt(n_to_kn(vs_max), 2)} kN",
    )
    if vs_req > vs_max * (1.0 + 1.0e-9):
        report.status = "NG"
        report.summary = [
            "NG: the section is too small. Vs,req exceeds Vs,max. Increase b, d, or f'c.",
            f"Vs,req = {fmt(n_to_kn(vs_req), 2)} kN > Vs,max = {fmt(n_to_kn(vs_max), 2)} kN",
        ]
        report.add(
            "Section is inadequate",
            "Vs,req ≤ Vs,max",
            f"{fmt(n_to_kn(vs_req), 2)} kN ≤ {fmt(n_to_kn(vs_max), 2)} kN",
            "Not satisfied.",
        )
        return ShearResult("NG", report, True, av, 0.0, n_to_kn(phi * (vc + vs_max)), vu_n / (phi * (vc + vs_max)))

    s_strength = av * fyt_use * d_mm / vs_req if vs_req > 0 else float("inf")
    s_min_steel = av / av_per_s if av_per_s > 0 else float("inf")
    s_max = code.s_max_mm(d_mm, vs_req, fc_mpa, b_mm)
    s_gov = min(s_strength, s_min_steel, s_max)
    s_prov = _round_down(s_gov, spacing_step_mm)
    report.add(
        "Stirrup spacing",
        "s = min( Av fyt d / Vs,req , Av / (Av,min/s) , s_max ), then rounded down",
        (
            f"s by strength = {fmt(s_strength, 1)} mm; s by Av,min = {fmt(s_min_steel, 1)} mm; "
            f"s_max = {fmt(s_max, 1)} mm; Av,min/s = {fmt(av_per_s, 4)} mm²/mm"
        ),
        f"s governing before rounding = {fmt(s_gov, 1)} mm; s provided = {fmt(s_prov, 1)} mm",
    )
    if s_prov <= 0:
        report.status = "NG"
        report.summary = ["NG: required spacing rounds to zero. Use a larger stirrup bar."]
        return ShearResult("NG", report, True, av, 0.0, 0.0, float("inf"))

    vs_prov = av * fyt_use * d_mm / s_prov
    phi_vn = phi * (vc + vs_prov)
    utilization = vu_n / phi_vn
    status = "OK" if utilization <= 1.001 else "NG"
    note = ""
    if s_prov < 50:
        note = "Spacing below 50 mm is tight to place. Confirm detailing; this is not a code failure by itself."
        if status == "OK":
            status = "REVIEW"
    report.status = status
    report.add(
        "Strength with the provided spacing",
        "Vs = Av fyt d / s;  φVn = φ (Vc + Vs);  utilization = Vu / φVn",
        f"s = {fmt(s_prov, 1)} mm, Av = {fmt(av, 1)} mm², fyt = {fmt(fyt_use, 1)} MPa, d = {fmt(d_mm, 1)} mm",
        (
            f"Vs = {fmt(n_to_kn(vs_prov), 2)} kN; φVn = {fmt(n_to_kn(phi_vn), 2)} kN; "
            f"Vu/φVn = {fmt(utilization, 3)}"
        ),
        note=note,
    )
    report.summary = [
        f"{legs}-leg DB{fmt(stirrup_diameter_mm, 0)} at {fmt(s_prov, 0)} mm",
        f"φVn = {fmt(n_to_kn(phi_vn), 2)} kN, Vu = {fmt(vu_kN, 2)} kN, utilization = {fmt(utilization, 3)}",
        "Station shear Vu is an input. v0.1 does not decide the critical section by itself.",
    ]
    return ShearResult(status, report, True, av, s_prov, n_to_kn(phi_vn), utilization)


def _round_down(value_mm: float, step_mm: float) -> float:
    if math.isinf(value_mm):
        return value_mm
    if step_mm <= 0:
        return value_mm
    steps = math.floor(value_mm / step_mm + 1.0e-9)
    return steps * step_mm
