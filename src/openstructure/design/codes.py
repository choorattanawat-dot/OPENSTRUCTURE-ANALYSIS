"""Design-code plugins.

ACI 318-19 factors below are the SI expressions of that code. EIT factors that
are confirmed by published Thai sources are marked as such. Factors that follow
the ACI 318-89 generation, which วสท. 1008 is based on, but that were not read
from a printed clause of วสท. 1008, are marked ASSUMPTION. See
``docs/CODE_ASSUMPTIONS.md``.
"""

from __future__ import annotations

import math
from dataclasses import dataclass

from openstructure.errors import OpenStructureError
from openstructure.units import kgf_to_n, mpa_to_ksc

PSI_PER_MPA = 1.0 / 0.006894757
ES_MPA = 200_000.0
EPS_CU = 0.003


def phi_transition_from_strain(
    eps_t: float,
    fy_mpa: float,
    phi_compression: float,
    phi_tension: float = 0.90,
    es_mpa: float = ES_MPA,
) -> float:
    """ACI 318-19 §21.2.2. εt ≥ 0.005 is tension-controlled."""
    eps_y = fy_mpa / es_mpa
    if eps_t <= eps_y:
        return phi_compression
    if eps_t >= 0.005:
        return phi_tension
    return phi_compression + (phi_tension - phi_compression) * (eps_t - eps_y) / (0.005 - eps_y)


@dataclass(frozen=True)
class CodeFactors:
    """Resolved numeric factors for one code and one quality-control case."""

    code_id: str
    code_name: str
    phi_flexure: float
    phi_shear: float
    phi_tied: float
    phi_spiral: float
    phi_bearing: float
    assumptions: tuple[str, ...]


class ACI318:
    """ACI 318-19 strength design, SI form, normalweight concrete (λ = 1)."""

    code_id = "ACI318-19"

    def __init__(self, lambda_lightweight: float = 1.0) -> None:
        if not 0.75 <= lambda_lightweight <= 1.0:
            raise OpenStructureError("ACI lightweight factor λ must be between 0.75 and 1.0.")
        self.lambda_lightweight = lambda_lightweight

    @property
    def code_name(self) -> str:
        return "ACI 318-19"

    def factors(self) -> CodeFactors:
        return CodeFactors(
            code_id=self.code_id,
            code_name=self.code_name,
            phi_flexure=0.90,
            phi_shear=0.75,
            phi_tied=0.65,
            phi_spiral=0.75,
            phi_bearing=0.65,
            assumptions=self.assumptions(),
        )

    def assumptions(self) -> tuple[str, ...]:
        return (
            "ACI 318-19 §21.2.1: φ = 0.90 tension-controlled flexure, 0.75 shear, "
            "0.65 tied compression-controlled, 0.75 spiral compression-controlled.",
            "ACI 318-19 §21.2.2: φ transitions linearly with εt from εy to 0.005.",
            "ACI 318-19 §22.2.2.4: Whitney block 0.85 f'c, a = β1 c, εcu = 0.003, "
            "β1 = 0.85 for f'c ≤ 28 MPa, then −0.05 per 7 MPa, not below 0.65.",
            "ACI 318-19 §9.6.1.2: As,min = max(0.25√f'c / fy, 1.4 / fy) bw d.",
            "Beam maximum reinforcement is the tension-controlled limit εt = 0.005 "
            "(§9.3.3.1 also requires εt ≥ 0.004 for a flexural member).",
            "ACI 318-19 §22.5.5.1 equation (a): Vc = 0.17 λ √f'c bw d, used for beams "
            "with Nu = 0. Equation (b), which may give a higher Vc, is not used in v0.1.",
            "√f'c used in Vc is capped at 8.3 MPa unless Av ≥ Av,min "
            "(ACI 318-19 §22.5.3.1; confirm the printed exception).",
            "ACI 318-19 §22.5.1.2: Vs ≤ 0.66 √f'c bw d.",
            "ACI 318-19 §9.6.3: Av,min/s = max(0.062√f'c, 0.35) bw / fyt. "
            "Av,min is required where Vu > 0.083 φ λ √f'c bw d, the SI form of the "
            "psi trigger in §9.6.3.1. Confirm the printed SI coefficient 0.083.",
            "Stirrup spacing limits follow §9.7.6.2.2: min(d/2, 600 mm) when "
            "Vs ≤ 0.33√f'c bw d, otherwise min(d/4, 300 mm). 0.33 is the SI form of 4√f'c (psi).",
            "fyt used for Vs and Av,min is capped at 420 MPa (ACI 318-19 §20.2.2.4).",
            "ACI 318-19 §22.4.2: Po = 0.85 f'c (Ag − Ast) + fy Ast, "
            "Pn,max = 0.80 Po tied or 0.85 Po spiral.",
            "ACI 318-19 §10.6.1.1: 0.01 Ag ≤ Ast ≤ 0.08 Ag.",
            "Es = 200000 MPa (§20.2.2.2). λ = 1 unless a lightweight factor is set. "
            "The §5.3.3 live-load reduction to 0.5 is not applied automatically.",
            "v0.1 does not check slenderness, shear in columns, torsion, development "
            "length, crack control, or deflection.",
        )

    def beta1(self, fc_mpa: float) -> float:
        if fc_mpa <= 28.0:
            return 0.85
        return max(0.65, 0.85 - 0.05 * (fc_mpa - 28.0) / 7.0)

    def c_over_d_limit(self, fy_mpa: float) -> float:
        """Neutral-axis limit for a tension-controlled beam, εt = 0.005."""
        del fy_mpa
        return EPS_CU / (EPS_CU + 0.005)

    def phi_flexure_for_strain(self, eps_t: float, fy_mpa: float) -> float:
        return phi_transition_from_strain(eps_t, fy_mpa, phi_compression=0.65, phi_tension=0.90)

    def phi_column(self, eps_t: float, fy_mpa: float, tied: bool, pn_n: float, **_ignored: float) -> float:
        del pn_n
        phi_c = 0.65 if tied else 0.75
        return phi_transition_from_strain(eps_t, fy_mpa, phi_compression=phi_c, phi_tension=0.90)

    def as_min_beam_mm2(self, fc_mpa: float, fy_mpa: float, bw_mm: float, d_mm: float) -> float:
        rho = max(0.25 * math.sqrt(fc_mpa) / fy_mpa, 1.4 / fy_mpa)
        return rho * bw_mm * d_mm

    def vc_n(
        self,
        fc_mpa: float,
        bw_mm: float,
        d_mm: float,
        av_meets_min: bool,
    ) -> float:
        root = math.sqrt(fc_mpa)
        if not av_meets_min:
            root = min(root, 8.3)
        return 0.17 * self.lambda_lightweight * root * bw_mm * d_mm

    def vs_max_n(self, fc_mpa: float, bw_mm: float, d_mm: float, av_meets_min: bool) -> float:
        root = math.sqrt(fc_mpa)
        if not av_meets_min:
            root = min(root, 8.3)
        return 0.66 * root * bw_mm * d_mm

    def av_min_per_mm(self, fc_mpa: float, bw_mm: float, fyt_mpa: float) -> float:
        fyt = min(fyt_mpa, 420.0)
        return max(0.062 * math.sqrt(fc_mpa) * bw_mm / fyt, 0.35 * bw_mm / fyt)

    def av_min_trigger_n(self, fc_mpa: float, bw_mm: float, d_mm: float) -> float:
        root = min(math.sqrt(fc_mpa), 8.3)
        return 0.75 * self.lambda_lightweight * 0.083 * root * bw_mm * d_mm

    def s_max_mm(self, d_mm: float, vs_n: float, fc_mpa: float, bw_mm: float) -> float:
        threshold = 0.33 * math.sqrt(fc_mpa) * bw_mm * d_mm
        if vs_n <= threshold:
            return min(d_mm / 2.0, 600.0)
        return min(d_mm / 4.0, 300.0)

    def fyt_for_design(self, fyt_mpa: float) -> float:
        return min(fyt_mpa, 420.0)

    def pn_max_factor(self, tied: bool) -> float:
        return 0.80 if tied else 0.85

    def phi_compression(self, tied: bool) -> float:
        return 0.65 if tied else 0.75


class EIT1008:
    """วสท. 1008 strength design, with the confirmed factors and flagged assumptions.

    ``quality_control="good"`` is the normal case (φ flexure 0.90, tied 0.70,
    spiral 0.75, shear 0.85). ``quality_control="ordinary"`` uses the reduced
    set quoted from a secondary Thai summary of the standard (not from a clause
    read in the printed book): flexure 0.75, tied 0.60, spiral 0.625, shear 0.70.
    """

    code_id = "EIT1008"

    def __init__(self, quality_control: str = "good") -> None:
        if quality_control not in ("good", "ordinary"):
            raise OpenStructureError("EIT quality_control must be 'good' or 'ordinary'.")
        self.quality_control = quality_control
        if quality_control == "good":
            self._phi_flexure = 0.90
            self._phi_tied = 0.70
            self._phi_spiral = 0.75
            self._phi_shear = 0.85
            self._phi_bearing = 0.70
        else:
            self._phi_flexure = 0.75
            self._phi_tied = 0.60
            self._phi_spiral = 0.625
            self._phi_shear = 0.70
            self._phi_bearing = 0.60

    @property
    def code_name(self) -> str:
        label = "good QC" if self.quality_control == "good" else "ordinary QC (reduced φ)"
        return f"วสท. 1008 / EIT 1008 ({label})"

    def factors(self) -> CodeFactors:
        return CodeFactors(
            code_id=self.code_id,
            code_name=self.code_name,
            phi_flexure=self._phi_flexure,
            phi_shear=self._phi_shear,
            phi_tied=self._phi_tied,
            phi_spiral=self._phi_spiral,
            phi_bearing=self._phi_bearing,
            assumptions=self.assumptions(),
        )

    def assumptions(self) -> tuple[str, ...]:
        return (
            "CONFIRMED in secondary Thai sources (TCA journal on ว.ส.ท. 1008 vs ACI 318-05, "
            "and วสท. 1008 column teaching notes): φ flexure = 0.90, φ shear = 0.85, "
            "φ tied = 0.70, φ spiral = 0.75 when workmanship and materials are controlled; "
            "U = 1.4D + 1.7L; Po formula with Pn,max = 0.80 Po (tied) or 0.85 Po (spiral); "
            "0.01 Ag ≤ Ast ≤ 0.08 Ag.",
            "ASSUMPTION: ordinary-QC φ (flexure 0.75, tied 0.60, spiral 0.625, shear 0.70, "
            "bearing 0.60) is taken from a university summary of the two φ sets in the Thai "
            "standard, not from the printed clause. Default is good QC.",
            "ASSUMPTION: β1 = 0.85 for f'c ≤ 280 ksc, then −0.05 per 70 ksc, not below 0.65 "
            "(Thai textbook rounding of the ACI 318-89 rule). The parent rule uses 4000 psi.",
            "ASSUMPTION: beam ρmax = 0.75 ρb (ACI 318-89), not the later εt = 0.005 limit. "
            "φ for flexure stays at the flexure value when the tension steel yields.",
            "ASSUMPTION: column φ transitions from the compression value to the flexure value "
            "as φPn decreases from 0.10 f'c Ag to zero for symmetric sections with "
            "(distance between extreme bar layers)/h ≥ 0.70 and fy ≤ 420 MPa. Otherwise the "
            "start of the ramp is min(0.10 f'c Ag, φ Pb). This is the ACI 318-99 §9.3.2 "
            "rule. Confirm the printed วสท. 1008 transition.",
            "ASSUMPTION: As,min is the ACI 318-89 expression max(3√f'c / fy, 200 / fy) bw d "
            "evaluated in psi and converted. Textbooks often round this to 14/fy in ksc.",
            "ASSUMPTION: Vc = 0.53 √f'c b d in ksc, cm, and kgf (the usual rounded form of "
            "2√f'c in psi). Vs,max = 2.12 √f'c b d. Spacing breaks at 1.06 √f'c b d. "
            "Av,min/s = max(0.20 √f'c, 3.5) b / fy. Av,min is required when Vu > 0.5 φ Vc. "
            "These are the ksc textbook forms. Confirm the printed coefficients.",
            "EIT stirrup design in v0.1 does not cap fyt at 420 MPa. Confirm whether "
            "วสท. 1008 limits fyt.",
            "Es = 200000 MPa and εcu = 0.003, shared with ACI. Wind and earthquake "
            "combinations in the EIT preset are assumptions; see that module.",
            "v0.1 does not check slenderness, column shear, torsion, development, "
            "crack control, or deflection. The section check is for a short member.",
        )

    def beta1(self, fc_mpa: float) -> float:
        fc_ksc = mpa_to_ksc(fc_mpa)
        if fc_ksc <= 280.0:
            return 0.85
        return max(0.65, 0.85 - 0.05 * (fc_ksc - 280.0) / 70.0)

    def c_over_d_limit(self, fy_mpa: float) -> float:
        eps_y = fy_mpa / ES_MPA
        c_balanced = EPS_CU / (EPS_CU + eps_y)
        return 0.75 * c_balanced

    def phi_flexure_for_strain(self, eps_t: float, fy_mpa: float) -> float:
        """Flexure φ if the tension steel yields; compression-controlled beams are rejected."""
        del eps_t, fy_mpa
        return self._phi_flexure

    def phi_column(
        self,
        eps_t: float,
        fy_mpa: float,
        tied: bool,
        pn_n: float,
        fc_mpa: float = 0.0,
        ag_mm2: float = 0.0,
        pb_n: float = 0.0,
        symmetric: bool = True,
        layer_ratio: float = 1.0,
    ) -> float:
        del eps_t
        phi_c = self._phi_tied if tied else self._phi_spiral
        phi_t = self._phi_flexure
        if pn_n <= 0.0:
            return phi_t
        limit = 0.10 * fc_mpa * ag_mm2
        special = symmetric and fy_mpa <= 420.0 and layer_ratio >= 0.70 - 1.0e-9
        if not special:
            limit = min(limit, phi_c * pb_n)
        if limit <= 0.0:
            return phi_c
        if pn_n >= limit / phi_c:
            return phi_c
        denominator = 1.0 - (phi_c - phi_t) * pn_n / limit
        if denominator <= 1.0e-9:
            return phi_c
        return phi_t / denominator

    def as_min_beam_mm2(self, fc_mpa: float, fy_mpa: float, bw_mm: float, d_mm: float) -> float:
        fc_psi = fc_mpa * PSI_PER_MPA
        fy_psi = fy_mpa * PSI_PER_MPA
        rho = max(3.0 * math.sqrt(fc_psi) / fy_psi, 200.0 / fy_psi)
        return rho * bw_mm * d_mm

    def vc_n(self, fc_mpa: float, bw_mm: float, d_mm: float, av_meets_min: bool) -> float:
        del av_meets_min
        fc_ksc = mpa_to_ksc(fc_mpa)
        vc_kgf = 0.53 * math.sqrt(fc_ksc) * (bw_mm / 10.0) * (d_mm / 10.0)
        return kgf_to_n(vc_kgf)

    def vs_max_n(self, fc_mpa: float, bw_mm: float, d_mm: float, av_meets_min: bool) -> float:
        del av_meets_min
        fc_ksc = mpa_to_ksc(fc_mpa)
        vs_kgf = 2.12 * math.sqrt(fc_ksc) * (bw_mm / 10.0) * (d_mm / 10.0)
        return kgf_to_n(vs_kgf)

    def av_min_per_mm(self, fc_mpa: float, bw_mm: float, fyt_mpa: float) -> float:
        fc_ksc = mpa_to_ksc(fc_mpa)
        fy_ksc = mpa_to_ksc(fyt_mpa)
        b_cm = bw_mm / 10.0
        # cm² per cm, converted to mm² per mm (the ratio is unchanged: *100/10 = *10, and
        # mm²/mm = 10 * cm²/cm). 1 cm²/cm = 10 mm²/mm.
        per_cm = max(0.20 * math.sqrt(fc_ksc), 3.5) * b_cm / fy_ksc
        return per_cm * 10.0

    def av_min_trigger_n(self, fc_mpa: float, bw_mm: float, d_mm: float) -> float:
        return 0.5 * self._phi_shear * self.vc_n(fc_mpa, bw_mm, d_mm, True)

    def s_max_mm(self, d_mm: float, vs_n: float, fc_mpa: float, bw_mm: float) -> float:
        fc_ksc = mpa_to_ksc(fc_mpa)
        threshold_n = kgf_to_n(1.06 * math.sqrt(fc_ksc) * (bw_mm / 10.0) * (d_mm / 10.0))
        if vs_n <= threshold_n:
            return min(d_mm / 2.0, 600.0)
        return min(d_mm / 4.0, 300.0)

    def fyt_for_design(self, fyt_mpa: float) -> float:
        return fyt_mpa

    def pn_max_factor(self, tied: bool) -> float:
        return 0.80 if tied else 0.85

    def phi_compression(self, tied: bool) -> float:
        return self._phi_tied if tied else self._phi_spiral


DesignCode = ACI318 | EIT1008


def get_code(name: str, **kwargs: object) -> DesignCode:
    key = name.strip().lower().replace("_", "").replace("-", "").replace(" ", "")
    if key in ("aci", "aci318", "aci31819", "aci318-19"):
        lambda_value = float(kwargs.get("lambda_lightweight", 1.0))  # type: ignore[arg-type]
        return ACI318(lambda_lightweight=lambda_value)
    if key in ("eit", "eit1008", "วสท", "วสท1008"):
        quality = str(kwargs.get("quality_control", "good"))
        return EIT1008(quality_control=quality)
    raise OpenStructureError(f"Unknown design code {name!r}. Use 'ACI318-19' or 'EIT1008'.")
