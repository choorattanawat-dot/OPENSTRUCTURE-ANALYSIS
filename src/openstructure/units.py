"""Unit conversions.

Analysis inside OpenSees uses a consistent set: N, mm, MPa (N/mm²), N·mm.
One kN/m equals one N/mm, which keeps distributed loads numerically the same
when converting from kN/m to N/mm.

EIT textbook formulas written in ksc (kgf/cm²) and centimetres are converted
with standard gravity: 1 kgf/cm² = 0.0980665 MPa.
"""

from __future__ import annotations

MPA_PER_KSC = 0.0980665
KSC_PER_MPA = 1.0 / MPA_PER_KSC
CM_PER_MM = 0.1
N_PER_KGF = 9.80665

# Public engineering units <-> internal N, mm.
N_PER_KN = 1000.0
MM_PER_M = 1000.0
NMM_PER_KNM = 1.0e6  # 1 kN·m = 10^6 N·mm


def kn_to_n(value_kN: float) -> float:
    return value_kN * N_PER_KN


def n_to_kn(value_n: float) -> float:
    return value_n / N_PER_KN


def m_to_mm(value_m: float) -> float:
    return value_m * MM_PER_M


def mm_to_m(value_mm: float) -> float:
    return value_mm / MM_PER_M


def knm_to_nmm(value_kNm: float) -> float:
    return value_kNm * NMM_PER_KNM


def nmm_to_knm(value_nmm: float) -> float:
    return value_nmm / NMM_PER_KNM


def kn_per_m_to_n_per_mm(value: float) -> float:
    """kN/m and N/mm are the same number."""
    return value


def mpa_to_ksc(fc_mpa: float) -> float:
    return fc_mpa * KSC_PER_MPA


def ksc_to_mpa(fc_ksc: float) -> float:
    return fc_ksc * MPA_PER_KSC


def kgf_to_n(value_kgf: float) -> float:
    return value_kgf * N_PER_KGF


def n_to_kgf(value_n: float) -> float:
    return value_n / N_PER_KGF
