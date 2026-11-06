"""Nominal axial-flexure strength of a rectangular concrete section.

The Whitney rectangular stress block is used (0.85 f'c over a depth a = β1 c,
not past the far face). Bar stress is Es·ε, limited to ±fy. A bar inside the
stress block contributes (fs − 0.85 f'c) As so the concrete stress is not
counted twice.

Sign convention (compression positive for axial force):

* y is measured from the geometric centroid, positive toward one face.
* Positive moment compresses the positive-y face.
* εt is the strain in the bar farthest from the chosen compression face,
  positive when that bar is in tension.

Moments are about the geometric centroid. For symmetric reinforcement that is
also the plastic centroid. v0.1 does not shift the origin to the plastic
centroid of an unsymmetric layout; the report says so when the bars are not
symmetric.

This routine is the section kernel shared by the column interaction diagram
and by the beam capacity check. Beam proportioning itself is algebraic and is
written out in the calculation report.
"""

from __future__ import annotations

from dataclasses import dataclass

from openstructure.errors import OpenStructureError

Bar = tuple[float, float]  # (area_mm2, y_mm from centroid)


@dataclass(frozen=True)
class NominalResult:
    c_mm: float
    pn_n: float
    mn_nmm: float
    eps_t: float
    compression_sign: int
    a_mm: float
    forces_n: tuple[float, ...]


def pure_tension(
    bars: list[Bar],
    fy_mpa: float,
) -> NominalResult:
    axial = 0.0
    moment = 0.0
    for area, y_mm in bars:
        force = -fy_mpa * area
        axial += force
        moment += force * y_mm
    return NominalResult(
        c_mm=0.0,
        pn_n=axial,
        mn_nmm=moment,
        eps_t=1.0,
        compression_sign=0,
        a_mm=0.0,
        forces_n=tuple(-fy_mpa * area for area, _y in bars),
    )


def pure_compression(
    b_mm: float,
    h_mm: float,
    bars: list[Bar],
    fc_mpa: float,
    fy_mpa: float,
) -> NominalResult:
    """Concentric nominal strength Po, moment about the geometric centroid.

    Po = 0.85 f'c (Ag − Ast) + fy Ast. Valid as a strain-diagram point when
    fy ≤ 0.003 Es (600 MPa), because εcu = 0.003 then yields every bar.
    """
    ast = sum(area for area, _y in bars)
    ag = b_mm * h_mm
    po = 0.85 * fc_mpa * (ag - ast) + fy_mpa * ast
    moment = 0.0
    forces: list[float] = []
    for area, y_mm in bars:
        force = (fy_mpa - 0.85 * fc_mpa) * area
        forces.append(force)
        moment += force * y_mm
    return NominalResult(
        c_mm=float("inf"),
        pn_n=po,
        mn_nmm=moment,
        eps_t=-0.003,
        compression_sign=0,
        a_mm=h_mm,
        forces_n=tuple(forces),
    )


def nominal_at_c(
    b_mm: float,
    h_mm: float,
    bars: list[Bar],
    c_mm: float,
    fc_mpa: float,
    fy_mpa: float,
    beta1: float,
    compression_sign: int,
    eps_cu: float = 0.003,
    es_mpa: float = 200_000.0,
) -> NominalResult:
    """Nominal Pn, Mn for a neutral-axis depth measured from one face.

    ``compression_sign`` is +1 to put the compression face on +y and −1 to put
    it on −y.
    """
    if compression_sign not in (1, -1):
        raise OpenStructureError("compression_sign must be +1 or -1.")
    if c_mm <= 0:
        raise OpenStructureError("Neutral-axis depth c must be positive.")
    if b_mm <= 0 or h_mm <= 0:
        raise OpenStructureError("Section dimensions must be positive.")

    a_mm = min(beta1 * c_mm, h_mm)
    concrete = 0.85 * fc_mpa * b_mm * a_mm
    forces: list[float] = []
    steel_moment_raw = 0.0
    d_t = 0.0
    for area, y_mm in bars:
        y_from_comp = h_mm / 2.0 - compression_sign * y_mm
        d_t = max(d_t, y_from_comp)
        strain = eps_cu * (c_mm - y_from_comp) / c_mm
        stress = max(-fy_mpa, min(fy_mpa, strain * es_mpa))
        if y_from_comp <= a_mm + 1.0e-9:
            stress -= 0.85 * fc_mpa
        force = stress * area
        forces.append(force)
        steel_moment_raw += force * (compression_sign * y_mm)

    axial = concrete + sum(forces)
    # Lever of the concrete block is positive toward the compression face.
    concrete_moment = concrete * (h_mm / 2.0 - a_mm / 2.0)
    moment_raw = concrete_moment + steel_moment_raw
    moment = moment_raw if compression_sign == 1 else -moment_raw
    eps_t = eps_cu * (d_t - c_mm) / c_mm if d_t > 0 else -eps_cu
    return NominalResult(
        c_mm=c_mm,
        pn_n=axial,
        mn_nmm=moment,
        eps_t=eps_t,
        compression_sign=compression_sign,
        a_mm=a_mm,
        forces_n=tuple(forces),
    )


def solve_neutral_axis(
    b_mm: float,
    h_mm: float,
    bars: list[Bar],
    fc_mpa: float,
    fy_mpa: float,
    beta1: float,
    compression_sign: int = 1,
    axial_target_n: float = 0.0,
    eps_cu: float = 0.003,
    es_mpa: float = 200_000.0,
) -> NominalResult:
    """Find c such that the nominal axial force equals ``axial_target_n``.

    Used for pure flexure (target 0). Bisection on c. The residual is Pn(c)
    minus the target: negative when tension wins (c too small).
    """

    def residual(c_mm: float) -> float:
        point = nominal_at_c(
            b_mm, h_mm, bars, c_mm, fc_mpa, fy_mpa, beta1, compression_sign, eps_cu, es_mpa
        )
        return point.pn_n - axial_target_n

    low = h_mm / 1.0e4
    high = h_mm * 5.0
    r_low = residual(low)
    r_high = residual(high)
    expand = 0
    while r_low > 0 and low > h_mm / 1.0e8 and expand < 20:
        low *= 0.5
        r_low = residual(low)
        expand += 1
    expand = 0
    while r_high < 0 and high < h_mm * 500 and expand < 20:
        high *= 2.0
        r_high = residual(high)
        expand += 1
    if r_low > 0 or r_high < 0:
        raise OpenStructureError(
            "Could not bracket the neutral axis. Check the reinforcement and the axial force."
        )
    for _ in range(80):
        mid = 0.5 * (low + high)
        if residual(mid) > 0:
            high = mid
        else:
            low = mid
    return nominal_at_c(
        b_mm,
        h_mm,
        bars,
        0.5 * (low + high),
        fc_mpa,
        fy_mpa,
        beta1,
        compression_sign,
        eps_cu,
        es_mpa,
    )


def bars_are_symmetric(bars: list[Bar], tolerance_mm: float = 1.0) -> bool:
    """True when every bar at +y has a matching area at −y. Bars on the axis match themselves."""
    positives = [(area, y) for area, y in bars if y > tolerance_mm]
    negatives = [(area, y) for area, y in bars if y < -tolerance_mm]
    if len(positives) != len(negatives):
        return False
    used = [False] * len(negatives)
    for area, y_mm in positives:
        found = False
        for index, (other_area, other_y) in enumerate(negatives):
            if used[index]:
                continue
            same_place = abs(other_y + y_mm) <= tolerance_mm
            same_area = abs(other_area - area) <= 1.0e-6 * max(area, 1.0)
            if same_place and same_area:
                used[index] = True
                found = True
                break
        if not found:
            return False
    return True


def layer_ratio(h_mm: float, bars: list[Bar]) -> float:
    """Distance between extreme bar layers divided by h."""
    if not bars or h_mm <= 0:
        return 0.0
    ys = [y for _area, y in bars]
    return (max(ys) - min(ys)) / h_mm


def c_grid(h_mm: float) -> list[float]:
    values: list[float] = []
    c_mm = h_mm / 400.0
    while c_mm < h_mm * 40.0:
        values.append(c_mm)
        c_mm *= 1.15
    return values
