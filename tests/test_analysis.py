"""Linear-static benchmarks.

Closed-form beam results are the classical Euler-Bernoulli solutions tabulated
in structural-analysis textbooks (for example Hibbeler, Structural Analysis,
beam-deflection tables):

* simply supported beam, uniform load: δ = 5 w L^4 / 384 EI, M = w L^2 / 8
* fixed beam, uniform load: δ = w L^4 / 384 EI, M_end = w L^2 / 12, M_mid = w L^2 / 24
* simply supported beam, midspan point load: δ = P L^3 / 48 EI, M = P L / 4

The portal frame is checked two ways. A direct-stiffness solution written in
this file (axial plus bending, no shear deformation) must match OpenSees,
because both implement the same element. The rigid-girder closed form
Δ = P h^3 / 24 EI and M = P h / 4 neglects axial deformation; it is compared
with a loose tolerance and the difference is stated below.
"""

import math

import pytest

from openstructure.analysis import AnalysisOptions, LinearStaticEngine
from openstructure.model import Model

E_MPA = 25_000.0
B_MM = 300.0
H_MM = 600.0


def _inertia() -> float:
    return B_MM * H_MM**3 / 12.0


def _beam(length_m: float, support: str) -> Model:
    model = Model(title="beam")
    model.add_material("C", ec_mpa=E_MPA, fc_mpa=28, fy_mpa=400)
    model.add_rectangular_section("S", "C", B_MM, H_MM)
    model.add_node("A", 0.0, 0.0)
    model.add_node("B", length_m, 0.0)
    if support == "simple":
        model.add_support("A", pinned=True)
        model.add_support("B", roller_x=True)
    else:
        model.add_support("A", fixed=True)
        model.add_support("B", fixed=True)
    model.add_member("AB", "A", "B", "S")
    model.add_load_case("DL", "dead")
    return model


def test_simply_supported_udl_matches_beam_theory():
    length_m = 6.0
    w = 12.0  # kN/m downward
    model = _beam(length_m, "simple")
    model.add_member_udl_global("AB", "DL", wy_kN_per_m=-w)
    model.add_combination("DL", {"DL": 1.0})
    result = LinearStaticEngine().analyze(model).by_combination["DL"]
    ei = E_MPA * _inertia()  # N·mm²
    length_mm = length_m * 1000.0
    w_n_per_mm = w
    delta = 5.0 * w_n_per_mm * length_mm**4 / (384.0 * ei)
    moment = w * length_m**2 / 8.0
    reaction = w * length_m / 2.0
    assert result.displacements["A"].uy_mm == pytest.approx(0.0, abs=1e-8)
    # One node at midspan is not in the model. The member diagram carries M.
    mid = max(result.members["AB"].stations, key=lambda station: station.m_kNm)
    assert mid.m_kNm == pytest.approx(moment, rel=1e-6)
    assert mid.x_m == pytest.approx(length_m / 2.0, abs=1e-3)
    assert result.reactions["A"].fy_kN == pytest.approx(reaction, rel=1e-6)
    assert result.reactions["B"].fy_kN == pytest.approx(reaction, rel=1e-6)
    # Deflection of a simple beam is recovered by adding a mid node.
    mid_model = _beam_with_mid_node(length_m, w)
    mid_result = LinearStaticEngine().analyze(mid_model).by_case["DL"]
    assert mid_result.displacements["M"].uy_mm == pytest.approx(-delta, rel=1e-6)
    fx, fy, mz = mid_result.equilibrium_kN
    assert fx == pytest.approx(0.0, abs=1e-6)
    assert fy == pytest.approx(0.0, abs=1e-6)
    assert mz == pytest.approx(0.0, abs=1e-4)


def test_fixed_beam_udl_matches_beam_theory():
    length_m = 4.0
    w = 10.0
    model = _beam(length_m, "fixed")
    model.add_member_udl_global("AB", "DL", wy_kN_per_m=-w)
    model.add_combination("DL", {"DL": 1.0})
    # Both ends fixed: zero free equations. The engine must still return fixed-end forces.
    result = LinearStaticEngine().analyze(model).by_combination["DL"]
    moment_end = w * length_m**2 / 12.0
    moment_mid = w * length_m**2 / 24.0
    reaction = w * length_m / 2.0
    member = result.members["AB"]
    assert member.stations[0].m_kNm == pytest.approx(-moment_end, rel=1e-6)
    assert member.stations[-1].m_kNm == pytest.approx(-moment_end, rel=1e-6)
    mid = max(member.stations, key=lambda station: station.m_kNm)
    assert mid.m_kNm == pytest.approx(moment_mid, rel=1e-6)
    assert result.reactions["A"].fy_kN == pytest.approx(reaction, rel=1e-6)
    assert result.reactions["A"].mz_kNm == pytest.approx(moment_end, rel=1e-6)
    assert result.reactions["B"].mz_kNm == pytest.approx(-moment_end, rel=1e-6)
    assert result.displacements["A"].uy_mm == pytest.approx(0.0, abs=1e-9)


def test_simply_supported_point_load_matches_beam_theory():
    length_m = 8.0
    load_kN = 24.0
    model = _beam(length_m, "simple")
    model.add_member_point_load_global("AB", "DL", 0.5, fy_kN=-load_kN)
    model.add_combination("DL", {"DL": 1.0})
    checked = _beam_with_mid_node(length_m, 0.0, point_kN=load_kN)
    # The single-member model checks the moment; the two-element model checks deflection.
    one = LinearStaticEngine().analyze(model).by_combination["DL"]
    two = LinearStaticEngine().analyze(checked).by_case["DL"]
    ei = E_MPA * _inertia()
    length_mm = length_m * 1000.0
    delta = (load_kN * 1000.0) * length_mm**3 / (48.0 * ei)
    assert max(s.m_kNm for s in one.members["AB"].stations) == pytest.approx(load_kN * length_m / 4.0, rel=1e-6)
    assert two.displacements["M"].uy_mm == pytest.approx(-delta, rel=1e-6)
    assert one.reactions["A"].fy_kN == pytest.approx(load_kN / 2.0, rel=1e-6)


def test_load_combination_scales_a_linear_case():
    model = _beam(5.0, "simple")
    model.add_member_udl_global("AB", "DL", wy_kN_per_m=-8.0)
    model.add_combination("1.4DL", {"DL": 1.4})
    result = LinearStaticEngine().analyze(model)
    base = max(s.m_kNm for s in result.by_case["DL"].members["AB"].stations)
    factored = max(s.m_kNm for s in result.by_combination["1.4DL"].members["AB"].stations)
    assert factored == pytest.approx(1.4 * base, rel=1e-9)


def test_portal_frame_matches_direct_stiffness_and_rigid_girder_formula():
    """Fixed portal, horizontal load at the left knee.

    Direct stiffness (this file) includes axial deformation. The rigid-girder
    formula Δ = P h³ / 24 EI, M = Ph/4 does not; the beam inertia is set very
    large and the comparison tolerance covers the axial shortening OpenSees keeps.
    """
    height_m = 3.0
    span_m = 6.0
    load_kN = 10.0
    model = _portal(height_m, span_m, beam_i_factor=1.0)
    model.add_nodal_load("B", "WX", fx_kN=load_kN)
    model.add_combination("WX", {"WX": 1.0})
    opens = LinearStaticEngine().analyze(model).by_combination["WX"]
    reference = _stiffness_portal(height_m, span_m, load_kN, beam_i_factor=1.0)
    assert opens.displacements["B"].ux_mm == pytest.approx(reference["ux_b"], rel=1e-5)
    assert opens.displacements["C"].ux_mm == pytest.approx(reference["ux_c"], rel=1e-5)
    assert opens.reactions["A"].mz_kNm == pytest.approx(reference["m_a"], rel=1e-5)
    assert opens.reactions["D"].mz_kNm == pytest.approx(reference["m_d"], rel=1e-5)

    stiff = _portal(height_m, span_m, beam_i_factor=1.0e5)
    stiff.add_nodal_load("B", "WX", fx_kN=load_kN)
    stiff.add_combination("WX", {"WX": 1.0})
    rigid = LinearStaticEngine().analyze(stiff, AnalysisOptions(n_stations=5)).by_combination["WX"]
    i_col = 400.0 * 400.0**3 / 12.0
    ei = E_MPA * i_col
    height_mm = height_m * 1000.0
    delta_mm = (load_kN * 1000.0) * height_mm**3 / (24.0 * ei)
    moment = load_kN * height_m / 4.0
    # Axial deformation of the girder makes the two knees differ by a few percent
    # of the bending sway. The average sway and the column moments stay close to
    # the rigid-girder values.
    average_sway = 0.5 * (rigid.displacements["B"].ux_mm + rigid.displacements["C"].ux_mm)
    assert average_sway == pytest.approx(delta_mm, rel=0.05)
    assert rigid.reactions["A"].mz_kNm == pytest.approx(moment, rel=0.05)
    assert rigid.reactions["D"].mz_kNm == pytest.approx(moment, rel=0.05)


def _beam_with_mid_node(length_m: float, w_kN_per_m: float, point_kN: float = 0.0) -> Model:
    model = Model()
    model.add_material("C", ec_mpa=E_MPA, fc_mpa=28, fy_mpa=400)
    model.add_rectangular_section("S", "C", B_MM, H_MM)
    model.add_node("A", 0.0, 0.0)
    model.add_node("M", length_m / 2.0, 0.0)
    model.add_node("B", length_m, 0.0)
    model.add_support("A", pinned=True)
    model.add_support("B", roller_x=True)
    model.add_member("AM", "A", "M", "S")
    model.add_member("MB", "M", "B", "S")
    model.add_load_case("DL", "dead")
    if w_kN_per_m:
        model.add_member_udl_global("AM", "DL", wy_kN_per_m=-w_kN_per_m)
        model.add_member_udl_global("MB", "DL", wy_kN_per_m=-w_kN_per_m)
    if point_kN:
        model.add_nodal_load("M", "DL", fy_kN=-point_kN)
    return model


def _portal(height_m: float, span_m: float, beam_i_factor: float) -> Model:
    model = Model(title="portal")
    model.add_material("C", ec_mpa=E_MPA, fc_mpa=28, fy_mpa=400)
    model.add_rectangular_section("COL", "C", 400, 400)
    model.add_rectangular_section("BEAM", "C", 300, 600, i_factor=beam_i_factor)
    model.add_node("A", 0.0, 0.0)
    model.add_node("B", 0.0, height_m)
    model.add_node("C", span_m, height_m)
    model.add_node("D", span_m, 0.0)
    model.add_support("A", fixed=True)
    model.add_support("D", fixed=True)
    model.add_member("AB", "A", "B", "COL")
    model.add_member("BC", "B", "C", "BEAM")
    model.add_member("DC", "D", "C", "COL")
    model.add_load_case("WX", "wind", reversible=False)
    return model


def _stiffness_portal(height_m: float, span_m: float, load_kN: float, beam_i_factor: float) -> dict[str, float]:
    """Independent 2D beam stiffness solution in N and mm. Returns mm and kN·m."""
    nodes = {
        "A": (0.0, 0.0),
        "B": (0.0, height_m * 1000.0),
        "C": (span_m * 1000.0, height_m * 1000.0),
        "D": (span_m * 1000.0, 0.0),
    }
    members = [
        ("A", "B", 400.0 * 400.0, 400.0**4 / 12.0),
        ("B", "C", 300.0 * 600.0, beam_i_factor * 300.0 * 600.0**3 / 12.0),
        ("D", "C", 400.0 * 400.0, 400.0**4 / 12.0),
    ]
    order = ["A", "B", "C", "D"]
    index = {name: 3 * i for i, name in enumerate(order)}
    size = 12
    stiffness = [[0.0] * size for _ in range(size)]
    load = [0.0] * size
    load[index["B"]] = load_kN * 1000.0
    for i_name, j_name, area, inertia in members:
        _add_beam(stiffness, nodes[i_name], nodes[j_name], index[i_name], index[j_name], area, inertia)
    fixed = list(range(3)) + list(range(index["D"], index["D"] + 3))
    free = [dof for dof in range(size) if dof not in fixed]
    disp = [0.0] * size
    solved = _solve([row[:] for row in stiffness], load[:], free)
    for dof, value in zip(free, solved, strict=True):
        disp[dof] = value
    reactions = _multiply(stiffness, disp)
    return {
        "ux_b": disp[index["B"]],
        "ux_c": disp[index["C"]],
        "m_a": (reactions[index["A"] + 2] - load[index["A"] + 2]) / 1.0e6,
        "m_d": (reactions[index["D"] + 2] - load[index["D"] + 2]) / 1.0e6,
    }


def _add_beam(stiffness, i_xy, j_xy, i_dof, j_dof, area, inertia) -> None:
    dx = j_xy[0] - i_xy[0]
    dy = j_xy[1] - i_xy[1]
    length = math.hypot(dx, dy)
    cosine, sine = dx / length, dy / length
    ea_l = E_MPA * area / length
    ei_l = E_MPA * inertia / length
    ei_l2 = ei_l / length
    ei_l3 = ei_l2 / length
    local = [
        [ea_l, 0.0, 0.0, -ea_l, 0.0, 0.0],
        [0.0, 12 * ei_l3, 6 * ei_l2, 0.0, -12 * ei_l3, 6 * ei_l2],
        [0.0, 6 * ei_l2, 4 * ei_l, 0.0, -6 * ei_l2, 2 * ei_l],
        [-ea_l, 0.0, 0.0, ea_l, 0.0, 0.0],
        [0.0, -12 * ei_l3, -6 * ei_l2, 0.0, 12 * ei_l3, -6 * ei_l2],
        [0.0, 6 * ei_l2, 2 * ei_l, 0.0, -6 * ei_l2, 4 * ei_l],
    ]
    rotation = [
        [cosine, sine, 0.0],
        [-sine, cosine, 0.0],
        [0.0, 0.0, 1.0],
    ]
    transform = [[0.0] * 6 for _ in range(6)]
    for block_row in (0, 1):
        for block_col in (0, 1):
            for r in range(3):
                for c in range(3):
                    transform[3 * block_row + r][3 * block_col + c] = rotation[r][c] if block_row == block_col else 0.0
    # k_global = T^T k T, with T rotating global DOF into local DOF.
    kt = _matmul(local, transform)
    global_k = _matmul(_transpose(transform), kt)
    dofs = [i_dof, i_dof + 1, i_dof + 2, j_dof, j_dof + 1, j_dof + 2]
    for row, dof_row in enumerate(dofs):
        for col, dof_col in enumerate(dofs):
            stiffness[dof_row][dof_col] += global_k[row][col]


def _solve(matrix, rhs, free):
    reduced = [[matrix[i][j] for j in free] for i in free]
    vector = [rhs[i] for i in free]
    return _gauss(reduced, vector)


def _gauss(matrix, vector):
    n = len(vector)
    for col in range(n):
        pivot = max(range(col, n), key=lambda row: abs(matrix[row][col]))
        matrix[col], matrix[pivot] = matrix[pivot], matrix[col]
        vector[col], vector[pivot] = vector[pivot], vector[col]
        scale = matrix[col][col]
        if abs(scale) < 1e-12:
            raise AssertionError("Reference frame is singular.")
        for row in range(col + 1, n):
            factor = matrix[row][col] / scale
            for k in range(col, n):
                matrix[row][k] -= factor * matrix[col][k]
            vector[row] -= factor * vector[col]
    answer = [0.0] * n
    for row in range(n - 1, -1, -1):
        total = vector[row] - sum(matrix[row][k] * answer[k] for k in range(row + 1, n))
        answer[row] = total / matrix[row][row]
    return answer


def _matmul(a, b):
    rows, cols, inner = len(a), len(b[0]), len(b)
    out = [[0.0] * cols for _ in range(rows)]
    for i in range(rows):
        for k in range(inner):
            for j in range(cols):
                out[i][j] += a[i][k] * b[k][j]
    return out


def _transpose(a):
    return [list(row) for row in zip(*a, strict=True)]


def _multiply(matrix, vector):
    return [sum(value * vector[col] for col, value in enumerate(row)) for row in matrix]
