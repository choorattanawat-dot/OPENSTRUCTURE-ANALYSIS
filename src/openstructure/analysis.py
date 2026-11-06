"""Linear-static analysis of a 2D frame with OpenSeesPy.

Only ``kind="linear_static"`` is implemented. Load combinations are applied by
superposition, which is valid while the analysis is linear. A nonlinear engine
should keep the same model and combination objects and apply each combination
as one load pattern.

Sign convention
---------------
Global axes: x to the right, y upward, moment counter-clockwise.

OpenSees ``localForce`` is stored as the force exerted by the joint on the
member, in the member local axes (local x from node i to node j, local y 90°
counter-clockwise from x).

Section forces along the member, used for design and for the diagrams:

* axial force positive in tension
* shear positive when the force exerted by the left side on the right side
  points in the local +y direction (for a left-to-right beam, upward shear at
  the left end of a downward-loaded simple beam is positive)
* bending moment positive when it compresses the local +y face (sagging is
  positive for a left-to-right beam)

The fixed-end moment of a downward uniform load is therefore negative.
"""

from __future__ import annotations

from dataclasses import dataclass, field

from openstructure.errors import OpenStructureError
from openstructure.model import (
    Combination,
    MemberPointLoad,
    MemberUDL,
    Model,
)
from openstructure.units import n_to_kn, nmm_to_knm

_DUMMY_NODE = 9_000_001
_DUMMY_ELEMENT = 9_000_001


@dataclass
class AnalysisOptions:
    """Options for the analysis engine.

    ``geometric_transformation`` is passed to OpenSees. ``Linear`` is the
    verified v0.1 setting. ``PDelta`` is accepted so a second-order run can be
    requested without a new model class; it is not part of the v0.1 benchmarks.
    ``kind`` other than ``linear_static`` is rejected with a clear error.
    """

    kind: str = "linear_static"
    geometric_transformation: str = "Linear"
    n_stations: int = 21


@dataclass
class NodeDisplacement:
    ux_mm: float
    uy_mm: float
    rz_rad: float


@dataclass
class NodeReaction:
    """Support action on the node. Forces in kN, moment in kN·m, global axes."""

    fx_kN: float
    fy_kN: float
    mz_kNm: float


@dataclass
class EndActions:
    """Joint-on-member end forces in local axes. Forces in kN, moments in kN·m."""

    fx_i_kN: float
    fy_i_kN: float
    m_i_kNm: float
    fx_j_kN: float
    fy_j_kN: float
    m_j_kNm: float


@dataclass
class LocalPointLoad:
    x_over_l: float
    px_kN: float
    py_kN: float


@dataclass
class StationForce:
    """Internal section forces. See the module docstring for signs."""

    x_m: float
    n_kN: float
    v_kN: float
    m_kNm: float


@dataclass
class MemberForces:
    end: EndActions
    wx_kN_per_m: float
    wy_kN_per_m: float
    point_loads: list[LocalPointLoad]
    stations: list[StationForce]
    length_m: float

    def extreme_moment(self) -> tuple[StationForce, StationForce]:
        """Return the stations of minimum and maximum moment (sagging positive)."""
        return min(self.stations, key=lambda s: s.m_kNm), max(self.stations, key=lambda s: s.m_kNm)

    def extreme_shear(self) -> StationForce:
        return max(self.stations, key=lambda s: abs(s.v_kN))


@dataclass
class CaseResult:
    name: str
    displacements: dict[str, NodeDisplacement]
    reactions: dict[str, NodeReaction]
    members: dict[str, MemberForces]
    equilibrium_kN: tuple[float, float, float] = (0.0, 0.0, 0.0)
    reference: str = ""
    expression: str = ""


@dataclass
class AnalysisResult:
    by_case: dict[str, CaseResult] = field(default_factory=dict)
    by_combination: dict[str, CaseResult] = field(default_factory=dict)
    options: AnalysisOptions = field(default_factory=AnalysisOptions)


class LinearStaticEngine:
    """OpenSeesPy linear-static solver.

    The OpenSees model builder is a process-wide global. This engine calls
    ``wipe()`` at the start of every case and is not safe to run from two
    threads at once.
    """

    def analyze(
        self,
        model: Model,
        options: AnalysisOptions | None = None,
        combinations: list[str] | None = None,
    ) -> AnalysisResult:
        options = options or AnalysisOptions()
        if options.kind != "linear_static":
            raise OpenStructureError(
                f"Analysis kind {options.kind!r} is not implemented in v0.1. "
                "Nodes, members, load cases, and combinations are independent of "
                "the solver, so a nonlinear engine (pushover, fiber sections) can "
                "be added without changing the model."
            )
        if options.geometric_transformation not in ("Linear", "PDelta"):
            raise OpenStructureError(
                "geometric_transformation must be 'Linear' or 'PDelta' in v0.1."
            )
        if options.n_stations < 2:
            raise OpenStructureError("n_stations must be at least 2.")
        model.validate()
        if not model.load_cases:
            raise OpenStructureError("The model has no load cases.")

        result = AnalysisResult(options=options)
        for case_id in model.load_cases:
            result.by_case[case_id] = self._analyze_case(model, case_id, options)

        selected = combinations
        if selected is None:
            selected = list(model.combinations)
        for name in selected:
            if name not in model.combinations:
                raise OpenStructureError(f"Unknown combination {name!r}.")
            result.by_combination[name] = _superpose(
                model, result.by_case, model.combinations[name], options.n_stations
            )
        return result

    def _analyze_case(self, model: Model, case_id: str, options: AnalysisOptions) -> CaseResult:
        import openseespy.opensees as ops

        ops.wipe()
        ops.model("basic", "-ndm", 2, "-ndf", 3)
        node_tag = {node_id: index for index, node_id in enumerate(model.nodes, start=1)}
        member_tag = {member_id: index for index, member_id in enumerate(model.members, start=1)}

        for node in model.nodes.values():
            ops.node(node_tag[node.id], node.x_mm, node.y_mm)
        for node_id, support in model.supports.items():
            ops.fix(
                node_tag[node_id],
                int(support.ux),
                int(support.uy),
                int(support.rz),
            )

        ops.geomTransf(options.geometric_transformation, 1)
        for member in model.members.values():
            section = model.sections[member.section]
            material = model.materials[section.material]
            ops.element(
                "elasticBeamColumn",
                member_tag[member.id],
                node_tag[member.i_node],
                node_tag[member.j_node],
                section.area_mm2,
                material.ec_mpa,
                section.inertia_mm4,
                1,
            )

        _add_dummy_dof_if_needed(ops, model, node_tag)

        ops.timeSeries("Linear", 1)
        ops.pattern("Plain", 1, 1)
        _apply_nodal_loads(ops, model, case_id, node_tag)
        local_loads = {}
        for member_id in model.members:
            wx, wy, points = _member_local_loads(model, member_id, case_id)
            local_loads[member_id] = (wx, wy, points)
            tag = member_tag[member_id]
            if abs(wx) > 0 or abs(wy) > 0:
                ops.eleLoad("-ele", tag, "-type", "-beamUniform", wy, wx)
            for point in points:
                ops.eleLoad(
                    "-ele",
                    tag,
                    "-type",
                    "-beamPoint",
                    point.py_kN * 1000.0,
                    point.x_over_l,
                    point.px_kN * 1000.0,
                )

        ops.system("BandGeneral")
        ops.numberer("RCM")
        ops.constraints("Plain")
        ops.test("NormDispIncr", 1.0e-8, 10, 0)
        ops.algorithm("Linear")
        ops.integrator("LoadControl", 1.0)
        ops.analysis("Static")
        status = ops.analyze(1)
        if status != 0:
            ops.wipe()
            raise OpenStructureError(
                f"OpenSees did not solve load case {case_id!r} (status {status}). "
                "The frame may be a mechanism, or a support may be missing."
            )
        ops.reactions()

        displacements = {
            node_id: NodeDisplacement(
                ux_mm=float(ops.nodeDisp(node_tag[node_id], 1)),
                uy_mm=float(ops.nodeDisp(node_tag[node_id], 2)),
                rz_rad=float(ops.nodeDisp(node_tag[node_id], 3)),
            )
            for node_id in model.nodes
        }
        reactions = {}
        for node_id in model.supports:
            raw = ops.nodeReaction(node_tag[node_id])
            reactions[node_id] = NodeReaction(
                fx_kN=n_to_kn(float(raw[0])),
                fy_kN=n_to_kn(float(raw[1])),
                mz_kNm=nmm_to_knm(float(raw[2])),
            )
        members: dict[str, MemberForces] = {}
        for member_id, tag in member_tag.items():
            raw = [float(value) for value in ops.eleResponse(tag, "localForce")]
            end = EndActions(
                fx_i_kN=n_to_kn(raw[0]),
                fy_i_kN=n_to_kn(raw[1]),
                m_i_kNm=nmm_to_knm(raw[2]),
                fx_j_kN=n_to_kn(raw[3]),
                fy_j_kN=n_to_kn(raw[4]),
                m_j_kNm=nmm_to_knm(raw[5]),
            )
            wx, wy, points = local_loads[member_id]
            length_m = model.member_length_mm(model.members[member_id].i_node, model.members[member_id].j_node) / 1000.0
            members[member_id] = MemberForces(
                end=end,
                wx_kN_per_m=wx,
                wy_kN_per_m=wy,
                point_loads=points,
                stations=section_stations(length_m, end, wx, wy, points, options.n_stations),
                length_m=length_m,
            )

        case = CaseResult(
            name=case_id,
            displacements=displacements,
            reactions=reactions,
            members=members,
            expression=case_id,
            reference="Single load case, factor 1.0, linear static OpenSeesPy.",
        )
        case.equilibrium_kN = equilibrium_residual(model, case, {case_id: 1.0})
        ops.wipe()
        return case


def section_stations(
    length_m: float,
    end: EndActions,
    wx_kN_per_m: float,
    wy_kN_per_m: float,
    point_loads: list[LocalPointLoad],
    n_stations: int,
) -> list[StationForce]:
    """Internal forces from joint end actions and member loads, by equilibrium."""
    xs = {0.0, length_m}
    if n_stations >= 2:
        step = length_m / (n_stations - 1)
        for index in range(n_stations):
            xs.add(index * step)
    for point in point_loads:
        x = point.x_over_l * length_m
        gap = min(1.0e-6, length_m * 1.0e-9) if length_m else 1.0e-6
        xs.add(max(0.0, x - gap))
        xs.add(min(length_m, x + gap))
    for zero in _shear_zeros(length_m, end, wy_kN_per_m, point_loads):
        xs.add(zero)
    return [
        section_force_at(x, length_m, end, wx_kN_per_m, wy_kN_per_m, point_loads) for x in sorted(xs)
    ]


def section_force_at(
    x_m: float,
    length_m: float,
    end: EndActions,
    wx_kN_per_m: float,
    wy_kN_per_m: float,
    point_loads: list[LocalPointLoad],
) -> StationForce:
    """Internal N, V, M at ``x_m`` from node i.

    Point loads exactly at ``x_m`` are not included; use a station just to the
    right to see the shear after the load. Moment is continuous.
    """
    px = 0.0
    py = 0.0
    moment_from_points = 0.0
    for point in point_loads:
        xp = point.x_over_l * length_m
        if xp < x_m:
            px += point.px_kN
            py += point.py_kN
            moment_from_points += (x_m - xp) * point.py_kN
    axial = -end.fx_i_kN - wx_kN_per_m * x_m - px
    shear = end.fy_i_kN + wy_kN_per_m * x_m + py
    moment = (
        -end.m_i_kNm
        + x_m * end.fy_i_kN
        + wy_kN_per_m * x_m * x_m / 2.0
        + moment_from_points
    )
    return StationForce(x_m=x_m, n_kN=axial, v_kN=shear, m_kNm=moment)


def equilibrium_residual(
    model: Model,
    result: CaseResult,
    factors: dict[str, float],
) -> tuple[float, float, float]:
    """Residual global Fx (kN), Fy (kN), and moment about the origin (kN·m).

    A solved linear model should leave a residual near zero. Moment is
    counter-clockwise positive: ``x * Fy - y * Fx`` with x, y in metres.
    """
    fx = 0.0
    fy = 0.0
    mz = 0.0
    for node_id, reaction in result.reactions.items():
        node = model.nodes[node_id]
        x_m = node.x_mm / 1000.0
        y_m = node.y_mm / 1000.0
        fx += reaction.fx_kN
        fy += reaction.fy_kN
        mz += reaction.mz_kNm + x_m * reaction.fy_kN - y_m * reaction.fx_kN
    for load in model.nodal_loads:
        factor = factors.get(load.case, 0.0)
        if factor == 0.0:
            continue
        node = model.nodes[load.node]
        x_m = node.x_mm / 1000.0
        y_m = node.y_mm / 1000.0
        load_fx = factor * load.fx_n / 1000.0
        load_fy = factor * load.fy_n / 1000.0
        load_mz = factor * load.mz_nmm / 1.0e6
        fx += load_fx
        fy += load_fy
        mz += load_mz + x_m * load_fy - y_m * load_fx
    for udl in model.member_udls:
        factor = factors.get(udl.case, 0.0)
        if factor == 0.0:
            continue
        gx, gy = _udl_global_resultant_kN(model, udl, factor)
        mx, my = _udl_centroid_m(model, udl.member)
        fx += gx
        fy += gy
        mz += mx * gy - my * gx
    for point in model.member_point_loads:
        factor = factors.get(point.case, 0.0)
        if factor == 0.0:
            continue
        gx, gy, x_m, y_m = _point_global_kN(model, point, factor)
        fx += gx
        fy += gy
        mz += x_m * gy - y_m * gx
    return (fx, fy, mz)


def _superpose(
    model: Model,
    by_case: dict[str, CaseResult],
    combination: Combination,
    n_stations: int,
) -> CaseResult:
    displacements: dict[str, NodeDisplacement] = {}
    for node_id in model.nodes:
        ux = uy = rz = 0.0
        for case_id, factor in combination.factors.items():
            disp = by_case[case_id].displacements[node_id]
            ux += factor * disp.ux_mm
            uy += factor * disp.uy_mm
            rz += factor * disp.rz_rad
        displacements[node_id] = NodeDisplacement(ux, uy, rz)
    reactions: dict[str, NodeReaction] = {}
    for node_id in model.supports:
        fx = fy = mz = 0.0
        for case_id, factor in combination.factors.items():
            reaction = by_case[case_id].reactions[node_id]
            fx += factor * reaction.fx_kN
            fy += factor * reaction.fy_kN
            mz += factor * reaction.mz_kNm
        reactions[node_id] = NodeReaction(fx, fy, mz)
    members: dict[str, MemberForces] = {}
    for member_id, member in model.members.items():
        fx_i = fy_i = m_i = fx_j = fy_j = m_j = 0.0
        wx = wy = 0.0
        points: dict[float, list[float]] = {}
        for case_id, factor in combination.factors.items():
            source = by_case[case_id].members[member_id]
            fx_i += factor * source.end.fx_i_kN
            fy_i += factor * source.end.fy_i_kN
            m_i += factor * source.end.m_i_kNm
            fx_j += factor * source.end.fx_j_kN
            fy_j += factor * source.end.fy_j_kN
            m_j += factor * source.end.m_j_kNm
            wx += factor * source.wx_kN_per_m
            wy += factor * source.wy_kN_per_m
            for point in source.point_loads:
                bucket = points.setdefault(point.x_over_l, [0.0, 0.0])
                bucket[0] += factor * point.px_kN
                bucket[1] += factor * point.py_kN
        end = EndActions(fx_i, fy_i, m_i, fx_j, fy_j, m_j)
        local_points = [
            LocalPointLoad(x_over_l, pair[0], pair[1]) for x_over_l, pair in sorted(points.items())
        ]
        length_m = model.member_length_mm(member.i_node, member.j_node) / 1000.0
        members[member_id] = MemberForces(
            end=end,
            wx_kN_per_m=wx,
            wy_kN_per_m=wy,
            point_loads=local_points,
            stations=section_stations(length_m, end, wx, wy, local_points, n_stations),
            length_m=length_m,
        )
    combined = CaseResult(
        name=combination.name,
        displacements=displacements,
        reactions=reactions,
        members=members,
        expression=combination.expression,
        reference=combination.reference,
    )
    combined.equilibrium_kN = equilibrium_residual(model, combined, combination.factors)
    return combined


def _add_dummy_dof_if_needed(ops, model: Model, node_tag: dict[str, int]) -> None:
    """Give OpenSees at least one equation.

    A member with every end degree of freedom fixed has no free equation, and
    the solver rejects that system even though the fixed-end forces are well
    defined. A free dummy node tied to a fixed joint does not load the frame.
    """
    free = 0
    for node_id in model.nodes:
        support = model.supports.get(node_id)
        restrained = (
            (support.ux, support.uy, support.rz) if support is not None else (False, False, False)
        )
        free += sum(1 for flag in restrained if not flag)
    if free > 0:
        return
    first = next(iter(model.nodes))
    anchor = model.nodes[first]
    # Offset the dummy so the stabilizer element has a real length. It connects
    # to a fixed joint, so it does not change the frame.
    ops.node(_DUMMY_NODE, anchor.x_mm, anchor.y_mm + 1000.0)
    ops.element(
        "elasticBeamColumn",
        _DUMMY_ELEMENT,
        node_tag[first],
        _DUMMY_NODE,
        1.0,
        1.0,
        1.0,
        1,
    )


def _apply_nodal_loads(ops, model: Model, case_id: str, node_tag: dict[str, int]) -> None:
    totals: dict[str, list[float]] = {}
    for load in model.nodal_loads:
        if load.case != case_id:
            continue
        bucket = totals.setdefault(load.node, [0.0, 0.0, 0.0])
        bucket[0] += load.fx_n
        bucket[1] += load.fy_n
        bucket[2] += load.mz_nmm
    for node_id, (fx, fy, mz) in totals.items():
        ops.load(node_tag[node_id], fx, fy, mz)


def _member_local_loads(
    model: Model, member_id: str, case_id: str
) -> tuple[float, float, list[LocalPointLoad]]:
    _length, cosine, sine, _dx, _dy = model.member_geometry(member_id)
    wx = 0.0
    wy = 0.0
    for load in model.member_udls:
        if load.member != member_id or load.case != case_id:
            continue
        local_x, local_y = _to_local(load, cosine, sine)
        wx += local_x
        wy += local_y
    points: list[LocalPointLoad] = []
    for load in model.member_point_loads:
        if load.member != member_id or load.case != case_id:
            continue
        local_x, local_y = _to_local_point(load, cosine, sine)
        points.append(LocalPointLoad(load.x_over_l, n_to_kn(local_x), n_to_kn(local_y)))
    return wx, wy, points


def _to_local(load: MemberUDL, cosine: float, sine: float) -> tuple[float, float]:
    if load.system == "local":
        return load.wx, load.wy
    return _global_to_local(load.wx, load.wy, cosine, sine)


def _to_local_point(load: MemberPointLoad, cosine: float, sine: float) -> tuple[float, float]:
    if load.system == "local":
        return load.px, load.py
    return _global_to_local(load.px, load.py, cosine, sine)


def _global_to_local(fx: float, fy: float, cosine: float, sine: float) -> tuple[float, float]:
    """Project a global vector onto local x (i→j) and local y (90° CCW)."""
    return fx * cosine + fy * sine, -fx * sine + fy * cosine


def _local_to_global(local_x: float, local_y: float, cosine: float, sine: float) -> tuple[float, float]:
    return local_x * cosine + local_y * (-sine), local_x * sine + local_y * cosine


def _udl_global_resultant_kN(model: Model, load: MemberUDL, factor: float) -> tuple[float, float]:
    length_mm, cosine, sine, _dx, _dy = model.member_geometry(load.member)
    length_m = length_mm / 1000.0
    if load.system == "global":
        return factor * load.wx * length_m, factor * load.wy * length_m
    gx, gy = _local_to_global(load.wx, load.wy, cosine, sine)
    return factor * gx * length_m, factor * gy * length_m


def _udl_centroid_m(model: Model, member_id: str) -> tuple[float, float]:
    member = model.members[member_id]
    i_node = model.nodes[member.i_node]
    j_node = model.nodes[member.j_node]
    return (i_node.x_mm + j_node.x_mm) / 2000.0, (i_node.y_mm + j_node.y_mm) / 2000.0


def _point_global_kN(
    model: Model, load: MemberPointLoad, factor: float
) -> tuple[float, float, float, float]:
    member = model.members[load.member]
    i_node = model.nodes[member.i_node]
    _length, cosine, sine, dx, dy = model.member_geometry(load.member)
    x_m = (i_node.x_mm + load.x_over_l * dx) / 1000.0
    y_m = (i_node.y_mm + load.x_over_l * dy) / 1000.0
    if load.system == "global":
        return factor * n_to_kn(load.px), factor * n_to_kn(load.py), x_m, y_m
    gx, gy = _local_to_global(load.px, load.py, cosine, sine)
    return factor * n_to_kn(gx), factor * n_to_kn(gy), x_m, y_m


def _segments(length_m: float, point_loads: list[LocalPointLoad]) -> list[tuple[float, float]]:
    cuts = [0.0, length_m]
    for point in point_loads:
        cuts.append(point.x_over_l * length_m)
    cuts = sorted(set(cuts))
    return [(cuts[i], cuts[i + 1]) for i in range(len(cuts) - 1) if cuts[i + 1] > cuts[i]]


def _shear_zeros(
    length_m: float,
    end: EndActions,
    wy: float,
    point_loads: list[LocalPointLoad],
) -> list[float]:
    """Stations where shear crosses zero. Moment peaks there."""
    if abs(wy) < 1.0e-12 or length_m <= 0:
        return []
    zeros: list[float] = []
    for x0, x1 in _segments(length_m, point_loads):
        py = 0.0
        for point in point_loads:
            xp = point.x_over_l * length_m
            if xp <= x0 + 1.0e-12:
                py += point.py_kN
        x = -(end.fy_i_kN + py) / wy
        if x0 + 1.0e-9 < x < x1 - 1.0e-9:
            zeros.append(x)
    return zeros
