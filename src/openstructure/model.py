"""Structural model: nodes, members, sections, loads, and combinations.

Coordinates are stored in millimetres. The methods that build a model take
metres, kilonewtons, and kN/m, which is the unit system used in input files.
"""

from __future__ import annotations

from dataclasses import dataclass, field

from openstructure.errors import OpenStructureError
from openstructure.units import kn_to_n, knm_to_nmm, m_to_mm

LOAD_KINDS = (
    "dead",
    "live",
    "roof_live",
    "snow",
    "rain",
    "wind",
    "earthquake",
    "other",
)
DIRECTIONAL_KINDS = ("wind", "earthquake")


@dataclass
class Material:
    """Concrete (and, later, steel) material.

    ``ec_mpa`` is the modulus used by the linear-elastic analysis. Reinforcement
    yield ``fy_mpa`` is not used by the elastic analysis; design routines read it.
    """

    id: str
    ec_mpa: float
    fc_mpa: float | None = None
    fy_mpa: float | None = None


@dataclass
class RectangularSection:
    """Gross rectangular section for in-plane bending.

    ``b_mm`` is the width perpendicular to the frame (the width of the concrete
    compression block). ``h_mm`` is the dimension along the element local y-axis,
    which is the depth that enters ``I = b h^3 / 12``.

    For a horizontal beam whose local y-axis points upward, ``h_mm`` is the
    vertical depth. For a vertical column whose i-node is below the j-node,
    local y points in the negative global x-direction, so ``h_mm`` is the
    horizontal in-plane dimension. A square column is unaffected.
    """

    id: str
    material: str
    b_mm: float
    h_mm: float
    i_factor: float = 1.0

    @property
    def area_mm2(self) -> float:
        return self.b_mm * self.h_mm

    @property
    def inertia_mm4(self) -> float:
        return self.i_factor * self.b_mm * self.h_mm**3 / 12.0


@dataclass
class Node:
    id: str
    x_mm: float
    y_mm: float


@dataclass
class Support:
    """True means the degree of freedom is restrained. Global axes, y upward."""

    node: str
    ux: bool = False
    uy: bool = False
    rz: bool = False


@dataclass
class Member:
    id: str
    i_node: str
    j_node: str
    section: str


@dataclass
class LoadCase:
    """A linear load case.

    Wind and earthquake cases are reversible by default: preset combinations
    apply the pattern in both senses. Set ``reversible`` false when the case
    already points in the only direction that should be used.
    """

    id: str
    kind: str = "other"
    reversible: bool = False


@dataclass
class NodalLoad:
    node: str
    case: str
    fx_n: float = 0.0
    fy_n: float = 0.0
    mz_nmm: float = 0.0


@dataclass
class MemberUDL:
    """Uniform load over the full member.

    ``system`` is ``global`` or ``local``. Components are force per length in
    N/mm (numerically equal to kN/m). In the local system, x runs from node i
    to node j and y is 90° counter-clockwise from x.
    """

    member: str
    case: str
    system: str
    wx: float = 0.0
    wy: float = 0.0


@dataclass
class MemberPointLoad:
    """Point load on a member. ``x_over_l`` is the fraction of length from node i."""

    member: str
    case: str
    system: str
    x_over_l: float
    px: float = 0.0
    py: float = 0.0


@dataclass
class Combination:
    """Linear combination of load-case results.

    Superposition is valid for the linear-static engine. A later nonlinear
    engine should apply ``factors`` to the loads and analyse the combination
    as one pattern instead of superposing results.
    """

    name: str
    factors: dict[str, float]
    expression: str = ""
    reference: str = ""


@dataclass
class Model:
    """A 2D frame in the global x–y plane. Global y is upward. Moments are CCW."""

    title: str = ""
    description: str = ""
    materials: dict[str, Material] = field(default_factory=dict)
    sections: dict[str, RectangularSection] = field(default_factory=dict)
    nodes: dict[str, Node] = field(default_factory=dict)
    supports: dict[str, Support] = field(default_factory=dict)
    members: dict[str, Member] = field(default_factory=dict)
    load_cases: dict[str, LoadCase] = field(default_factory=dict)
    nodal_loads: list[NodalLoad] = field(default_factory=list)
    member_udls: list[MemberUDL] = field(default_factory=list)
    member_point_loads: list[MemberPointLoad] = field(default_factory=list)
    combinations: dict[str, Combination] = field(default_factory=dict)

    def add_material(
        self,
        id: str,
        ec_mpa: float,
        fc_mpa: float | None = None,
        fy_mpa: float | None = None,
    ) -> Material:
        self._new_id(id, self.materials, "material")
        if ec_mpa <= 0:
            raise OpenStructureError(f"Material {id}: Ec must be positive (MPa).")
        material = Material(id=id, ec_mpa=ec_mpa, fc_mpa=fc_mpa, fy_mpa=fy_mpa)
        self.materials[id] = material
        return material

    def add_rectangular_section(
        self,
        id: str,
        material: str,
        b_mm: float,
        h_mm: float,
        i_factor: float = 1.0,
    ) -> RectangularSection:
        self._new_id(id, self.sections, "section")
        if material not in self.materials:
            raise OpenStructureError(f"Section {id}: unknown material {material!r}.")
        if b_mm <= 0 or h_mm <= 0:
            raise OpenStructureError(f"Section {id}: b and h must be positive (mm).")
        if i_factor <= 0:
            raise OpenStructureError(f"Section {id}: i_factor must be positive.")
        section = RectangularSection(
            id=id, material=material, b_mm=b_mm, h_mm=h_mm, i_factor=i_factor
        )
        self.sections[id] = section
        return section

    def add_node(self, id: str, x_m: float, y_m: float) -> Node:
        self._new_id(id, self.nodes, "node")
        node = Node(id=id, x_mm=m_to_mm(x_m), y_mm=m_to_mm(y_m))
        self.nodes[id] = node
        return node

    def add_support(
        self,
        node: str,
        ux: bool = False,
        uy: bool = False,
        rz: bool = False,
        *,
        fixed: bool = False,
        pinned: bool = False,
        roller_x: bool = False,
    ) -> Support:
        """Restrain a node.

        ``roller_x`` restrains global y only (the node can move in x).
        ``pinned`` restrains x and y. ``fixed`` restrains x, y, and rotation.
        """
        self._require_node(node)
        if node in self.supports:
            raise OpenStructureError(f"Node {node} already has a support.")
        if fixed:
            ux, uy, rz = True, True, True
        elif pinned:
            ux, uy, rz = True, True, False
        elif roller_x:
            ux, uy, rz = False, True, False
        support = Support(node=node, ux=ux, uy=uy, rz=rz)
        self.supports[node] = support
        return support

    def add_member(self, id: str, i_node: str, j_node: str, section: str) -> Member:
        self._new_id(id, self.members, "member")
        self._require_node(i_node)
        self._require_node(j_node)
        if i_node == j_node:
            raise OpenStructureError(f"Member {id}: the two nodes must differ.")
        if section not in self.sections:
            raise OpenStructureError(f"Member {id}: unknown section {section!r}.")
        if self.member_length_mm(i_node, j_node) <= 0:
            raise OpenStructureError(f"Member {id}: length is zero.")
        member = Member(id=id, i_node=i_node, j_node=j_node, section=section)
        self.members[id] = member
        return member

    def add_load_case(self, id: str, kind: str = "other", reversible: bool | None = None) -> LoadCase:
        self._new_id(id, self.load_cases, "load case")
        if kind not in LOAD_KINDS:
            raise OpenStructureError(
                f"Load case {id}: kind {kind!r} is not one of {', '.join(LOAD_KINDS)}."
            )
        if reversible is None:
            reversible = kind in DIRECTIONAL_KINDS
        case = LoadCase(id=id, kind=kind, reversible=reversible)
        self.load_cases[id] = case
        return case

    def add_nodal_load(
        self,
        node: str,
        case: str,
        fx_kN: float = 0.0,
        fy_kN: float = 0.0,
        mz_kNm: float = 0.0,
    ) -> NodalLoad:
        self._require_node(node)
        self._require_case(case)
        load = NodalLoad(
            node=node,
            case=case,
            fx_n=kn_to_n(fx_kN),
            fy_n=kn_to_n(fy_kN),
            mz_nmm=knm_to_nmm(mz_kNm),
        )
        self.nodal_loads.append(load)
        return load

    def add_member_udl_global(
        self,
        member: str,
        case: str,
        wx_kN_per_m: float = 0.0,
        wy_kN_per_m: float = 0.0,
    ) -> MemberUDL:
        """Uniform load in global axes. Positive wy is upward. kN/m equals N/mm."""
        self._require_member(member)
        self._require_case(case)
        load = MemberUDL(
            member=member,
            case=case,
            system="global",
            wx=wx_kN_per_m,
            wy=wy_kN_per_m,
        )
        self.member_udls.append(load)
        return load

    def add_member_udl_local(
        self,
        member: str,
        case: str,
        wx_kN_per_m: float = 0.0,
        wy_kN_per_m: float = 0.0,
    ) -> MemberUDL:
        """Uniform load in member local axes. Positive local y is 90° CCW from i→j."""
        self._require_member(member)
        self._require_case(case)
        load = MemberUDL(
            member=member,
            case=case,
            system="local",
            wx=wx_kN_per_m,
            wy=wy_kN_per_m,
        )
        self.member_udls.append(load)
        return load

    def add_member_point_load_global(
        self,
        member: str,
        case: str,
        x_over_l: float,
        fx_kN: float = 0.0,
        fy_kN: float = 0.0,
    ) -> MemberPointLoad:
        self._check_member_point(member, case, x_over_l)
        load = MemberPointLoad(
            member=member,
            case=case,
            system="global",
            x_over_l=x_over_l,
            px=kn_to_n(fx_kN),
            py=kn_to_n(fy_kN),
        )
        self.member_point_loads.append(load)
        return load

    def add_member_point_load_local(
        self,
        member: str,
        case: str,
        x_over_l: float,
        px_kN: float = 0.0,
        py_kN: float = 0.0,
    ) -> MemberPointLoad:
        self._check_member_point(member, case, x_over_l)
        load = MemberPointLoad(
            member=member,
            case=case,
            system="local",
            x_over_l=x_over_l,
            px=kn_to_n(px_kN),
            py=kn_to_n(py_kN),
        )
        self.member_point_loads.append(load)
        return load

    def add_combination(
        self,
        name: str,
        factors: dict[str, float],
        expression: str = "",
        reference: str = "",
    ) -> Combination:
        if name in self.combinations:
            raise OpenStructureError(f"Combination {name!r} already exists.")
        if not factors:
            raise OpenStructureError(f"Combination {name}: no load cases.")
        for case_id in factors:
            self._require_case(case_id)
        combination = Combination(
            name=name,
            factors=dict(factors),
            expression=expression or _format_factors(factors),
            reference=reference,
        )
        self.combinations[name] = combination
        return combination

    def member_length_mm(self, i_node: str, j_node: str) -> float:
        a = self.nodes[i_node]
        b = self.nodes[j_node]
        return ((b.x_mm - a.x_mm) ** 2 + (b.y_mm - a.y_mm) ** 2) ** 0.5

    def member_geometry(self, member_id: str) -> tuple[float, float, float, float, float]:
        """Return length (mm), cos, sin of the i→j axis, and (dx, dy) in mm."""
        member = self.members[member_id]
        a = self.nodes[member.i_node]
        b = self.nodes[member.j_node]
        dx = b.x_mm - a.x_mm
        dy = b.y_mm - a.y_mm
        length = (dx * dx + dy * dy) ** 0.5
        return length, dx / length, dy / length, dx, dy

    def validate(self) -> None:
        if not self.nodes:
            raise OpenStructureError("The model has no nodes.")
        if not self.members:
            raise OpenStructureError("The model has no members.")
        for member in self.members.values():
            section = self.sections[member.section]
            material = self.materials[section.material]
            if material.ec_mpa <= 0:
                raise OpenStructureError(f"Material {material.id}: Ec must be positive.")
            if self.member_length_mm(member.i_node, member.j_node) <= 0:
                raise OpenStructureError(f"Member {member.id}: length is zero.")

    def _new_id(self, id: str, store: dict, what: str) -> None:
        if not id or not str(id).strip():
            raise OpenStructureError(f"A {what} needs an id.")
        if id in store:
            raise OpenStructureError(f"Duplicate {what} id {id!r}.")

    def _require_node(self, node: str) -> None:
        if node not in self.nodes:
            raise OpenStructureError(f"Unknown node {node!r}.")

    def _require_member(self, member: str) -> None:
        if member not in self.members:
            raise OpenStructureError(f"Unknown member {member!r}.")

    def _require_case(self, case: str) -> None:
        if case not in self.load_cases:
            raise OpenStructureError(f"Unknown load case {case!r}.")

    def _check_member_point(self, member: str, case: str, x_over_l: float) -> None:
        self._require_member(member)
        self._require_case(case)
        if not 0.0 < x_over_l < 1.0:
            raise OpenStructureError(
                f"Member point load on {member}: x_over_l must be between 0 and 1 "
                "(nodal loads belong on nodes)."
            )


def _format_factors(factors: dict[str, float]) -> str:
    parts = []
    for case_id, factor in factors.items():
        parts.append(f"{factor:g}·{case_id}")
    return " + ".join(parts)
