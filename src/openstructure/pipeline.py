"""Run a model file: linear analysis, diagrams, and the design checks it asks for.

Design demands are read from the analysis result. The design routines
themselves stay free of file and plotting code.
"""

from __future__ import annotations

import json
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any

from openstructure.analysis import (
    AnalysisOptions,
    AnalysisResult,
    CaseResult,
    section_force_at,
)
from openstructure.design import (
    check_column_uniaxial,
    design_beam_flexure,
    design_beam_shear,
    get_code,
    perimeter_bars,
)
from openstructure.design.column import ColumnResult
from openstructure.design.flexure import FlexureResult
from openstructure.design.shear import ShearResult
from openstructure.errors import OpenStructureError
from openstructure.io import load_model
from openstructure.model import Model
from openstructure.plotting import plot_member_forces


@dataclass
class DesignPacket:
    check_id: str
    code_id: str
    kind: str
    status: str
    summary: str
    report_markdown: str


@dataclass
class RunOutput:
    model: Model
    analysis: AnalysisResult
    designs: list[DesignPacket] = field(default_factory=list)
    warnings: list[str] = field(default_factory=list)


def run_file(path: str | Path, out_dir: str | Path | None = None, plots: bool = True) -> RunOutput:
    model, raw = load_model(path)
    options = _options(raw.get("analysis") or {})
    engine_name = options.kind
    if engine_name != "linear_static":
        raise OpenStructureError(
            f"Analysis kind {engine_name!r} is not implemented. Use linear_static."
        )
    from openstructure.analysis import LinearStaticEngine

    analysis = LinearStaticEngine().analyze(model, options)
    designs = _run_designs(model, analysis, raw.get("design") or [])
    output = RunOutput(model=model, analysis=analysis, designs=designs)
    if out_dir is not None:
        _write(output, Path(out_dir), plots)
    return output


def _options(block: dict[str, Any]) -> AnalysisOptions:
    return AnalysisOptions(
        kind=str(block.get("kind", "linear_static")),
        geometric_transformation=str(block.get("transformation", "Linear")),
        n_stations=int(block.get("n_stations", 21)),
    )


def _codes(entry: dict[str, Any]) -> list:
    raw = entry.get("code", "ACI318-19")
    names = raw if isinstance(raw, list) else [raw]
    quality = str(entry.get("quality_control", "good"))
    return [get_code(str(name), quality_control=quality) for name in names]


def _combinations(entry: dict[str, Any], analysis: AnalysisResult) -> list[str]:
    requested = entry.get("combinations", "all")
    if requested == "all":
        if not analysis.by_combination:
            raise OpenStructureError(f"Design check {entry.get('id')}: the model has no combinations.")
        return list(analysis.by_combination)
    if isinstance(requested, str):
        return [requested]
    return [str(name) for name in requested]


def _run_designs(model: Model, analysis: AnalysisResult, checks: list[dict[str, Any]]) -> list[DesignPacket]:
    packets: list[DesignPacket] = []
    for entry in checks:
        kind = str(entry.get("type", ""))
        if kind == "beam_flexure":
            packets.extend(_flexure(model, analysis, entry))
        elif kind == "beam_shear":
            packets.extend(_shear(model, analysis, entry))
        elif kind == "column_pm":
            packets.extend(_column(model, analysis, entry))
        else:
            raise OpenStructureError(
                f"Unknown design type {kind!r}. Use beam_flexure, beam_shear, or column_pm."
            )
    return packets


def _section_material(model: Model, member_id: str) -> tuple[float, float, float, float]:
    member = model.members[member_id]
    section = model.sections[member.section]
    material = model.materials[section.material]
    if material.fc_mpa is None or material.fy_mpa is None:
        raise OpenStructureError(
            f"Material {material.id} needs fc_mpa and fy_mpa for design of member {member_id}."
        )
    return section.b_mm, section.h_mm, material.fc_mpa, material.fy_mpa


def _flexure(model: Model, analysis: AnalysisResult, entry: dict[str, Any]) -> list[DesignPacket]:
    member_id = str(entry["member"])
    b_mm, h_mm, fc, fy = _section_material(model, member_id)
    b_mm = float(entry.get("b_mm", b_mm))
    h_mm = float(entry.get("h_mm", h_mm))
    d_mm = float(entry["d_mm"])
    d_prime = entry.get("d_prime_mm")
    d_prime_mm = None if d_prime is None else float(d_prime)
    packets: list[DesignPacket] = []
    for code in _codes(entry):
        candidates: list[tuple[float, str, FlexureResult]] = []
        lines = [
            f"### {entry.get('id', member_id)} flexure — {code.code_name}",
            "",
            "Positive moment is sagging (compression on local +y). "
            "Negative moment is hogging. Mu is the extreme centre-line moment, "
            "including the value at the node; it is not reduced to the column face. "
            "Each row is one load combination.",
            "",
            "| Combination | Face | Mu (kN·m) | As (mm²) | As' (mm²) | φMn (kN·m) | Status |",
            "|---|---|---:|---:|---:|---:|---|",
        ]
        for name in _combinations(entry, analysis):
            forces = analysis.by_combination[name].members[member_id]
            low, high = forces.extreme_moment()
            for label, station in (("sagging", high), ("hogging", low)):
                if abs(station.m_kNm) < 1.0e-8:
                    continue
                result = design_beam_flexure(
                    station.m_kNm,
                    b_mm,
                    d_mm,
                    h_mm,
                    fc,
                    fy,
                    code,
                    d_prime_mm=d_prime_mm,
                )
                candidates.append((result.as_mm2, f"{name} / {label}", result))
                lines.append(
                    f"| {name} | {label} | {station.m_kNm:.2f} | {result.as_mm2:.1f} | "
                    f"{result.as_prime_mm2:.1f} | {result.phi_mn_kNm:.2f} | {result.status} |"
                )
        if not candidates:
            raise OpenStructureError(f"Flexure check {entry.get('id')}: no moment demand.")
        _as, label, governing = max(candidates, key=lambda item: (item[0], abs(item[2].phi_mn_kNm)))
        lines.extend(["", f"Governing case for tension steel: {label}.", "", governing.report.to_markdown()])
        status = "NG" if any(item[2].status == "NG" for item in candidates) else governing.status
        packets.append(
            DesignPacket(
                check_id=str(entry.get("id", member_id)),
                code_id=code.code_id,
                kind="beam_flexure",
                status=status,
                summary=f"{label}: As = {governing.as_mm2:.1f} mm², status {governing.status}",
                report_markdown="\n".join(lines),
            )
        )
    return packets


def _shear(model: Model, analysis: AnalysisResult, entry: dict[str, Any]) -> list[DesignPacket]:
    member_id = str(entry["member"])
    b_mm, _h, fc, fy = _section_material(model, member_id)
    b_mm = float(entry.get("b_mm", b_mm))
    d_mm = float(entry["d_mm"])
    face = float(entry.get("face_offset_mm", 0.0))
    diameter = float(entry.get("stirrup_diameter_mm", 9.0))
    legs = int(entry.get("legs", 2))
    step = float(entry.get("spacing_step_mm", 25.0))
    fyt = float(entry.get("fyt_mpa", fy))
    packets: list[DesignPacket] = []
    for code in _codes(entry):
        rows = [
            f"### {entry.get('id', member_id)} shear — {code.code_name}",
            "",
            "Vu is taken at `face_offset_mm + d` from each member end (the node). "
            "Set face_offset_mm to half the support depth when the model uses centre-lines. "
            "v0.1 does not do this offset by itself.",
            "",
            "| Combination | x (m) | Vu (kN) | s (mm) | φVn (kN) | Status |",
            "|---|---:|---:|---:|---:|---|",
        ]
        chosen: tuple[float, str, ShearResult] | None = None
        for name in _combinations(entry, analysis):
            forces = analysis.by_combination[name].members[member_id]
            station, note = _critical_shear(forces, d_mm, face)
            result = design_beam_shear(
                abs(station.v_kN),
                b_mm,
                d_mm,
                fc,
                fyt,
                code,
                stirrup_diameter_mm=diameter,
                legs=legs,
                spacing_step_mm=step,
            )
            rows.append(
                f"| {name} | {station.x_m:.3f} | {abs(station.v_kN):.2f} | "
                f"{result.s_provided_mm:.0f} | {result.phi_vn_kN:.2f} | {result.status} |"
            )
            if note:
                rows.append(f"|  |  | {note} |  |  |  |")
            if chosen is None or abs(station.v_kN) > chosen[0]:
                chosen = (abs(station.v_kN), name, result)
        assert chosen is not None
        rows.extend(["", f"Full calculation for the largest Vu ({chosen[1]}).", "", chosen[2].report.to_markdown()])
        packets.append(
            DesignPacket(
                check_id=str(entry.get("id", member_id)),
                code_id=code.code_id,
                kind="beam_shear",
                status=chosen[2].status,
                summary=(
                    f"{chosen[1]}: Vu = {chosen[0]:.2f} kN, "
                    f"s = {chosen[2].s_provided_mm:.0f} mm, status {chosen[2].status}"
                ),
                report_markdown="\n".join(rows),
            )
        )
    return packets


def _critical_shear(forces, d_mm: float, face_offset_mm: float):
    x1 = (face_offset_mm + d_mm) / 1000.0
    x2 = forces.length_m - x1
    if x1 <= 0 or x2 <= x1:
        return forces.extreme_shear(), "d and the face offset cover the span; max |V| was used."
    left = section_force_at(
        x1, forces.length_m, forces.end, forces.wx_kN_per_m, forces.wy_kN_per_m, forces.point_loads
    )
    right = section_force_at(
        x2, forces.length_m, forces.end, forces.wx_kN_per_m, forces.wy_kN_per_m, forces.point_loads
    )
    return (left if abs(left.v_kN) >= abs(right.v_kN) else right), ""


def _column(model: Model, analysis: AnalysisResult, entry: dict[str, Any]) -> list[DesignPacket]:
    member_id = str(entry["member"])
    b_mm, h_mm, fc, fy = _section_material(model, member_id)
    b_mm = float(entry.get("b_mm", b_mm))
    h_mm = float(entry.get("h_mm", h_mm))
    bars = _bars(entry, b_mm, h_mm)
    tied = str(entry.get("confinement", "tied")) != "spiral"
    packets: list[DesignPacket] = []
    for code in _codes(entry):
        rows = [
            f"### {entry.get('id', member_id)} column P–M — {code.code_name}",
            "",
            "Pu is compression positive, converted from the analysis axial force "
            "(analysis tension is positive, so Pu = −N). Both ends are checked.",
            "",
            "| Combination | End | Pu (kN) | Mu (kN·m) | Utilization | Status |",
            "|---|---|---:|---:|---:|---|",
        ]
        chosen: tuple[float, str, ColumnResult] | None = None
        for name in _combinations(entry, analysis):
            forces = analysis.by_combination[name].members[member_id]
            for x_m, end_name in ((0.0, "i"), (forces.length_m, "j")):
                station = section_force_at(
                    x_m,
                    forces.length_m,
                    forces.end,
                    forces.wx_kN_per_m,
                    forces.wy_kN_per_m,
                    forces.point_loads,
                )
                pu = -station.n_kN
                mu = station.m_kNm
                result = check_column_uniaxial(pu, mu, b_mm, h_mm, bars, fc, fy, code, tied=tied)
                rows.append(
                    f"| {name} | {end_name} | {pu:.2f} | {mu:.2f} | {result.utilization:.3f} | {result.status} |"
                )
                if chosen is None or result.utilization > chosen[0]:
                    chosen = (result.utilization, f"{name} end {end_name}", result)
        assert chosen is not None
        rows.extend(
            ["", f"Full calculation for the governing end ({chosen[1]}).", "", chosen[2].report.to_markdown()]
        )
        packets.append(
            DesignPacket(
                check_id=str(entry.get("id", member_id)),
                code_id=code.code_id,
                kind="column_pm",
                status=chosen[2].status,
                summary=f"{chosen[1]}: utilization {chosen[0]:.3f}, status {chosen[2].status}",
                report_markdown="\n".join(rows),
            )
        )
    return packets


def _bars(entry: dict[str, Any], b_mm: float, h_mm: float) -> list[tuple[float, float]]:
    reinforcement = entry.get("reinforcement")
    if reinforcement:
        layout = reinforcement.get("layout", "perimeter")
        if layout != "perimeter":
            raise OpenStructureError("v0.1 column layouts: 'perimeter' or an explicit bars list.")
        return perimeter_bars(
            b_mm,
            h_mm,
            int(reinforcement["n_along_b"]),
            int(reinforcement["n_along_h"]),
            float(reinforcement["bar_diameter_mm"]),
            float(reinforcement["cover_to_center_mm"]),
        )
    explicit = entry.get("bars") or []
    if not explicit:
        raise OpenStructureError("Column check needs reinforcement or a bars list.")
    return [(float(bar["area_mm2"]), float(bar["y_mm"])) for bar in explicit]


def _write(output: RunOutput, out_dir: Path, plots: bool) -> None:
    out_dir.mkdir(parents=True, exist_ok=True)
    (out_dir / "results.json").write_text(
        json.dumps(_result_payload(output), indent=2),
        encoding="utf-8",
    )
    report = ["# OPENSTRUCTURE-ANALYSIS results", "", f"## {output.model.title}", ""]
    report.append("Units: displacements in mm, forces in kN, moments in kN·m, stress in MPa.")
    report.append("")
    report.append("## Equilibrium residual (applied loads + reactions)")
    report.append("")
    report.append("| Combination | Fx (kN) | Fy (kN) | M about origin (kN·m) |")
    report.append("|---|---:|---:|---:|")
    for name, case in output.analysis.by_combination.items():
        fx, fy, mz = case.equilibrium_kN
        report.append(f"| {name} | {fx:.4g} | {fy:.4g} | {mz:.4g} |")
    report.append("")
    for packet in output.designs:
        report.append(packet.report_markdown)
        report.append("")
    (out_dir / "report.md").write_text("\n".join(report), encoding="utf-8")
    if plots:
        for name, case in output.analysis.by_combination.items():
            for member_id, forces in case.members.items():
                safe = _safe(f"{member_id}_{name}")
                plot_member_forces(member_id, forces, name, out_dir / "diagrams" / f"{safe}.png")
                import matplotlib.pyplot as plt

                plt.close("all")


def _result_payload(output: RunOutput) -> dict[str, Any]:
    def case_payload(case: CaseResult) -> dict[str, Any]:
        return {
            "expression": case.expression,
            "reference": case.reference,
            "equilibrium_residual": {
                "fx_kN": case.equilibrium_kN[0],
                "fy_kN": case.equilibrium_kN[1],
                "mz_kNm": case.equilibrium_kN[2],
                "note": "Sum of reactions and applied loads. Near zero means global equilibrium.",
            },
            "displacements_mm_rad": {
                node: {"ux_mm": d.ux_mm, "uy_mm": d.uy_mm, "rz_rad": d.rz_rad}
                for node, d in case.displacements.items()
            },
            "reactions": {
                node: {"fx_kN": r.fx_kN, "fy_kN": r.fy_kN, "mz_kNm": r.mz_kNm}
                for node, r in case.reactions.items()
            },
            "members": {
                member_id: {
                    "length_m": forces.length_m,
                    "end_joint_on_member_local": {
                        "fx_i_kN": forces.end.fx_i_kN,
                        "fy_i_kN": forces.end.fy_i_kN,
                        "m_i_kNm": forces.end.m_i_kNm,
                        "fx_j_kN": forces.end.fx_j_kN,
                        "fy_j_kN": forces.end.fy_j_kN,
                        "m_j_kNm": forces.end.m_j_kNm,
                    },
                    "stations": [
                        {"x_m": s.x_m, "N_kN": s.n_kN, "V_kN": s.v_kN, "M_kNm": s.m_kNm}
                        for s in forces.stations
                    ],
                    "sign_convention": {
                        "N_kN": "positive tension",
                        "V_kN": "positive when the left side pushes the right side toward local +y",
                        "M_kNm": "positive when it compresses the local +y face (sagging on a left-to-right beam)",
                    },
                }
                for member_id, forces in case.members.items()
            },
        }

    return {
        "title": output.model.title,
        "units": {
            "displacement": "mm",
            "rotation": "rad",
            "force": "kN",
            "moment": "kN.m",
            "stress": "MPa",
            "length_in_diagrams": "m",
        },
        "load_cases": {name: case_payload(case) for name, case in output.analysis.by_case.items()},
        "combinations": {
            name: case_payload(case) for name, case in output.analysis.by_combination.items()
        },
        "design": [
            {
                "id": packet.check_id,
                "code": packet.code_id,
                "type": packet.kind,
                "status": packet.status,
                "summary": packet.summary,
            }
            for packet in output.designs
        ],
    }


def _safe(name: str) -> str:
    return "".join(char if char.isalnum() or char in "-_" else "_" for char in name)
