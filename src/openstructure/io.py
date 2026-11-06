"""Read a model from YAML or JSON.

The computation core does not import this module. File keys use the units in
their names: metres for coordinates, millimetres for sections, kN, kN/m, kN·m,
and MPa.
"""

from __future__ import annotations

import json
from pathlib import Path
from typing import Any

from openstructure.errors import OpenStructureError
from openstructure.model import Model

try:
    import yaml
except ImportError as exc:  # pragma: no cover
    yaml = None
    _YAML_ERROR = exc
else:
    _YAML_ERROR = None


def load_model(path: str | Path) -> tuple[Model, dict[str, Any]]:
    """Return the model and the raw mapping (needed for analysis and design blocks)."""
    file_path = Path(path)
    if not file_path.is_file():
        raise OpenStructureError(f"Model file not found: {file_path}")
    text = file_path.read_text(encoding="utf-8")
    suffix = file_path.suffix.lower()
    if suffix in (".yaml", ".yml"):
        if yaml is None:
            raise OpenStructureError(f"PyYAML is required to read {file_path}.") from _YAML_ERROR
        raw = yaml.safe_load(text)
    elif suffix == ".json":
        raw = json.loads(text)
    else:
        raise OpenStructureError("Model file must be .yaml, .yml, or .json.")
    if not isinstance(raw, dict):
        raise OpenStructureError("The model file must be a mapping at the top level.")
    return _build_model(raw), raw


def _build_model(raw: dict[str, Any]) -> Model:
    meta = raw.get("metadata") or {}
    model = Model(title=str(meta.get("title", "")), description=str(meta.get("description", "")))
    for item in raw.get("materials") or []:
        fc = item.get("fc_mpa")
        ec = item.get("ec_mpa")
        if ec is None:
            if fc is None:
                raise OpenStructureError(f"Material {item.get('id')}: give ec_mpa or fc_mpa.")
            ec = 4700.0 * float(fc) ** 0.5
        model.add_material(
            str(item["id"]),
            ec_mpa=float(ec),
            fc_mpa=None if fc is None else float(fc),
            fy_mpa=None if item.get("fy_mpa") is None else float(item["fy_mpa"]),
        )
    for item in raw.get("sections") or []:
        if item.get("type", "rectangular") != "rectangular":
            raise OpenStructureError(
                f"Section {item.get('id')}: v0.1 supports rectangular sections only."
            )
        model.add_rectangular_section(
            str(item["id"]),
            str(item["material"]),
            b_mm=float(item["b_mm"]),
            h_mm=float(item["h_mm"]),
            i_factor=float(item.get("i_factor", 1.0)),
        )
    for item in raw.get("nodes") or []:
        model.add_node(str(item["id"]), float(item["x_m"]), float(item["y_m"]))
    for item in raw.get("supports") or []:
        model.add_support(
            str(item["node"]),
            ux=bool(item.get("ux", False)),
            uy=bool(item.get("uy", False)),
            rz=bool(item.get("rz", False)),
            fixed=bool(item.get("fixed", False)),
            pinned=bool(item.get("pinned", False)),
            roller_x=bool(item.get("roller_x", False)),
        )
    for item in raw.get("members") or []:
        model.add_member(str(item["id"]), str(item["i"]), str(item["j"]), str(item["section"]))
    for item in raw.get("load_cases") or []:
        model.add_load_case(
            str(item["id"]),
            kind=str(item.get("kind", "other")),
            reversible=item.get("reversible"),
        )
    loads = raw.get("loads") or {}
    for item in loads.get("nodal") or []:
        model.add_nodal_load(
            str(item["node"]),
            str(item["case"]),
            fx_kN=float(item.get("fx_kN", 0.0)),
            fy_kN=float(item.get("fy_kN", 0.0)),
            mz_kNm=float(item.get("mz_kNm", 0.0)),
        )
    for item in loads.get("member_udl") or []:
        system = str(item.get("system", "global"))
        wx = float(item.get("wx_kN_per_m", 0.0))
        wy = float(item.get("wy_kN_per_m", 0.0))
        if system == "global":
            model.add_member_udl_global(str(item["member"]), str(item["case"]), wx, wy)
        elif system == "local":
            model.add_member_udl_local(str(item["member"]), str(item["case"]), wx, wy)
        else:
            raise OpenStructureError(f"member_udl system must be global or local, not {system!r}.")
    for item in loads.get("member_point") or []:
        system = str(item.get("system", "global"))
        kwargs = {
            "x_over_l": float(item["x_over_l"]),
        }
        if system == "global":
            model.add_member_point_load_global(
                str(item["member"]),
                str(item["case"]),
                fx_kN=float(item.get("fx_kN", 0.0)),
                fy_kN=float(item.get("fy_kN", 0.0)),
                **kwargs,
            )
        elif system == "local":
            model.add_member_point_load_local(
                str(item["member"]),
                str(item["case"]),
                px_kN=float(item.get("px_kN", 0.0)),
                py_kN=float(item.get("py_kN", 0.0)),
                **kwargs,
            )
        else:
            raise OpenStructureError(f"member_point system must be global or local, not {system!r}.")

    from openstructure.combinations import apply_preset

    combos = raw.get("combinations") or {}
    preset = combos.get("preset")
    if preset:
        apply_preset(model, str(preset))
    for item in combos.get("explicit") or []:
        model.add_combination(
            str(item["name"]),
            {str(key): float(value) for key, value in item["factors"].items()},
            expression=str(item.get("expression", "")),
            reference=str(item.get("reference", "")),
        )
    return model
