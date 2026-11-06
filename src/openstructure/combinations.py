"""Ultimate load-combination presets for ACI 318-19 and EIT-era strength design.

Factors are applied to load cases by ``LoadCase.kind``. Cases of the same kind
are all included. Wind and earthquake are not applied together. Each reversible
wind or earthquake case is expanded in both senses.

The live-load factor in ACI 318-19 §5.3.1(c), (d), and (e) is kept at 1.0.
Section 5.3.3 allows 0.5 for some occupancies; that reduction is not applied
automatically.
"""

from __future__ import annotations

from dataclasses import dataclass

from openstructure.errors import OpenStructureError
from openstructure.model import DIRECTIONAL_KINDS, Combination, Model

ACI318 = "ACI318-19"
EIT1008 = "EIT1008"


@dataclass(frozen=True)
class _Template:
    name: str
    expression: str
    reference: str
    factors_by_kind: dict[str, float]
    # If this kind is absent the whole template is skipped.
    requires_kind: str | None = None
    # If any of these kinds exist, skip this template (a more specific one is used).
    skip_if_kinds: tuple[str, ...] = ()
    # Kinds that may be absent; the term is then omitted rather than rejecting the template.
    optional_kinds: tuple[str, ...] = ()
    directional_kind: str | None = None
    assumption: bool = False


def _aci_templates() -> list[_Template]:
    ref = "ACI 318-19 Section 5.3.1"
    live = 1.0  # §5.3.3 reduction is not automatic; see module docstring.
    return [
        _Template("ACI-1.4D", "1.4D", f"{ref}(a)", {"dead": 1.4}),
        _Template(
            "ACI-1.2D+1.6L",
            "1.2D + 1.6L",
            f"{ref}(b), with no roof live, snow, or rain",
            {"dead": 1.2, "live": 1.6},
            skip_if_kinds=("roof_live", "snow", "rain"),
        ),
        _Template(
            "ACI-1.2D+1.6L+0.5Lr",
            "1.2D + 1.6L + 0.5Lr",
            f"{ref}(b)",
            {"dead": 1.2, "live": 1.6, "roof_live": 0.5},
            requires_kind="roof_live",
        ),
        _Template(
            "ACI-1.2D+1.6L+0.5S",
            "1.2D + 1.6L + 0.5S",
            f"{ref}(b)",
            {"dead": 1.2, "live": 1.6, "snow": 0.5},
            requires_kind="snow",
        ),
        _Template(
            "ACI-1.2D+1.6L+0.5R",
            "1.2D + 1.6L + 0.5R",
            f"{ref}(b)",
            {"dead": 1.2, "live": 1.6, "rain": 0.5},
            requires_kind="rain",
        ),
        _Template(
            "ACI-1.2D+1.6Lr+1.0L",
            "1.2D + 1.6Lr + 1.0L",
            f"{ref}(c); L factor is 1.0 (§5.3.3 not applied)",
            {"dead": 1.2, "roof_live": 1.6, "live": live},
            requires_kind="roof_live",
            optional_kinds=("live",),
        ),
        _Template(
            "ACI-1.2D+1.6S+1.0L",
            "1.2D + 1.6S + 1.0L",
            f"{ref}(c); L factor is 1.0 (§5.3.3 not applied)",
            {"dead": 1.2, "snow": 1.6, "live": live},
            requires_kind="snow",
            optional_kinds=("live",),
        ),
        _Template(
            "ACI-1.2D+1.6R+1.0L",
            "1.2D + 1.6R + 1.0L",
            f"{ref}(c); L factor is 1.0 (§5.3.3 not applied)",
            {"dead": 1.2, "rain": 1.6, "live": live},
            requires_kind="rain",
            optional_kinds=("live",),
        ),
        _Template(
            "ACI-1.2D+1.0L+1.0W",
            "1.2D + 1.0W + 1.0L",
            f"{ref}(d) without Lr, S, or R; L factor is 1.0",
            {"dead": 1.2, "live": live, "wind": 1.0},
            requires_kind="wind",
            skip_if_kinds=("roof_live", "snow", "rain"),
            optional_kinds=("live",),
            directional_kind="wind",
        ),
        _Template(
            "ACI-1.2D+1.0L+1.0W+0.5Lr",
            "1.2D + 1.0W + 1.0L + 0.5Lr",
            f"{ref}(d)",
            {"dead": 1.2, "live": live, "wind": 1.0, "roof_live": 0.5},
            requires_kind="roof_live",
            optional_kinds=("live",),
            directional_kind="wind",
        ),
        _Template(
            "ACI-1.2D+1.0L+1.0W+0.5S",
            "1.2D + 1.0W + 1.0L + 0.5S",
            f"{ref}(d)",
            {"dead": 1.2, "live": live, "wind": 1.0, "snow": 0.5},
            requires_kind="snow",
            optional_kinds=("live",),
            directional_kind="wind",
        ),
        _Template(
            "ACI-1.2D+1.0L+1.0W+0.5R",
            "1.2D + 1.0W + 1.0L + 0.5R",
            f"{ref}(d)",
            {"dead": 1.2, "live": live, "wind": 1.0, "rain": 0.5},
            requires_kind="rain",
            optional_kinds=("live",),
            directional_kind="wind",
        ),
        _Template(
            "ACI-1.2D+1.0L+1.0E",
            "1.2D + 1.0E + 1.0L",
            f"{ref}(e) without snow; L factor is 1.0. E is the case as defined by the user.",
            {"dead": 1.2, "live": live, "earthquake": 1.0},
            requires_kind="earthquake",
            skip_if_kinds=("snow",),
            optional_kinds=("live",),
            directional_kind="earthquake",
        ),
        _Template(
            "ACI-1.2D+1.0L+1.0E+0.2S",
            "1.2D + 1.0E + 1.0L + 0.2S",
            f"{ref}(e)",
            {"dead": 1.2, "live": live, "earthquake": 1.0, "snow": 0.2},
            requires_kind="snow",
            optional_kinds=("live",),
            directional_kind="earthquake",
        ),
        _Template(
            "ACI-0.9D+1.0W",
            "0.9D + 1.0W",
            f"{ref}(f)",
            {"dead": 0.9, "wind": 1.0},
            requires_kind="wind",
            directional_kind="wind",
        ),
        _Template(
            "ACI-0.9D+1.0E",
            "0.9D + 1.0E",
            f"{ref}(g)",
            {"dead": 0.9, "earthquake": 1.0},
            requires_kind="earthquake",
            directional_kind="earthquake",
        ),
    ]


def _eit_templates() -> list[_Template]:
    # 0.75 * (1.4, 1.7, 1.7) and the 1.1E substitution used with ACI 318-99.
    d_w = 0.75 * 1.4
    l_w = 0.75 * 1.7
    w_w = 0.75 * 1.7
    e_in_w = 0.75 * 1.7 * 1.1  # 1.87 = 1.1 * 1.7, then * 0.75
    e_in_09 = 1.3 * 1.1
    gravity_ref = (
        "Strength-design gravity combination U = 1.4D + 1.7L used with "
        "วสท. 1008 (Thai Concrete Association journal comparing ว.ส.ท. 1008 and ACI 318-05; "
        "also the worked column examples in standard วสท. 1008 teaching notes)."
    )
    wind_ref = (
        "ASSUMPTION, confirm in the printed วสท. 1008: wind combination taken from the "
        "ACI 318-99 Section 9.2 form that this generation of strength design used, "
        "U = 0.75(1.4D + 1.7L + 1.7W) and U = 0.9D + 1.3W. Not given an EIT clause number here."
    )
    seismic_ref = (
        "ASSUMPTION, confirm in the printed วสท. 1008 and the applicable ministerial "
        "regulation: earthquake is introduced by replacing W with 1.1E in the ACI 318-99 "
        "wind combinations, giving 0.75(1.4D + 1.7L + 1.87E) and 0.9D + 1.43E. "
        "Not given an EIT clause number here."
    )
    return [
        _Template("EIT-1.4D+1.7L", "1.4D + 1.7L", gravity_ref, {"dead": 1.4, "live": 1.7}),
        _Template(
            "EIT-0.75(1.4D+1.7L+1.7W)",
            "0.75(1.4D + 1.7L + 1.7W)",
            wind_ref,
            {"dead": d_w, "live": l_w, "wind": w_w},
            requires_kind="wind",
            optional_kinds=("live",),
            directional_kind="wind",
            assumption=True,
        ),
        _Template(
            "EIT-0.9D+1.3W",
            "0.9D + 1.3W",
            wind_ref,
            {"dead": 0.9, "wind": 1.3},
            requires_kind="wind",
            directional_kind="wind",
            assumption=True,
        ),
        _Template(
            "EIT-0.75(1.4D+1.7L+1.87E)",
            "0.75(1.4D + 1.7L + 1.87E)",
            seismic_ref,
            {"dead": d_w, "live": l_w, "earthquake": e_in_w},
            requires_kind="earthquake",
            optional_kinds=("live",),
            directional_kind="earthquake",
            assumption=True,
        ),
        _Template(
            "EIT-0.9D+1.43E",
            "0.9D + 1.43E",
            seismic_ref,
            {"dead": 0.9, "earthquake": e_in_09},
            requires_kind="earthquake",
            directional_kind="earthquake",
            assumption=True,
        ),
    ]


PRESETS: dict[str, list[_Template]] = {
    ACI318: _aci_templates(),
    EIT1008: _eit_templates(),
}


def apply_preset(model: Model, preset: str, *, replace: bool = False) -> list[str]:
    """Add preset combinations to ``model`` and return their names."""
    if preset not in PRESETS:
        known = ", ".join(sorted(PRESETS))
        raise OpenStructureError(f"Unknown combination preset {preset!r}. Known presets: {known}.")
    if replace:
        model.combinations.clear()
    names: list[str] = []
    for template in PRESETS[preset]:
        for combination in _expand(model, template):
            if combination.name in model.combinations:
                raise OpenStructureError(
                    f"Combination {combination.name!r} already exists. "
                    "Rename the explicit combination or build presets first."
                )
            model.combinations[combination.name] = combination
            names.append(combination.name)
    if not names:
        raise OpenStructureError(
            f"Preset {preset} produced no combinations. "
            "Add at least a dead or live load case before applying it."
        )
    return names


def _cases_of_kind(model: Model, kind: str) -> list[str]:
    return [case.id for case in model.load_cases.values() if case.kind == kind]


def _expand(model: Model, template: _Template) -> list[Combination]:
    if template.requires_kind and not _cases_of_kind(model, template.requires_kind):
        return []
    if any(_cases_of_kind(model, kind) for kind in template.skip_if_kinds):
        return []
    for kind in template.factors_by_kind:
        if kind in template.optional_kinds:
            continue
        if not _cases_of_kind(model, kind):
            return []
    active = {
        kind: factor
        for kind, factor in template.factors_by_kind.items()
        if _cases_of_kind(model, kind)
    }
    if not active:
        return []

    directional_ids: list[str] = []
    if template.directional_kind:
        directional_ids = _cases_of_kind(model, template.directional_kind)
        if not directional_ids:
            return []

    if not directional_ids:
        factors = _factors_for(model, active, directional_id=None, sign=1.0)
        if not factors:
            return []
        return [
            Combination(
                name=template.name,
                factors=factors,
                expression=template.expression,
                reference=template.reference,
            )
        ]

    combinations: list[Combination] = []
    for case_id in directional_ids:
        signs = (1.0, -1.0) if model.load_cases[case_id].reversible else (1.0,)
        for sign in signs:
            factors = _factors_for(model, active, directional_id=case_id, sign=sign)
            sense = "+" if sign > 0 else "-"
            name = f"{template.name}[{sense}{case_id}]"
            expression = template.expression.replace(
                "W" if template.directional_kind == "wind" else "E",
                f"{sense}{case_id}",
            )
            combinations.append(
                Combination(
                    name=name,
                    factors=factors,
                    expression=expression,
                    reference=template.reference,
                )
            )
    return combinations


def _factors_for(
    model: Model,
    active: dict[str, float],
    directional_id: str | None,
    sign: float,
) -> dict[str, float]:
    factors: dict[str, float] = {}
    for kind, factor in active.items():
        case_ids = _cases_of_kind(model, kind)
        if kind in DIRECTIONAL_KINDS:
            if directional_id is None:
                continue
            factors[directional_id] = factor * sign
            continue
        for case_id in case_ids:
            factors[case_id] = factor
    return factors
