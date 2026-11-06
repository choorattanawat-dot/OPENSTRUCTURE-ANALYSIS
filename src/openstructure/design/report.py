"""Step-by-step calculation reports.

A report is the record a reviewing engineer reads. Each step keeps the formula,
the substitution with units, and the result. Nothing in the design routines
should hide an intermediate value that the step list does not show.
"""

from __future__ import annotations

from dataclasses import dataclass, field


def fmt(value: float, digits: int = 3) -> str:
    if value is None:
        return "—"
    magnitude = abs(value)
    if magnitude != 0 and (magnitude >= 1.0e5 or magnitude < 1.0e-3):
        return f"{value:.4g}"
    return f"{value:.{digits}f}"


@dataclass
class CalcStep:
    title: str
    formula: str
    substitution: str
    result: str
    reference: str = ""
    note: str = ""


@dataclass
class CalcReport:
    title: str
    status: str
    code: str
    steps: list[CalcStep] = field(default_factory=list)
    assumptions: list[str] = field(default_factory=list)
    summary: list[str] = field(default_factory=list)

    def add(
        self,
        title: str,
        formula: str,
        substitution: str,
        result: str,
        reference: str = "",
        note: str = "",
    ) -> None:
        self.steps.append(
            CalcStep(
                title=title,
                formula=formula,
                substitution=substitution,
                result=result,
                reference=reference,
                note=note,
            )
        )

    def to_markdown(self) -> str:
        lines = [
            f"# {self.title}",
            "",
            f"**Status:** {self.status}",
            "",
            f"**Code:** {self.code}",
            "",
        ]
        if self.summary:
            lines.append("## Summary")
            lines.append("")
            for item in self.summary:
                lines.append(f"- {item}")
            lines.append("")
        for index, step in enumerate(self.steps, start=1):
            lines.append(f"## Step {index}. {step.title}")
            lines.append("")
            if step.reference:
                lines.append(f"Reference: {step.reference}")
                lines.append("")
            lines.append("Formula:")
            lines.append("")
            lines.append(f"    {step.formula}")
            lines.append("")
            lines.append("Substitution:")
            lines.append("")
            lines.append(f"    {step.substitution}")
            lines.append("")
            lines.append(f"Result: {step.result}")
            lines.append("")
            if step.note:
                lines.append(step.note)
                lines.append("")
        if self.assumptions:
            lines.append("## Assumptions to confirm")
            lines.append("")
            for item in self.assumptions:
                lines.append(f"- {item}")
            lines.append("")
        return "\n".join(lines)
