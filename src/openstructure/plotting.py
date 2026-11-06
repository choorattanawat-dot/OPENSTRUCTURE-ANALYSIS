"""Matplotlib diagrams of internal forces.

Positive moment is drawn on the side that is in tension for a horizontal
member when local y is upward (sagging moment plotted downward of the axis).
The figure title states the sign convention so the plot can be read without
this module.
"""

from __future__ import annotations

from pathlib import Path

from openstructure.analysis import MemberForces


def plot_member_forces(
    member_id: str,
    forces: MemberForces,
    combination: str,
    path: Path | None = None,
):
    """Draw axial, shear, and moment diagrams. Returns the matplotlib figure."""
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    xs = [station.x_m for station in forces.stations]
    axial = [station.n_kN for station in forces.stations]
    shear = [station.v_kN for station in forces.stations]
    moment = [station.m_kNm for station in forces.stations]

    fig, axes = plt.subplots(3, 1, figsize=(8.0, 8.0), sharex=True)
    _panel(axes[0], xs, axial, "Axial force N (kN)", "positive = tension")
    _panel(axes[1], xs, shear, "Shear V (kN)", "positive = upward on a left-to-right member at the left end")
    _panel(
        axes[2],
        xs,
        [-value for value in moment],
        "Moment, plotted on the tension side (kN·m)",
        "sagging positive in the data; curve is drawn toward the tension face",
    )
    axes[2].set_xlabel("Distance from node i (m)")
    fig.suptitle(
        f"{member_id}  ·  {combination}\n"
        "Local x from i to j. Local y is 90° counter-clockwise from x."
    )
    fig.tight_layout()
    if path is not None:
        path.parent.mkdir(parents=True, exist_ok=True)
        fig.savefig(path, dpi=120)
    return fig


def _panel(axis, xs: list[float], values: list[float], title: str, note: str) -> None:
    axis.axhline(0.0, color="0.4", linewidth=0.8)
    axis.fill_between(xs, values, 0.0, color="0.75")
    axis.plot(xs, values, color="0.1", linewidth=1.4)
    axis.set_ylabel(title)
    axis.set_title(note, fontsize=9, loc="left")
    axis.grid(True, linewidth=0.3)
