"""The example frame runs from the file and from the command line."""

import json
from pathlib import Path

from openstructure.cli import main
from openstructure.pipeline import run_file

EXAMPLE = Path(__file__).resolve().parents[1] / "examples" / "two_bay_two_storey.yaml"


def test_example_frame_equilibrium_and_design(tmp_path):
    output = run_file(EXAMPLE, tmp_path, plots=True)
    assert output.model.title.startswith("Two-bay")
    for case in output.analysis.by_combination.values():
        fx, fy, mz = case.equilibrium_kN
        assert abs(fx) < 1e-4
        assert abs(fy) < 1e-4
        assert abs(mz) < 1e-3
    # Gravity moment in the first-floor beam is sagging somewhere in the span.
    gravity = output.analysis.by_combination["EIT-1.4D+1.7L"].members["B10"]
    assert max(station.m_kNm for station in gravity.stations) > 0
    statuses = {packet.status for packet in output.designs}
    assert statuses <= {"OK", "NG", "REVIEW"}
    assert any(packet.kind == "beam_flexure" and packet.code_id == "EIT1008" for packet in output.designs)
    assert any(packet.kind == "column_pm" and packet.code_id == "ACI318-19" for packet in output.designs)
    report = (tmp_path / "report.md").read_text(encoding="utf-8")
    assert "Whitney" in report or "φMn" in report
    payload = json.loads((tmp_path / "results.json").read_text(encoding="utf-8"))
    assert payload["units"]["force"] == "kN"
    assert any((tmp_path / "diagrams").glob("*.png"))


def test_cli_run(tmp_path):
    code = main(["run", str(EXAMPLE), "--out", str(tmp_path), "--no-plots"])
    assert code == 0
    assert (tmp_path / "results.json").is_file()
    assert (tmp_path / "report.md").is_file()
