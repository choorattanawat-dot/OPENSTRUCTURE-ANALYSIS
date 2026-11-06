"""The example notebook's code cells run as written."""

import json
from pathlib import Path


def test_example_notebook_cells_run(monkeypatch):
    notebook = Path(__file__).resolve().parents[1] / "examples" / "two_bay_two_storey.ipynb"
    nb = json.loads(notebook.read_text(encoding="utf-8"))
    monkeypatch.chdir(notebook.parent)
    namespace: dict = {"__name__": "__notebook__"}
    for cell in nb["cells"]:
        if cell["cell_type"] != "code":
            continue
        source = "".join(cell["source"])
        exec(compile(source, "<notebook>", "exec"), namespace)
    assert "1281" in namespace["flexure"].report_markdown
    import matplotlib.pyplot as plt

    plt.close("all")
