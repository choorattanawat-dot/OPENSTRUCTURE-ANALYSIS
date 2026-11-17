"""Regression against ETABS for a real project model.

Project models are not stored in this (public) repository. Point the two environment
variables at the files to run it:
    OPENSTRUCTURE_REF_E2K     = model text export (.e2k)
    OPENSTRUCTURE_REF_TABLES  = ETABS result tables exported as text
Reference building (ETABS 23.1.1, kgf-m, PT flat slab 250, columns 40x80, core wall).
"""
import os

import pytest

E2K = os.environ.get("OPENSTRUCTURE_REF_E2K")
TABLES = os.environ.get("OPENSTRUCTURE_REF_TABLES")
pytestmark = pytest.mark.skipif(not (E2K and TABLES and os.path.exists(E2K) and os.path.exists(TABLES)),
                                reason="reference model not available")


def test_matches_etabs():
    from openstructure import Model
    from openstructure.verify import verify
    v = verify(Model.from_e2k(E2K), TABLES, 6)
    rows = {r["item"]: r for r in v["rows"]}
    for pat in ("DL", "SDL", "LL"):
        assert abs(rows[f"Base FZ {pat}"]["diff_pct"]) < 0.1
    for k in ("T1", "T2", "T3", "T4"):
        assert abs(rows[k]["diff_pct"]) < 1.0
    for k in ("T5", "T6"):
        assert abs(rows[k]["diff_pct"]) < 3.5
