"""OPENSTRUCTURE-ANALYSIS: 2D frame analysis and RC member design.

The computation core (``model``, ``analysis``, ``design``) does not read files
or print reports. Input/output lives in ``io``, ``pipeline``, ``plotting``
and ``cli`` so a later interface can call the same core.

Internal analysis units are newtons and millimetres. Public results use
kilonewtons, metres, and kN·m, and every result object names its units.
"""

from openstructure.analysis import AnalysisOptions, LinearStaticEngine
from openstructure.combinations import apply_preset
from openstructure.design import (
    ACI318,
    EIT1008,
    check_column_uniaxial,
    design_beam_flexure,
    design_beam_shear,
    get_code,
)
from openstructure.model import Model

__version__ = "0.1.0"

__all__ = [
    "ACI318",
    "AnalysisOptions",
    "EIT1008",
    "LinearStaticEngine",
    "Model",
    "__version__",
    "apply_preset",
    "check_column_uniaxial",
    "design_beam_flexure",
    "design_beam_shear",
    "get_code",
]
