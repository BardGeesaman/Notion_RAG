"""
Biophysical assay data processing module.

This module provides parsers and analysis tools for various biophysical techniques:
- SPR (Surface Plasmon Resonance): Biacore instruments
- MST (Microscale Thermophoresis): NanoTemper instruments  
- DSC (Differential Scanning Calorimetry): MicroCal instruments

Each technique has dedicated parsers for common instrument file formats.
"""

from .spr_parser import (
    SPRData,
    Sensorgram,
    parse_biacore_csv,
    parse_biacore_sensorgram_txt,
    extract_sensorgrams,
    detect_injection_phases,
    validate_array_size,
)
from .mst_parser import (
    MSTData,
    DosePoint,
    parse_nanotemper_xlsx,
    parse_mst_csv,
    extract_dose_response,
    calculate_fnorm,
)
from .dsc_parser import (
    DSCData,
    Scan,
    parse_microcal_csv,
    parse_ta_instruments_txt,
    baseline_correction,
)

__all__ = [
    # SPR
    "SPRData",
    "Sensorgram", 
    "parse_biacore_csv",
    "parse_biacore_sensorgram_txt",
    "extract_sensorgrams",
    "detect_injection_phases",
    "validate_array_size",
    # MST
    "MSTData",
    "DosePoint",
    "parse_nanotemper_xlsx",
    "parse_mst_csv",
    "extract_dose_response",
    "calculate_fnorm",
    # DSC
    "DSCData",
    "Scan",
    "parse_microcal_csv",
    "parse_ta_instruments_txt",
    "baseline_correction",
]
