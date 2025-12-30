"""
Flow Cytometry / FACS data processing module.

This module provides comprehensive tools for flow cytometry data analysis including:
- FCS file parsing and validation
- Data transformations (logicle, arcsinh)  
- Compensation matrix application
- Event subsampling
- Parquet storage for performance

Key components:
- fcs_parser: FCS file I/O and metadata extraction
- transforms: Mathematical transformations for flow data
"""

from amprenta_rag.flow_cytometry.fcs_parser import (
    load_fcs,
    extract_metadata,
    validate_fcs,
    save_events_parquet,
    load_events_parquet,
    FlowMetadata,
)

from amprenta_rag.flow_cytometry.transforms import (
    logicle_transform,
    arcsinh_transform,
    apply_compensation,
    auto_logicle_params,
    subsample_events,
)

__all__ = [
    # FCS Parser
    "load_fcs",
    "extract_metadata", 
    "validate_fcs",
    "save_events_parquet",
    "load_events_parquet",
    "FlowMetadata",
    
    # Transforms
    "logicle_transform",
    "arcsinh_transform", 
    "apply_compensation",
    "auto_logicle_params",
    "subsample_events",
]
