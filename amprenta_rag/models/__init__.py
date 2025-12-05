"""
Domain models for the multi-omics platform.

This package provides unified domain models that serve as:
1. Stable abstraction layer across the codebase
2. Foundation for Postgres schema design
3. Type-safe data structures

The models work with both Notion (current) and Postgres (future).
"""

from amprenta_rag.models.domain import (
    Dataset,
    Experiment,
    Feature,
    FeatureType,
    OmicsType,
    Program,
    Signature,
    SignatureComponent,
    SignatureDirection,
    feature_type_from_string,
    omics_type_from_string,
)

__all__ = [
    "Dataset",
    "Experiment",
    "Feature",
    "FeatureType",
    "OmicsType",
    "Program",
    "Signature",
    "SignatureComponent",
    "SignatureDirection",
    "feature_type_from_string",
    "omics_type_from_string",
]

