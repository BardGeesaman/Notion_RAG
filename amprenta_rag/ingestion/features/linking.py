"""
Compatibility wrapper for the refactored feature linking pipeline.

This module maintains backward compatibility by re-exporting all public APIs
from the new modular structure.

Scripts and other code importing from linking.py will continue to work
without changes.
"""

from __future__ import annotations

from typing import Any

from amprenta_rag.ingestion.features import general_linking as _gl  # type: ignore[import]
from amprenta_rag.ingestion.features import metabolite_linking as _ml  # type: ignore[import]

# Best-effort attribute lookups; silence missing attrs for mypy
_add_dataset_relation: Any = getattr(_gl, "_add_dataset_relation", None)  # type: ignore[attr-defined]
_find_or_create_feature_page: Any = getattr(_gl, "_find_or_create_feature_page", None)  # type: ignore[attr-defined]
link_feature: Any = getattr(_gl, "link_feature", None)  # type: ignore[attr-defined]

_add_relation_to_metabolite_page: Any = getattr(_ml, "_add_relation_to_metabolite_page", None)  # type: ignore[attr-defined]
_find_or_create_metabolite_page: Any = getattr(_ml, "_find_or_create_metabolite_page", None)  # type: ignore[attr-defined]
link_features_to_notion_items: Any = getattr(_ml, "link_features_to_notion_items", None)  # type: ignore[attr-defined]

__all__ = [
    "_find_or_create_metabolite_page",
    "_add_relation_to_metabolite_page",
    "_find_or_create_feature_page",
    "_add_dataset_relation",
    "link_feature",
    "link_features_to_notion_items",
]
