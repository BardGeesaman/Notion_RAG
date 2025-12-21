"""
Compatibility wrapper for the refactored feature linking pipeline.

This module maintains backward compatibility by re-exporting all public APIs
from the new modular structure.

Scripts and other code importing from linking.py will continue to work
without changes.
"""

from __future__ import annotations

import importlib
from typing import Any

_gl = importlib.import_module("amprenta_rag.ingestion.features.general_linking")
_ml = importlib.import_module("amprenta_rag.ingestion.features.metabolite_linking")

_add_dataset_relation: Any = getattr(_gl, "_add_dataset_relation", None)
_find_or_create_feature_page: Any = getattr(_gl, "_find_or_create_feature_page", None)
link_feature: Any = getattr(_gl, "link_feature", None)

_add_relation_to_metabolite_page: Any = getattr(_ml, "_add_relation_to_metabolite_page", None)
_find_or_create_metabolite_page: Any = getattr(_ml, "_find_or_create_metabolite_page", None)
link_features_to_notion_items: Any = getattr(_ml, "link_features_to_notion_items", None)

__all__ = [
    "_find_or_create_metabolite_page",
    "_add_relation_to_metabolite_page",
    "_find_or_create_feature_page",
    "_add_dataset_relation",
    "link_feature",
    "link_features_to_notion_items",
]
