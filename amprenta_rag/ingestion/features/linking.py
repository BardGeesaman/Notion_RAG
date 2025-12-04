"""
Compatibility wrapper for the refactored feature linking pipeline.

This module maintains backward compatibility by re-exporting all public APIs
from the new modular structure.

Scripts and other code importing from linking.py will continue to work
without changes.
"""

from __future__ import annotations

from amprenta_rag.ingestion.features.general_linking import (
    _add_dataset_relation,
    _find_or_create_feature_page,
    link_feature,
)
from amprenta_rag.ingestion.features.metabolite_linking import (
    _add_relation_to_metabolite_page,
    _find_or_create_metabolite_page,
    link_features_to_notion_items,
)

__all__ = [
    "_find_or_create_metabolite_page",
    "_add_relation_to_metabolite_page",
    "_find_or_create_feature_page",
    "_add_dataset_relation",
    "link_feature",
    "link_features_to_notion_items",
]
