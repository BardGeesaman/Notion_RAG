"""
Semantic metadata extraction from Notion pages.

This package extracts structured metadata from Notion pages for different
source types (Literature, Email, Experiments, Datasets), including:
- Disease classifications
- Matrix types
- Model systems
- Lipid signature relations
- Phenotype axes
- Other semantic annotations

Maintains backward compatibility by re-exporting all public functions.
"""

from __future__ import annotations

from amprenta_rag.ingestion.metadata.dataset_extraction import (
    get_dataset_semantic_metadata,
)
from amprenta_rag.ingestion.metadata.email_extraction import (
    get_email_semantic_metadata,
)
from amprenta_rag.ingestion.metadata.experiment_extraction import (
    get_experiment_semantic_metadata,
)
from amprenta_rag.ingestion.metadata.literature_extraction import (
    get_literature_semantic_metadata,
)

__all__ = [
    "get_literature_semantic_metadata",
    "get_email_semantic_metadata",
    "get_experiment_semantic_metadata",
    "get_dataset_semantic_metadata",
]

