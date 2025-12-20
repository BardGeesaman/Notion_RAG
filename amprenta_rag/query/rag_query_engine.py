# amprenta_rag/query/rag_query_engine.py

"""
Compatibility wrapper for the refactored query engine.

This module maintains backward compatibility by re-exporting all public APIs
from the new modular structure (pinecone_query.py and rag_engine.py).

Scripts and other code importing from rag_query_engine.py will continue to work
without changes.
"""

from __future__ import annotations

# Re-export Pinecone query functions from pinecone_query
from amprenta_rag.query.pinecone_query import (build_meta_filter, embed_query,
                                               query_pinecone)
# Re-export dataclasses and main query function from rag_engine
from amprenta_rag.query.rag_engine import (MatchSummary, RAGQueryResult,
                                           query_rag, signature_similarity_query)
# Re-export cross-omics reasoning functions
from amprenta_rag.query.cross_omics.dataset_summary_postgres import cross_omics_dataset_summary_postgres
from amprenta_rag.query.cross_omics.feature_summary_postgres import cross_omics_feature_summary_postgres
from amprenta_rag.query.cross_omics.program_summary_postgres import cross_omics_program_summary_postgres
from amprenta_rag.query.cross_omics.signature_summary_postgres import cross_omics_signature_summary_postgres

__all__ = [
    "MatchSummary",
    "RAGQueryResult",
    "query_rag",
    "signature_similarity_query",
    "cross_omics_program_summary_postgres",
    "cross_omics_signature_summary_postgres",
    "cross_omics_feature_summary_postgres",
    "cross_omics_dataset_summary_postgres",
    "build_meta_filter",
    "embed_query",
    "query_pinecone",
]
