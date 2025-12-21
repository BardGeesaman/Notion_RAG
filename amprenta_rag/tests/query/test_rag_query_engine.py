from __future__ import annotations

from amprenta_rag.query import rag_query_engine as rqe

def test_rag_query_engine_reexports():
    # Verify all expected symbols are available
    assert rqe.query_rag is not None
    assert rqe.signature_similarity_query is not None
    assert rqe.cross_omics_program_summary_postgres is not None
    assert rqe.cross_omics_signature_summary_postgres is not None
    assert rqe.cross_omics_feature_summary_postgres is not None
    assert rqe.cross_omics_dataset_summary_postgres is not None
    assert rqe.build_meta_filter is not None
    assert rqe.embed_query is not None
    assert rqe.query_pinecone is not None
    assert rqe.MatchSummary is not None
    assert rqe.RAGQueryResult is not None

