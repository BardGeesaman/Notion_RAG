from amprenta_rag.domain.rag import (
    RAGChunkCreate,
    RAGChunkMetadata,
    RAGQueryRequest,
    RAGQueryResult,
)


def test_chunk_metadata_defaults():
    meta = RAGChunkMetadata(source_type="doc", source_id="s1")
    assert meta.tags is None
    assert meta.extra is None


def test_chunk_create_and_result():
    meta = RAGChunkMetadata(source_type="doc", source_id="s1", tags=["t1"])
    chunk = RAGChunkCreate(content="text", metadata=meta)
    assert chunk.metadata.tags == ["t1"]

    res = RAGQueryResult(chunk_id="c1", content="text", score=0.9, metadata=meta)
    assert res.metadata.source_id == "s1"


def test_query_request_defaults():
    req = RAGQueryRequest(query="hello")
    assert req.top_k == 10
    assert req.filters is None

