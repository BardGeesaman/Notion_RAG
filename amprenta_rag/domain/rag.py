from typing import List, Optional

from pydantic import BaseModel


class RAGChunkMetadata(BaseModel):
    source_type: str
    source_id: str
    tags: Optional[List[str]] = None
    semantic_labels: Optional[dict] = None
    extra: Optional[dict] = None


class RAGChunkCreate(BaseModel):
    content: str
    metadata: RAGChunkMetadata


class RAGQueryRequest(BaseModel):
    query: str
    top_k: Optional[int] = 10
    filters: Optional[dict] = None


class RAGQueryResult(BaseModel):
    chunk_id: str
    content: str
    score: Optional[float] = None
    metadata: Optional[RAGChunkMetadata] = None
