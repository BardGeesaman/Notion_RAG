"""Vector store abstraction layer.

Supports both Pinecone and Postgres/pgvector backends via a shared interface.
"""

from __future__ import annotations

from typing import Any, Dict, List, Optional, Protocol, Sequence, TypedDict, cast

import sqlalchemy as sa

from amprenta_rag.config import get_config


class VectorMatch(TypedDict, total=False):
    id: str
    score: float
    metadata: Dict[str, Any]


class VectorRecord(TypedDict, total=False):
    id: str
    values: List[float]
    metadata: Dict[str, Any]


PineconeFilter = Dict[str, Any]


class VectorStore(Protocol):
    """Backend-agnostic vector store interface."""

    def query(
        self,
        vector: Sequence[float],
        top_k: int = 10,
        filter: Optional[PineconeFilter] = None,  # noqa: A002
        namespace: Optional[str] = None,
    ) -> List[VectorMatch]:
        ...

    def upsert(self, vectors: List[VectorRecord], namespace: Optional[str] = None) -> None:
        ...

    def delete(
        self,
        filter: Optional[PineconeFilter] = None,  # noqa: A002
        ids: Optional[List[str]] = None,
        namespace: Optional[str] = None,
    ) -> None:
        ...


class PineconeStore:
    """Pinecone-backed implementation of VectorStore."""

    def __init__(self, index: Any | None = None) -> None:
        if index is None:
            raise RuntimeError(
                "Pinecone backend is deprecated. Please use pgvector backend instead. "
                "Set VECTOR_BACKEND=pgvector in your environment."
            )
        self._index = index

    def query(
        self,
        vector: Sequence[float],
        top_k: int = 10,
        filter: Optional[PineconeFilter] = None,  # noqa: A002
        namespace: Optional[str] = None,
    ) -> List[VectorMatch]:
        cfg = get_config()
        ns = namespace or cfg.pinecone.namespace
        kwargs: Dict[str, Any] = {
            "vector": list(vector),
            "top_k": top_k,
            "include_metadata": True,
            "namespace": ns,
        }
        if filter:
            kwargs["filter"] = filter

        res = self._index.query(**kwargs)
        matches = getattr(res, "matches", None) or (res.get("matches", []) if isinstance(res, dict) else [])

        out: List[VectorMatch] = []
        for m in matches:
            if isinstance(m, dict):
                out.append(
                    {
                        "id": cast(str, m.get("id", "")),
                        "score": float(m.get("score", 0.0)),
                        "metadata": cast(Dict[str, Any], m.get("metadata", {}) or {}),
                    }
                )
            else:
                out.append(
                    {
                        "id": cast(str, getattr(m, "id", "")),
                        "score": float(getattr(m, "score", 0.0)),
                        "metadata": cast(Dict[str, Any], getattr(m, "metadata", {}) or {}),
                    }
                )
        return out

    def upsert(self, vectors: List[VectorRecord], namespace: Optional[str] = None) -> None:
        cfg = get_config()
        ns = namespace or cfg.pinecone.namespace
        self._index.upsert(vectors=vectors, namespace=ns)

    def delete(
        self,
        filter: Optional[PineconeFilter] = None,  # noqa: A002
        ids: Optional[List[str]] = None,
        namespace: Optional[str] = None,
    ) -> None:
        cfg = get_config()
        ns = namespace or cfg.pinecone.namespace
        kwargs: Dict[str, Any] = {"namespace": ns}
        if filter:
            kwargs["filter"] = filter
        if ids:
            kwargs["ids"] = ids
        self._index.delete(**kwargs)


def _filter_to_sql(where: Optional[PineconeFilter]) -> tuple[str, Dict[str, Any]]:
    """Translate a subset of Pinecone filter syntax to SQL WHERE over chunk_metadata.

    Supported operators:
    - {"field": {"$eq": value}}
    - {"field": {"$in": [v1, v2]}}
    - {"$and": [<filter>, ...]}
    - {"$or": [<filter>, ...]}

    Notes:
    - This is intentionally minimal for Phase 2.
    - Conditions are evaluated against `chunk_metadata` JSON/JSONB.
    """

    if not where:
        return "", {}

    params: Dict[str, Any] = {}
    counter = 0

    def next_param(val: Any) -> str:
        nonlocal counter
        counter += 1
        k = f"p{counter}"
        params[k] = val
        return k

    def render(obj: Any) -> str:
        if not isinstance(obj, dict):
            raise ValueError("Filter must be a dict")

        if "$and" in obj:
            parts = [render(x) for x in obj["$and"]]
            return "(" + " AND ".join(parts) + ")"
        if "$or" in obj:
            parts = [render(x) for x in obj["$or"]]
            return "(" + " OR ".join(parts) + ")"

        parts2: List[str] = []
        for field, cond in obj.items():
            if field in ("$and", "$or"):
                continue
            if isinstance(cond, dict) and "$eq" in cond:
                p = next_param(cond["$eq"])
                parts2.append(f"(chunk_metadata::jsonb ->> '{field}') = :{p}")
            elif isinstance(cond, dict) and "$in" in cond:
                # If metadata is an array, check membership via jsonb_array_elements_text;
                # otherwise compare scalar value via ANY(array).
                vals = list(cond["$in"] or [])
                p = next_param(vals)
                parts2.append(
                    "("
                    f"EXISTS ("
                    f"  SELECT 1 "
                    f"  FROM jsonb_array_elements_text(chunk_metadata::jsonb -> '{field}') AS elem(val) "
                    f"  WHERE elem.val = ANY(:{p})"
                    f") "
                    f"OR (chunk_metadata::jsonb ->> '{field}') = ANY(:{p})"
                    ")"
                )
            else:
                # treat scalar value as equality
                p = next_param(cond)
                parts2.append(f"(chunk_metadata::jsonb ->> '{field}') = :{p}")
        return "(" + " AND ".join(parts2) + ")"

    clause = render(where)
    return clause, params


class PgVectorStore:
    """Postgres/pgvector-backed vector store using rag_chunks.embedding."""

    def __init__(self, session_factory: Any | None = None) -> None:
        if session_factory is None:
            from amprenta_rag.database.session import db_session as session_factory  # type: ignore

        self._db_session = session_factory

    def query(
        self,
        vector: Sequence[float],
        top_k: int = 10,
        filter: Optional[PineconeFilter] = None,  # noqa: A002
        namespace: Optional[str] = None,
    ) -> List[VectorMatch]:
        where_sql, params = _filter_to_sql(filter)

        sql = """
        SELECT
          chunk_id,
          chunk_metadata,
          1 - (embedding <=> :vector) AS score
        FROM rag_chunks
        WHERE embedding IS NOT NULL
        """
        if namespace:
            params["namespace"] = namespace
            sql += " AND (chunk_metadata::jsonb ->> 'namespace') = :namespace"
        if where_sql:
            sql += f" AND {where_sql}"
        sql += " ORDER BY embedding <=> :vector LIMIT :top_k"

        params2: Dict[str, Any] = {"vector": list(vector), "top_k": top_k, **params}

        with self._db_session() as db:
            rows = db.execute(sa.text(sql), params2).mappings().all()

        out: List[VectorMatch] = []
        for r in rows:
            out.append(
                {
                    "id": str(r["chunk_id"]),
                    "score": float(r["score"]),
                    "metadata": cast(Dict[str, Any], r.get("chunk_metadata") or {}),
                }
            )
        return out

    def upsert(self, vectors: List[VectorRecord], namespace: Optional[str] = None) -> None:
        if not vectors:
            return

        with self._db_session() as db:
            for v in vectors:
                chunk_id = v.get("id")
                values = v.get("values")
                metadata = cast(Dict[str, Any], v.get("metadata") or {})
                if namespace:
                    metadata = {**metadata, "namespace": namespace}
                if not chunk_id or values is None:
                    continue

                embedding_model = metadata.get("embedding_model")
                # Update existing rag_chunks row (chunk_text is non-null; we assume rows exist).
                db.execute(
                    sa.text(
                        """
                        UPDATE rag_chunks
                        SET embedding = :embedding,
                            embedding_model = :embedding_model,
                            chunk_metadata = (COALESCE(chunk_metadata::jsonb, '{}'::jsonb) || :chunk_metadata::jsonb)
                        WHERE chunk_id = :chunk_id
                        """
                    ),
                    {
                        "embedding": list(values),
                        "embedding_model": embedding_model,
                        "chunk_metadata": metadata,
                        "chunk_id": chunk_id,
                    },
                )
            db.commit()

    def delete(
        self,
        filter: Optional[PineconeFilter] = None,  # noqa: A002
        ids: Optional[List[str]] = None,
        namespace: Optional[str] = None,
    ) -> None:
        # For parity with Pinecone delete, we clear embeddings rather than deleting rows.
        where_parts: List[str] = []
        params: Dict[str, Any] = {}

        if ids:
            where_parts.append("chunk_id = ANY(:ids)")
            params["ids"] = ids

        if namespace:
            where_parts.append("(chunk_metadata::jsonb ->> 'namespace') = :namespace")
            params["namespace"] = namespace

        where_sql, filter_params = _filter_to_sql(filter)
        if where_sql:
            where_parts.append(where_sql)
            params.update(filter_params)

        if not where_parts:
            return

        sql = f"""
        UPDATE rag_chunks
        SET embedding = NULL,
            embedding_model = NULL
        WHERE {' AND '.join(where_parts)}
        """

        with self._db_session() as db:
            db.execute(sa.text(sql), params)
            db.commit()


def get_vector_store() -> VectorStore:
    """Factory: return configured vector store backend."""
    cfg = get_config()
    backend = (cfg.vector_backend or "pgvector").lower()

    if backend == "pinecone":
        raise RuntimeError(
            "Pinecone backend is deprecated and has been removed. "
            "Please use pgvector backend instead. "
            "Set VECTOR_BACKEND=pgvector in your environment."
        )
    if backend in {"pgvector", "postgres", "postgresql"}:
        return PgVectorStore()

    raise ValueError(f"Unknown vector backend: {cfg.vector_backend!r}")


