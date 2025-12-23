#!/usr/bin/env python3
"""Migrate Pinecone vectors into Postgres/pgvector.

This script backfills `rag_chunks.embedding` and `rag_chunks.embedding_model`
from an existing Pinecone index.

Features:
- **Resumable** via a `migration_progress` table (created if missing)
- **Batch commits** (default: 1000 rows/transaction)
- **Parallel workers** for Pinecone fetch + DB update
- **Dry-run** mode that performs no network/DB calls

Example:
  python scripts/migrate_pinecone_to_pgvector.py --batch-size 1000 --workers 5
  python scripts/migrate_pinecone_to_pgvector.py --dry-run
  python scripts/migrate_pinecone_to_pgvector.py --resume
"""

from __future__ import annotations

import argparse
import logging
import os
import random
import time
from concurrent.futures import ThreadPoolExecutor, as_completed
from dataclasses import dataclass
from typing import Any, Dict, Iterable, List, Optional, Sequence, Tuple


LOGGER = logging.getLogger("migrate_pinecone_to_pgvector")


@dataclass
class BatchResult:
    total_ids: int
    skipped_done: int
    updated: int
    missing_chunk: int
    errors: int


def _vector_literal(values: Sequence[float]) -> str:
    """Convert a python list of floats into pgvector literal string: `[1,2,3]`."""
    # Keep it compact (pgvector accepts standard float text)
    return "[" + ",".join(f"{float(x):.8g}" for x in values) + "]"


def _sleep_backoff(attempt: int, base: float = 0.5, cap: float = 20.0) -> None:
    # Exponential backoff with jitter.
    delay = min(cap, base * (2**attempt))
    delay = delay * (0.7 + random.random() * 0.6)
    time.sleep(delay)


def _pinecone_fetch_with_retry(index: Any, ids: List[str], namespace: str, max_retries: int = 6) -> Dict[str, Any]:
    """Fetch vectors from Pinecone with basic retry/backoff."""
    last_err: Exception | None = None
    for attempt in range(max_retries + 1):
        try:
            return index.fetch(ids=ids, namespace=namespace)
        except Exception as e:  # noqa: BLE001 - depends on Pinecone client
            last_err = e
            msg = str(e).lower()
            # crude rate-limit / transient detection
            transient = any(k in msg for k in ("429", "rate", "timeout", "tempor", "throttle", "503", "502"))
            if attempt < max_retries and transient:
                _sleep_backoff(attempt)
                continue
            raise
    raise RuntimeError(f"Pinecone fetch failed after retries: {last_err}")


def _parse_list_response(resp: Any) -> Tuple[List[str], Optional[str]]:
    """Best-effort parsing of Pinecone Index.list responses across client versions."""
    # Common shapes:
    # - {"vectors": [{"id": "..."}, ...], "pagination": {"next": "..."}} (dict)
    # - {"ids": ["..."], "pagination_token": "..."} (dict)
    # - object with .vectors and .pagination / .pagination_token
    data: Any = resp
    if not isinstance(resp, dict) and hasattr(resp, "__dict__"):
        # pinecone-sdk response objects are often attrs-like
        try:
            data = resp.__dict__
        except Exception:
            data = resp

    ids: List[str] = []
    next_token: Optional[str] = None

    if isinstance(data, dict):
        raw_vectors = data.get("vectors")
        raw_ids = data.get("ids")
        if raw_ids and isinstance(raw_ids, list):
            ids = [str(x) for x in raw_ids]
        elif raw_vectors and isinstance(raw_vectors, list):
            if raw_vectors and isinstance(raw_vectors[0], dict) and "id" in raw_vectors[0]:
                ids = [str(v["id"]) for v in raw_vectors]
            else:
                ids = [str(v) for v in raw_vectors]

        # pagination
        if "pagination_token" in data:
            next_token = data.get("pagination_token")
        elif "next_page_token" in data:
            next_token = data.get("next_page_token")
        elif "pagination" in data and isinstance(data["pagination"], dict):
            next_token = data["pagination"].get("next") or data["pagination"].get("next_token")
    else:
        # fallback: try attributes
        raw_vectors = getattr(resp, "vectors", None)
        raw_ids = getattr(resp, "ids", None)
        if raw_ids:
            ids = [str(x) for x in raw_ids]
        elif raw_vectors:
            ids = [str(getattr(v, "id", v)) for v in raw_vectors]
        pag = getattr(resp, "pagination", None)
        if pag is not None:
            next_token = getattr(pag, "next", None) or getattr(pag, "next_token", None)
        next_token = next_token or getattr(resp, "pagination_token", None) or getattr(resp, "next_page_token", None)

    if not ids:
        raise RuntimeError(f"Unable to parse Pinecone list response: {type(resp)}")

    return ids, (str(next_token) if next_token else None)


def iter_pinecone_ids(index: Any, *, namespace: str, page_size: int) -> Iterable[List[str]]:
    """Yield ID pages from Pinecone index."""
    token: Optional[str] = None

    # Prefer list() API if available.
    if not hasattr(index, "list"):
        raise RuntimeError("Pinecone Index has no .list() method; cannot paginate IDs safely in this environment.")

    while True:
        kwargs: Dict[str, Any] = {"namespace": namespace, "limit": page_size}
        # token param name varies; try common ones.
        if token:
            kwargs["pagination_token"] = token
        resp = index.list(**kwargs)
        ids, next_token = _parse_list_response(resp)
        yield ids
        if not next_token:
            break
        token = next_token


def _ensure_progress_table(db: Any) -> None:
    import sqlalchemy as sa

    db.execute(
        sa.text(
            """
            CREATE TABLE IF NOT EXISTS migration_progress (
              chunk_id TEXT PRIMARY KEY,
              status TEXT NOT NULL,
              migrated_at TIMESTAMPTZ NULL,
              error TEXT NULL
            )
            """
        )
    )


def _filter_done_ids(db: Any, ids: List[str]) -> Tuple[List[str], int]:
    """Return (remaining_ids, skipped_count)."""
    import sqlalchemy as sa

    if not ids:
        return [], 0
    rows = db.execute(
        sa.text("SELECT chunk_id FROM migration_progress WHERE status = 'done' AND chunk_id = ANY(:ids)"),
        {"ids": ids},
    ).fetchall()
    done = {r[0] for r in rows}
    remaining = [i for i in ids if i not in done]
    return remaining, len(done)


def _upsert_progress(db: Any, chunk_id: str, status: str, error: Optional[str] = None) -> None:
    import sqlalchemy as sa

    db.execute(
        sa.text(
            """
            INSERT INTO migration_progress (chunk_id, status, migrated_at, error)
            VALUES (:chunk_id, :status, CASE WHEN :status = 'done' THEN NOW() ELSE NULL END, :error)
            ON CONFLICT (chunk_id) DO UPDATE
              SET status = EXCLUDED.status,
                  migrated_at = EXCLUDED.migrated_at,
                  error = EXCLUDED.error
            """
        ),
        {"chunk_id": chunk_id, "status": status, "error": error},
    )


def _update_chunk_embedding(db: Any, *, chunk_id: str, embedding: Sequence[float], embedding_model: str) -> bool:
    """Update rag_chunks row; returns True if updated, False if chunk_id not found."""
    import sqlalchemy as sa

    emb = _vector_literal(embedding)
    res = db.execute(
        sa.text(
            """
            UPDATE rag_chunks
            SET embedding = :embedding::vector,
                embedding_model = :embedding_model,
                updated_at = NOW()
            WHERE chunk_id = :chunk_id
            """
        ),
        {"chunk_id": chunk_id, "embedding": emb, "embedding_model": embedding_model},
    )
    return bool(getattr(res, "rowcount", 0))


def _process_id_batch(
    *,
    ids: List[str],
    namespace: str,
    resume: bool,
    embedding_model: str,
    dry_run: bool,
) -> BatchResult:
    """Worker function: optionally skip done, fetch vectors, write to Postgres, write progress."""
    if dry_run:
        return BatchResult(total_ids=len(ids), skipped_done=0, updated=0, missing_chunk=0, errors=0)

    # Lazy imports: keep --dry-run fast and dependency-light.
    from amprenta_rag.clients.pinecone_client import get_pinecone_index
    from amprenta_rag.database.session import db_session

    index = get_pinecone_index()

    # 1) Filter done IDs if resume enabled
    skipped = 0
    if resume:
        with db_session() as db:
            _ensure_progress_table(db)
            ids, skipped = _filter_done_ids(db, ids)
    if not ids:
        return BatchResult(total_ids=skipped, skipped_done=skipped, updated=0, missing_chunk=0, errors=0)

    # 2) Fetch vectors from Pinecone
    fetched = _pinecone_fetch_with_retry(index, ids, namespace=namespace)
    vectors_map: Dict[str, Any] = {}
    if isinstance(fetched, dict):
        # pinecone returns {"vectors": {id: {...}}, "namespace": "..."}
        vectors_map = fetched.get("vectors") or {}
    else:
        vectors_map = getattr(fetched, "vectors", {}) or {}

    updated = 0
    missing = 0
    errors = 0

    # 3) Update Postgres in one transaction for the batch
    with db_session() as db:
        _ensure_progress_table(db)

        for chunk_id in ids:
            try:
                v = vectors_map.get(chunk_id)
                if not v:
                    _upsert_progress(db, chunk_id, "error", error="missing_in_pinecone_fetch")
                    errors += 1
                    continue

                # v shape is usually {"id":..., "values":[...], "metadata":{...}}
                values = v.get("values") if isinstance(v, dict) else getattr(v, "values", None)
                if not values:
                    _upsert_progress(db, chunk_id, "error", error="missing_values")
                    errors += 1
                    continue

                ok = _update_chunk_embedding(db, chunk_id=chunk_id, embedding=values, embedding_model=embedding_model)
                if ok:
                    _upsert_progress(db, chunk_id, "done", error=None)
                    updated += 1
                else:
                    _upsert_progress(db, chunk_id, "missing", error="chunk_id_not_found_in_rag_chunks")
                    missing += 1
            except Exception as e:  # noqa: BLE001 - isolate batch failures
                _upsert_progress(db, chunk_id, "error", error=str(e)[:4000])
                errors += 1

    return BatchResult(total_ids=len(ids) + skipped, skipped_done=skipped, updated=updated, missing_chunk=missing, errors=errors)


def main(argv: Optional[List[str]] = None) -> int:
    parser = argparse.ArgumentParser(description="Migrate Pinecone vectors to Postgres/pgvector.")
    parser.add_argument("--batch-size", type=int, default=1000, help="IDs per batch/transaction (default: 1000)")
    parser.add_argument("--workers", type=int, default=1, help="Parallel worker threads (default: 1)")
    parser.add_argument("--dry-run", action="store_true", help="Do not connect to Pinecone/Postgres; only print plan.")
    parser.add_argument("--resume", action="store_true", help="Skip vectors already marked done in migration_progress.")
    parser.add_argument("--namespace", type=str, default=None, help="Override Pinecone namespace (defaults to config).")
    parser.add_argument("--log-level", type=str, default="INFO", help="Logging level (INFO, DEBUG, ...)")
    args = parser.parse_args(argv)

    logging.basicConfig(
        level=getattr(logging, args.log_level.upper(), logging.INFO),
        format="%(asctime)s %(levelname)s %(name)s - %(message)s",
    )

    if args.batch_size <= 0:
        LOGGER.error("--batch-size must be > 0")
        return 2
    if args.workers <= 0:
        LOGGER.error("--workers must be > 0")
        return 2

    if args.dry_run:
        # Keep this path free of config/network/db requirements.
        ns = args.namespace or os.getenv("PINECONE_NAMESPACE", "<from config>")
        LOGGER.info("Dry-run enabled. No network/DB calls will be made.")
        LOGGER.info("Would migrate Pinecone vectors from namespace=%s into rag_chunks.embedding", ns)
        LOGGER.info("Planned settings: batch_size=%d workers=%d resume=%s", args.batch_size, args.workers, args.resume)
        return 0

    # Lazy imports so --dry-run is safe even without deps/creds.
    from amprenta_rag.config import get_config
    from amprenta_rag.clients.pinecone_client import get_pinecone_index
    from amprenta_rag.clients.openai_client import get_default_models
    from amprenta_rag.database.session import db_session

    cfg = get_config()
    namespace = args.namespace or cfg.pinecone.namespace

    # Determine embedding model to store (best-effort).
    try:
        _, emb_model = get_default_models()
    except Exception:
        emb_model = "unknown"

    # Ensure progress table exists up-front.
    with db_session() as db:
        _ensure_progress_table(db)

    index = get_pinecone_index()

    # Estimate total for % logging.
    total_est: Optional[int] = None
    try:
        stats = index.describe_index_stats()
        if isinstance(stats, dict):
            ns_stats = (stats.get("namespaces") or {}).get(namespace) or {}
            total_est = ns_stats.get("vector_count")
        else:
            namespaces = getattr(stats, "namespaces", None) or {}
            ns_stats = namespaces.get(namespace) if isinstance(namespaces, dict) else None
            total_est = getattr(ns_stats, "vector_count", None) if ns_stats else None
    except Exception:
        total_est = None

    LOGGER.info("Starting migration: namespace=%s batch_size=%d workers=%d resume=%s", namespace, args.batch_size, args.workers, args.resume)
    if total_est is not None:
        LOGGER.info("Pinecone namespace vector_count≈%d", total_est)

    processed_ids = 0
    skipped_done = 0
    updated = 0
    missing = 0
    errors = 0

    # Produce ID pages (page size == batch size)
    id_pages = iter_pinecone_ids(index, namespace=namespace, page_size=args.batch_size)

    with ThreadPoolExecutor(max_workers=args.workers) as ex:
        futures = []
        for ids in id_pages:
            futures.append(
                ex.submit(
                    _process_id_batch,
                    ids=ids,
                    namespace=namespace,
                    resume=args.resume,
                    embedding_model=emb_model,
                    dry_run=False,
                )
            )

        for fut in as_completed(futures):
            r = fut.result()
            processed_ids += r.total_ids
            skipped_done += r.skipped_done
            updated += r.updated
            missing += r.missing_chunk
            errors += r.errors

            if total_est:
                pct = min(100.0, (processed_ids / max(1, total_est)) * 100.0)
                LOGGER.info(
                    "Progress: %.1f%% (%d/%d) updated=%d skipped=%d missing=%d errors=%d",
                    pct,
                    processed_ids,
                    total_est,
                    updated,
                    skipped_done,
                    missing,
                    errors,
                )
            else:
                LOGGER.info(
                    "Progress: processed=%d updated=%d skipped=%d missing=%d errors=%d",
                    processed_ids,
                    updated,
                    skipped_done,
                    missing,
                    errors,
                )

    LOGGER.info("Done. processed=%d updated=%d skipped=%d missing=%d errors=%d", processed_ids, updated, skipped_done, missing, errors)
    return 0 if errors == 0 else 1


if __name__ == "__main__":
    raise SystemExit(main())


