## Overview

This document describes how to migrate Amprenta RAG vector search from **Pinecone** to **Postgres + pgvector**.

### Why migrate
- **Cost control**: vector storage + search runs in Postgres (no separate managed vector DB bill)
- **Operational simplicity**: fewer external services and credentials
- **Lower latency potential**: colocate vector search near your Postgres/RDS deployment
- **Better observability**: CloudWatch/RDS metrics cover both data + vector search

---

## Prerequisites

- **PostgreSQL 15+** (recommended; align with existing RDS versioning)
- **pgvector extension available** in your database instance
  - On RDS, ensure the engine supports pgvector (extension `vector`)
- **Database migrations configured** (Alembic)
- Pinecone credentials available **only** if you are backfilling existing embeddings

---

## Migration Steps

### 1) Run Alembic migrations

1. Ensure your database is reachable via your configured `DATABASE_URL` / Postgres env vars.
2. Run:

```bash
alembic upgrade head
```

This should:
- enable the `vector` extension
- add `rag_chunks.embedding` and `rag_chunks.embedding_model`
- create the HNSW index (if your migration includes it)

### 2) Backfill embeddings from Pinecone into Postgres (optional but typical)

Use:

```bash
python scripts/migrate_pinecone_to_pgvector.py --batch-size 1000 --workers 5 --resume
```

Notes:
- The script creates/uses a progress table: `migration_progress`
- Use `--resume` to skip chunk IDs already marked `done`
- Use `--dry-run` to verify the script can start without making any calls

### 3) Switch query backend

Set:
- `VECTOR_BACKEND=pgvector`

Restart the API and dashboard processes after updating environment variables.

---

## Configuration

### `VECTOR_BACKEND`

Controls which vector backend is used by `amprenta_rag.clients.vector_store.get_vector_store()`:
- `VECTOR_BACKEND=pinecone` (default)
- `VECTOR_BACKEND=pgvector`

### Optional related environment variables

- `PINECONE_API_KEY`, `PINECONE_INDEX_NAME`, `PINECONE_NAMESPACE` (for Pinecone)
- `DATABASE_URL` (or Postgres component env vars used by `amprenta_rag.config`)

---

## Rollback Procedures

### Rollback Phase A (Query rollback)

If you switched to pgvector and need to revert quickly:
- set `VECTOR_BACKEND=pinecone`
- restart services

### Rollback Phase B (Stop backfill)

If `scripts/migrate_pinecone_to_pgvector.py` is failing mid-run:
- fix the underlying error
- re-run with `--resume`

### Rollback Phase C (Schema rollback)

If you must roll back schema changes:
- run Alembic downgrade to the pre-pgvector revision
- **warning**: you will lose `rag_chunks.embedding` contents unless you preserve them

---

## Monitoring

After switching to pgvector, watch:
- **RDS CPU and memory** (vector search can be CPU-intensive)
- **DB connections** (pooling / task count)
- **Slow queries**: `embedding <=> :vector` queries and HNSW index usage
- **API latency**: end-to-end response times for RAG endpoints
- **Error rates**: failed vector queries, timeouts, Postgres connection errors

Consider adding:
- a canary query endpoint or scheduled health check that runs a small vector search

---

## Troubleshooting

### `psycopg2.errors.FeatureNotSupported: extension "vector" is not available`

- The pgvector extension is not installed/available in the DB engine.
- Fix by enabling pgvector support in your Postgres distribution (or use an RDS version that supports it).

### `embedding is Text` / vectors not comparable

- Your environment may not have the `pgvector` Python package installed.
- Install `pgvector` and ensure the Alembic migration has run.

### Backfill script is slow / rate-limited

- Reduce `--workers`
- Reduce `--batch-size`
- Re-run with `--resume`

### Recall drops after switching

- Confirm embeddings were successfully backfilled
- Verify you are querying the same namespace/metadata constraints
- Run the benchmark script to compare overlap and latency:

```bash
python scripts/benchmark_vector_stores.py --queries 100 --top-k 10 --use-db-vectors
```


