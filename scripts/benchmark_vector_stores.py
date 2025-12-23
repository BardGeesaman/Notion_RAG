#!/usr/bin/env python3
"""Benchmark Pinecone vs pgvector query latency and retrieval overlap.

This script is intended for the pgvector migration to compare:
- **Latency**: p50 / p95 / p99 per backend across N queries
- **Recall@k**: overlap of result IDs between backends

Query vectors are sourced from Postgres `rag_chunks.embedding` when available;
otherwise random vectors are generated (dimension defaults to 3072).

Usage:
  python scripts/benchmark_vector_stores.py --queries 100 --top-k 10
  python scripts/benchmark_vector_stores.py --queries 50 --top-k 20 --output results.json
"""

from __future__ import annotations

import argparse
import json
import logging
import random
import statistics
import time
from dataclasses import asdict, dataclass
from typing import Any, List, Optional, Sequence, Tuple

import sqlalchemy as sa


LOGGER = logging.getLogger("benchmark_vector_stores")


DEFAULT_DIM = 3072


@dataclass
class LatencyStats:
    count: int
    p50_ms: float
    p95_ms: float
    p99_ms: float
    mean_ms: float
    min_ms: float
    max_ms: float


@dataclass
class BackendReport:
    name: str
    latency: LatencyStats


@dataclass
class BenchmarkReport:
    queries: int
    top_k: int
    namespace: Optional[str]
    vector_dim: int
    recall_at_k: float
    pinecone: BackendReport
    pgvector: BackendReport


def _quantile(xs: List[float], q: float) -> float:
    """Compute quantile with simple nearest-rank over sorted list."""
    if not xs:
        return 0.0
    xs2 = sorted(xs)
    idx = int((len(xs2) - 1) * q)
    return float(xs2[idx])


def _latency_stats(durations_ms: List[float]) -> LatencyStats:
    if not durations_ms:
        return LatencyStats(count=0, p50_ms=0.0, p95_ms=0.0, p99_ms=0.0, mean_ms=0.0, min_ms=0.0, max_ms=0.0)
    return LatencyStats(
        count=len(durations_ms),
        p50_ms=_quantile(durations_ms, 0.50),
        p95_ms=_quantile(durations_ms, 0.95),
        p99_ms=_quantile(durations_ms, 0.99),
        mean_ms=float(statistics.mean(durations_ms)),
        min_ms=float(min(durations_ms)),
        max_ms=float(max(durations_ms)),
    )


def _print_table(report: BenchmarkReport) -> None:
    def row(name: str, ls: LatencyStats) -> str:
        return (
            f"{name:<10}  "
            f"{ls.count:>6}  "
            f"{ls.p50_ms:>9.2f}  "
            f"{ls.p95_ms:>9.2f}  "
            f"{ls.p99_ms:>9.2f}  "
            f"{ls.mean_ms:>9.2f}  "
            f"{ls.min_ms:>9.2f}  "
            f"{ls.max_ms:>9.2f}"
        )

    print("")
    print("Vector Store Benchmark")
    print(f"- queries: {report.queries}")
    print(f"- top_k:   {report.top_k}")
    print(f"- ns:      {report.namespace or '<default>'}")
    print(f"- dim:     {report.vector_dim}")
    print(f"- recall@{report.top_k}: {report.recall_at_k:.3f}")
    print("")
    print(f"{'backend':<10}  {'count':>6}  {'p50(ms)':>9}  {'p95(ms)':>9}  {'p99(ms)':>9}  {'mean(ms)':>9}  {'min(ms)':>9}  {'max(ms)':>9}")
    print("-" * 86)
    print(row("pinecone", report.pinecone.latency))
    print(row("pgvector", report.pgvector.latency))
    print("")


def _load_query_vectors_from_db(limit: int) -> Tuple[List[List[float]], int]:
    """Return (vectors, dim) from rag_chunks.embedding when possible."""
    from amprenta_rag.database.session import db_session

    # We fetch random embeddings from DB; assumes pgvector is installed and column exists.
    sql = sa.text(
        """
        SELECT embedding
        FROM rag_chunks
        WHERE embedding IS NOT NULL
        ORDER BY random()
        LIMIT :limit
        """
    )
    with db_session() as db:
        rows = db.execute(sql, {"limit": limit}).fetchall()

    vectors: List[List[float]] = []
    for (emb,) in rows:
        # emb might be pgvector.Vector, list[float], or string literal depending on driver.
        if emb is None:
            continue
        if isinstance(emb, (list, tuple)):
            v = [float(x) for x in emb]
        else:
            s = str(emb).strip()
            # handle formats like "[1, 2, 3]" or "1,2,3"
            s = s.strip("[]")
            v = [float(x) for x in s.split(",") if x.strip()]
        if v:
            vectors.append(v)

    dim = len(vectors[0]) if vectors else DEFAULT_DIM
    return vectors, dim


def _random_vectors(n: int, dim: int, seed: int) -> List[List[float]]:
    r = random.Random(seed)
    # Normal-ish distribution; exact distribution not critical for timing.
    return [[r.uniform(-1.0, 1.0) for _ in range(dim)] for _ in range(n)]


def _time_query(store: Any, vector: Sequence[float], top_k: int, namespace: Optional[str]) -> Tuple[float, List[str]]:
    t0 = time.perf_counter()
    matches = store.query(vector=vector, top_k=top_k, filter=None, namespace=namespace)
    dt_ms = (time.perf_counter() - t0) * 1000.0
    ids = [str(m.get("id", "")) for m in matches if isinstance(m, dict)]
    return dt_ms, ids


def main(argv: Optional[List[str]] = None) -> int:
    parser = argparse.ArgumentParser(description="Benchmark Pinecone vs pgvector VectorStore backends.")
    parser.add_argument("--queries", type=int, default=100, help="Number of queries to run (default: 100)")
    parser.add_argument("--top-k", type=int, default=10, help="Top-k to retrieve (default: 10)")
    parser.add_argument("--namespace", type=str, default=None, help="Namespace override (defaults to config)")
    parser.add_argument("--vector-dim", type=int, default=DEFAULT_DIM, help="Vector dimension if using random vectors")
    parser.add_argument("--seed", type=int, default=1337, help="Random seed")
    parser.add_argument("--output", type=str, default=None, help="Write JSON report to this path (optional)")
    parser.add_argument("--use-db-vectors", action="store_true", help="Use Postgres embeddings as query vectors when possible")
    parser.add_argument("--log-level", type=str, default="INFO", help="Logging level")
    args = parser.parse_args(argv)

    logging.basicConfig(
        level=getattr(logging, args.log_level.upper(), logging.INFO),
        format="%(asctime)s %(levelname)s %(name)s - %(message)s",
    )

    if args.queries <= 0 or args.top_k <= 0:
        LOGGER.error("--queries and --top-k must be > 0")
        return 2

    # Load namespace from config if not provided
    namespace = args.namespace
    if namespace is None:
        try:
            from amprenta_rag.config import get_config

            namespace = get_config().pinecone.namespace
        except Exception:
            namespace = None

    # Init stores explicitly (avoid global VECTOR_BACKEND selection)
    from amprenta_rag.clients.vector_store import PineconeStore, PgVectorStore

    pinecone_store = PineconeStore()
    pgvector_store = PgVectorStore()

    # Prepare query vectors
    vectors: List[List[float]] = []
    dim = args.vector_dim
    if args.use_db_vectors:
        try:
            vectors, dim = _load_query_vectors_from_db(args.queries)
            if len(vectors) < args.queries:
                LOGGER.warning("Only found %d embeddings in DB; generating %d random vectors to fill.",
                               len(vectors), args.queries - len(vectors))
                vectors.extend(_random_vectors(args.queries - len(vectors), dim, args.seed))
        except Exception as e:
            LOGGER.warning("Failed to load DB vectors (%s). Falling back to random vectors.", e)
            vectors = _random_vectors(args.queries, dim, args.seed)
    else:
        vectors = _random_vectors(args.queries, dim, args.seed)

    pine_ms: List[float] = []
    pg_ms: List[float] = []
    recall_sum = 0.0

    for i, v in enumerate(vectors, start=1):
        dt_p, ids_p = _time_query(pinecone_store, v, args.top_k, namespace)
        dt_g, ids_g = _time_query(pgvector_store, v, args.top_k, namespace)
        pine_ms.append(dt_p)
        pg_ms.append(dt_g)

        # Recall@k: treat Pinecone results as baseline.
        set_p = set(ids_p)
        set_g = set(ids_g)
        denom = max(1, len(set_p))
        recall_sum += len(set_p & set_g) / denom

        if i % max(1, args.queries // 10) == 0:
            LOGGER.info("Progress %d/%d", i, args.queries)

    report = BenchmarkReport(
        queries=args.queries,
        top_k=args.top_k,
        namespace=namespace,
        vector_dim=dim,
        recall_at_k=recall_sum / max(1, args.queries),
        pinecone=BackendReport(name="pinecone", latency=_latency_stats(pine_ms)),
        pgvector=BackendReport(name="pgvector", latency=_latency_stats(pg_ms)),
    )

    # Human-readable table
    _print_table(report)

    # JSON
    out_json = json.dumps(asdict(report), indent=2, sort_keys=True)
    if args.output:
        with open(args.output, "w", encoding="utf-8") as f:
            f.write(out_json + "\n")
        LOGGER.info("Wrote JSON report to %s", args.output)
    else:
        print(out_json)

    return 0


if __name__ == "__main__":
    raise SystemExit(main())


