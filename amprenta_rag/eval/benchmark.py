from __future__ import annotations

import json
import os
from dataclasses import dataclass
from typing import Any, Dict, Iterable, List, Optional, Sequence


@dataclass(frozen=True)
class BenchmarkTask:
    """
    Minimal benchmark task schema used by the evaluator.
    """

    id: str
    role: str  # e.g., "rag_query", "chat_turn"
    input: Dict[str, Any]
    expected_contains: Optional[str] = None


def _load_tasks_from_json(json_path: str) -> List[BenchmarkTask]:
    with open(json_path, "r", encoding="utf-8") as f:
        raw = json.load(f)
    tasks: List[BenchmarkTask] = []
    for item in raw:
        tasks.append(
            BenchmarkTask(
                id=str(item.get("id") or item.get("name") or f"task_{len(tasks)}"),
                role=str(item.get("role", "rag_query")),
                input=dict(item.get("input", {})),
                expected_contains=item.get("expected_contains"),
            )
        )
    return tasks


def _maybe_find_default_json(base_dir: str) -> Optional[str]:
    candidates = [
        os.path.join(base_dir, "RAG", "experiments", "evolution", "benchmark_tasks.json"),
        os.path.join(base_dir, "experiments", "evolution", "benchmark_tasks.json"),
    ]
    for p in candidates:
        if os.path.isfile(p):
            return p
    return None


def _scan_tests_for_benchmark_markers(tests_dir: str) -> List[BenchmarkTask]:
    """
    Best-effort static scan for lines containing 'BENCHMARK:' markers.
    Each marker line should contain JSON after the marker.
    Example:
      # BENCHMARK: {"id":"ds_summary_1","role":"rag_query","input":{"question":"..."}}
    """
    tasks: List[BenchmarkTask] = []
    if not os.path.isdir(tests_dir):
        return tasks

    for root, _, files in os.walk(tests_dir):
        for name in files:
            if not name.endswith(".py"):
                continue
            full = os.path.join(root, name)
            try:
                with open(full, "r", encoding="utf-8") as f:
                    for line in f:
                        if "BENCHMARK:" in line:
                            try:
                                marker = line.split("BENCHMARK:", 1)[1].strip()
                                obj = json.loads(marker)
                                tasks.append(
                                    BenchmarkTask(
                                        id=str(
                                            obj.get("id")
                                            or obj.get("name")
                                            or f"{name}_{len(tasks)}"
                                        ),
                                        role=str(obj.get("role", "rag_query")),
                                        input=dict(obj.get("input", {})),
                                        expected_contains=obj.get("expected_contains"),
                                    )
                                )
                            except Exception:
                                # Skip malformed markers
                                continue
            except (OSError, UnicodeDecodeError):
                continue
    return tasks


def load_benchmark_tasks(base_dir: str) -> List[BenchmarkTask]:
    """
    Load benchmark tasks from a well-known JSON file or by scanning tests.

    Order of precedence:
      1) RAG/experiments/evolution/benchmark_tasks.json
      2) experiments/evolution/benchmark_tasks.json
      3) Scan amprenta_rag/tests for 'BENCHMARK:' markers
      4) Fallback: empty list
    """
    # 1/2) Well-known JSON locations
    default_json = _maybe_find_default_json(base_dir)
    if default_json:
        try:
            return _load_tasks_from_json(default_json)
        except Exception:
            pass

    # 3) Scan tests directory
    tests_dir_candidates: Sequence[str] = [
        os.path.join(base_dir, "RAG", "amprenta_rag", "tests"),
        os.path.join(base_dir, "amprenta_rag", "tests"),
        os.path.join(base_dir, "tests"),
    ]
    for tdir in tests_dir_candidates:
        tasks = _scan_tests_for_benchmark_markers(tdir)
        if tasks:
            return tasks

    # 4) Fallback: empty
    return []




