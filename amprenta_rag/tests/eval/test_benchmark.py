from __future__ import annotations

import json
import os
from pathlib import Path

from amprenta_rag.eval import benchmark as bm


def test_load_tasks_from_json(tmp_path: Path):
    data = [{"id": "t1", "role": "rag_query", "input": {"q": "hi"}, "expected_contains": "hello"}]
    json_path = tmp_path / "benchmark_tasks.json"
    json_path.write_text(json.dumps(data))

    tasks = bm._load_tasks_from_json(str(json_path))
    assert tasks[0].id == "t1"
    assert tasks[0].input["q"] == "hi"


def test_maybe_find_default_json(tmp_path: Path):
    default = tmp_path / "RAG" / "experiments" / "evolution"
    default.mkdir(parents=True)
    target = default / "benchmark_tasks.json"
    target.write_text("[]")
    found = bm._maybe_find_default_json(str(tmp_path))
    assert found == str(target)


def test_scan_tests_for_benchmark_markers(tmp_path: Path):
    tests_dir = tmp_path / "tests"
    tests_dir.mkdir()
    fpath = tests_dir / "sample.py"
    marker = '{"id":"m1","role":"chat","input":{"q":"x"}}'
    fpath.write_text(f"# BENCHMARK: {marker}\n")
    tasks = bm._scan_tests_for_benchmark_markers(str(tests_dir))
    assert tasks and tasks[0].id == "m1"


def test_load_benchmark_tasks_prefers_default_json(tmp_path: Path):
    # Write default JSON
    default = tmp_path / "RAG" / "experiments" / "evolution"
    default.mkdir(parents=True)
    (default / "benchmark_tasks.json").write_text('[{"id":"t1","role":"r","input":{}}]')
    tasks = bm.load_benchmark_tasks(str(tmp_path))
    assert tasks and tasks[0].id == "t1"

