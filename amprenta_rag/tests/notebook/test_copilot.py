from __future__ import annotations

import json

import pytest

from amprenta_rag.notebook.context import AnalysisContext, generate_context_cell
import amprenta_rag.notebook.copilot as copilot


def test_analysis_context_to_dict_from_dict_roundtrip():
    ctx = AnalysisContext(
        entity_type="dataset",
        entity_id="00000000-0000-0000-0000-000000000001",
        campaign_id="camp-1",
        plate_id="plate-1",
        compound_id="cmpd-1",
        version=2,
        timestamp="2020-01-01T00:00:00",
        metadata={"k": "v"},
    )

    d = ctx.to_dict()
    assert d["entity_type"] == "dataset"
    assert d["entity_id"] == "00000000-0000-0000-0000-000000000001"
    assert d["campaign_id"] == "camp-1"
    assert d["plate_id"] == "plate-1"
    assert d["compound_id"] == "cmpd-1"
    assert d["version"] == 2
    assert d["timestamp"] == "2020-01-01T00:00:00"
    assert d["metadata"] == {"k": "v"}

    ctx2 = AnalysisContext.from_dict(d)
    assert ctx2.to_dict() == d


def test_analysis_context_to_json_from_json_roundtrip():
    ctx = AnalysisContext(
        entity_type="experiment",
        entity_id="00000000-0000-0000-0000-000000000002",
        version=1,
        timestamp="2021-02-03T04:05:06",
        metadata={"a": 1, "b": True},
    )

    s = ctx.to_json()
    # sanity: it's valid json
    assert json.loads(s)["entity_id"] == "00000000-0000-0000-0000-000000000002"

    ctx2 = AnalysisContext.from_json(s)
    assert ctx2.to_dict() == ctx.to_dict()


def test_generate_context_cell_includes_expected_scaffold():
    ctx = AnalysisContext(
        entity_type="dataset",
        entity_id="00000000-0000-0000-0000-000000000003",
        version=1,
        timestamp="2022-01-01T00:00:00",
        metadata={"x": "y"},
    )

    cell = generate_context_cell(ctx)
    assert "from amprenta_rag.notebook.context import AnalysisContext" in cell
    assert "CONTEXT = AnalysisContext.from_dict(" in cell
    assert 'print(f"Analysis: {CONTEXT.entity_type} {CONTEXT.entity_id}")' in cell
    # Ensure the context dict is embedded.
    assert repr(ctx.to_dict()) in cell


def test_synthesize_cell_calls_llm_and_returns_code(monkeypatch: pytest.MonkeyPatch):
    calls = []

    def fake_call_model(model_name, messages, temperature=0.2):
        calls.append((model_name, messages, temperature))
        return "print('hello')"

    monkeypatch.setattr(copilot, "call_model", fake_call_model)
    monkeypatch.setattr(copilot, "DEFAULT_NOTEBOOK_MODEL", "test-model")

    ctx = AnalysisContext(entity_type="dataset", entity_id="ds-1", timestamp="2020-01-01T00:00:00")
    out = copilot.synthesize_cell("show dataset summary", ctx)

    assert out == "print('hello')"
    assert len(calls) == 1
    model_name, messages, temperature = calls[0]
    assert model_name == "test-model"
    assert temperature == 0.2
    assert messages[0]["role"] == "system"
    assert "Generate ONE Python notebook cell" in messages[0]["content"]
    assert messages[1]["role"] == "user"
    assert "Intent: show dataset summary" in messages[1]["content"]
    assert ctx.to_json() in messages[1]["content"]


def test_fix_cell_calls_llm_and_returns_code(monkeypatch: pytest.MonkeyPatch):
    calls = []

    def fake_call_model(model_name, messages, temperature=0.2):
        calls.append((model_name, messages, temperature))
        return "print('fixed')"

    monkeypatch.setattr(copilot, "call_model", fake_call_model)
    monkeypatch.setattr(copilot, "DEFAULT_NOTEBOOK_MODEL", "test-model")

    out = copilot.fix_cell(code="print(1/0)", error="ZeroDivisionError")
    assert out == "print('fixed')"
    assert len(calls) == 1
    model_name, messages, temperature = calls[0]
    assert model_name == "test-model"
    assert temperature == 0.1
    assert "The following cell failed:" in messages[1]["content"]
    assert "ZeroDivisionError" in messages[1]["content"]


def test_explain_cell_calls_llm_and_returns_text(monkeypatch: pytest.MonkeyPatch):
    calls = []

    def fake_call_model(model_name, messages, temperature=0.2):
        calls.append((model_name, messages, temperature))
        return "It prints a value."

    monkeypatch.setattr(copilot, "call_model", fake_call_model)
    monkeypatch.setattr(copilot, "DEFAULT_NOTEBOOK_MODEL", "test-model")

    out = copilot.explain_cell(code="print('x')")
    assert out == "It prints a value."
    assert len(calls) == 1
    model_name, messages, temperature = calls[0]
    assert model_name == "test-model"
    assert temperature == 0.2
    assert "Explain what this cell does:" in messages[1]["content"]


def test_summarize_notebook_valid_input(monkeypatch: pytest.MonkeyPatch):
    calls = []

    def fake_call_model(model_name, messages, temperature=0.2):
        calls.append((model_name, messages, temperature))
        return json.dumps(
            {
                "title": "Transcriptomics Analysis",
                "entity_summary": "Dataset ds-1 (transcriptomics).",
                "methods": "Load dataset, basic QC, differential expression.",
                "key_findings": "Top pathways suggest immune activation.",
            }
        )

    monkeypatch.setattr(copilot, "call_model", fake_call_model)
    monkeypatch.setattr(copilot, "DEFAULT_NOTEBOOK_MODEL", "test-model")

    nb = {
        "cells": [
            {"cell_type": "markdown", "source": "# Title\nSome intro"},
            {"cell_type": "code", "source": "print('hi')"},
        ]
    }
    out = copilot.summarize_notebook(nb)

    assert out["title"] == "Transcriptomics Analysis"
    assert out["entity_summary"]
    assert out["methods"]
    assert out["key_findings"]
    assert out["cell_count"] == 2

    assert len(calls) == 1
    model_name, messages, temperature = calls[0]
    assert model_name == "test-model"
    assert temperature == 0.3
    assert messages[0]["role"] == "system"
    assert "Return ONLY a JSON object" in messages[0]["content"]
    assert messages[1]["role"] == "user"
    assert "Notebook cell count: 2" in messages[1]["content"]


def test_summarize_notebook_empty_notebook(monkeypatch: pytest.MonkeyPatch):
    # Should return fallback without calling the LLM.
    called = {"n": 0}

    def fake_call_model(model_name, messages, temperature=0.2):
        called["n"] += 1
        return "{}"

    monkeypatch.setattr(copilot, "call_model", fake_call_model)
    monkeypatch.setattr(copilot, "DEFAULT_NOTEBOOK_MODEL", "test-model")

    out = copilot.summarize_notebook({"cells": []})
    assert out["title"] == "Notebook Summary"
    assert out["cell_count"] == 0
    assert out["entity_summary"] == ""
    assert out["methods"] == ""
    assert out["key_findings"] == ""
    assert called["n"] == 0


def test_summarize_notebook_llm_failure(monkeypatch: pytest.MonkeyPatch):
    def fake_call_model(model_name, messages, temperature=0.2):
        return "NOT JSON"

    monkeypatch.setattr(copilot, "call_model", fake_call_model)
    monkeypatch.setattr(copilot, "DEFAULT_NOTEBOOK_MODEL", "test-model")

    nb = {"cells": [{"cell_type": "code", "source": "x = 1"}]}
    out = copilot.summarize_notebook(nb)
    assert out["title"] == "Notebook Summary"
    assert out["cell_count"] == 1
    assert out["key_findings"] == "NOT JSON"

