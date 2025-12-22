from __future__ import annotations

import json

import nbformat
import pytest

import amprenta_rag.notebook.entity_extractor as entity_extractor
import amprenta_rag.notebook.notebook_planner as notebook_planner
import amprenta_rag.notebook.query_to_notebook as query_to_notebook
from amprenta_rag.notebook.context import AnalysisContext


def test_entity_extraction_resolves_and_returns_context(monkeypatch: pytest.MonkeyPatch):
    # Fake LLM output with names
    def fake_llm(query: str):
        return {
            "datasets": ["ALS dataset 1"],
            "experiments": [],
            "compounds": [],
            "campaigns": [],
            "primary": {"type": "dataset", "value": "ALS dataset 1"},
        }

    # Fake DB resolution
    def fake_resolve_by_name(db, model, names):
        if not names:
            return []
        assert names == ["ALS dataset 1"]
        return ["00000000-0000-0000-0000-000000000001"]

    class DummyDB:
        pass

    class DummySession:
        def __enter__(self):
            return DummyDB()

        def __exit__(self, exc_type, exc, tb):
            return False

    monkeypatch.setattr(entity_extractor, "_llm_extract_candidates", fake_llm)
    monkeypatch.setattr(entity_extractor, "_resolve_by_name", fake_resolve_by_name)
    monkeypatch.setattr(entity_extractor, "db_session", lambda: DummySession())

    out = entity_extractor.extract_entities("What does ALS dataset 1 show?")
    assert out["dataset_ids"] == ["00000000-0000-0000-0000-000000000001"]
    assert out["experiment_ids"] == []
    assert "context" in out
    ctx = AnalysisContext.from_dict(out["context"])
    assert ctx.entity_type == "dataset"
    assert ctx.entity_id == "00000000-0000-0000-0000-000000000001"


def test_notebook_planner_returns_valid_plan(monkeypatch: pytest.MonkeyPatch):
    def fake_call_model(model_name, messages, temperature=0.2):
        # Return JSON list
        return json.dumps(
            [
                {"type": "markdown", "section": "Context", "intent": "Describe the question."},
                {"type": "code", "section": "Data Loading", "intent": "Load the dataset and show head."},
            ]
        )

    monkeypatch.setattr(notebook_planner, "call_model", fake_call_model)
    monkeypatch.setattr(notebook_planner, "DEFAULT_PLANNER_MODEL", "test-model")

    ctx = AnalysisContext(entity_type="dataset", entity_id="ds-1", timestamp="2020-01-01T00:00:00")
    plan = notebook_planner.plan_notebook("q", ctx, rag_chunks=["c1", "c2"])
    assert isinstance(plan, list)
    assert plan[0]["type"] == "markdown"
    assert plan[1]["type"] == "code"
    assert plan[1]["section"] == "Data Loading"


def test_generate_notebook_orchestrates_with_mocks(monkeypatch: pytest.MonkeyPatch):
    # Mock entity extraction
    monkeypatch.setattr(
        query_to_notebook,
        "extract_entities",
        lambda q: {
            "dataset_ids": ["ds-1"],
            "experiment_ids": [],
            "compound_ids": [],
            "campaign_ids": [],
            "context": AnalysisContext(entity_type="dataset", entity_id="ds-1", timestamp="2020-01-01T00:00:00").to_dict(),
        },
    )

    # Mock RAG retrieval
    class DummyRAG:
        context_chunks = ["evidence 1", "evidence 2"]

    monkeypatch.setattr(query_to_notebook, "query_rag", lambda *a, **k: DummyRAG())

    # Mock planner
    monkeypatch.setattr(
        query_to_notebook,
        "plan_notebook",
        lambda query, context, rag_chunks: [
            {"type": "markdown", "section": "Context", "intent": "Explain context."},
            {"type": "code", "section": "Analysis", "intent": "Compute summary stats."},
        ],
    )

    # Mock code synthesis
    monkeypatch.setattr(query_to_notebook, "synthesize_cell", lambda intent, context: "print('ok')")

    nb = query_to_notebook.generate_notebook("test query")
    assert nb["nbformat"] == 4
    assert len(nb["cells"]) >= 3
    assert nb["cells"][0]["cell_type"] == "code"  # context cell

    # Ensure it's a valid nbformat structure
    nbformat.validate(nb)


