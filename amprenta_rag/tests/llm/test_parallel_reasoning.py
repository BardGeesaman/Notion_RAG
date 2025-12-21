from __future__ import annotations

import sys
from types import ModuleType
import pytest

from amprenta_rag.llm import parallel_reasoning as pr


def test_run_parallel_models_success(monkeypatch):
    def fake_call_model(model_name, messages, temperature=0.2):
        return f"Response from {model_name}"

    monkeypatch.setattr(pr, "call_model", fake_call_model)

    results = pr.run_parallel_models("prompt", ["gpt-4", "gpt-3.5-turbo"])
    assert len(results) == 2
    assert results[0]["model"] == "gpt-4"
    assert results[0]["response"] == "Response from gpt-4"
    assert results[0]["success"] is True
    assert results[1]["model"] == "gpt-3.5-turbo"
    assert results[1]["response"] == "Response from gpt-3.5-turbo"


def test_run_parallel_models_handles_errors(monkeypatch):
    def fake_call_model(model_name, messages, temperature=0.2):
        if model_name == "bad-model":
            raise ValueError("Model failed")
        return "ok"

    monkeypatch.setattr(pr, "call_model", fake_call_model)

    results = pr.run_parallel_models("prompt", ["good-model", "bad-model"])
    assert len(results) == 2
    
    good = next(r for r in results if r["model"] == "good-model")
    assert good["success"] is True
    
    bad = next(r for r in results if r["model"] == "bad-model")
    assert bad["success"] is False
    assert "Model failed" in bad["error"]


def test_synthesize_responses_no_successful(monkeypatch):
    responses = [{"success": False}, {"success": False}]
    res = pr.synthesize_responses("q", responses)
    assert res == "No successful model responses to synthesize."


def test_synthesize_responses_single_success(monkeypatch):
    responses = [{"success": True, "response": "Single answer"}, {"success": False}]
    res = pr.synthesize_responses("q", responses)
    assert res == "Single answer"


def test_synthesize_responses_multiple_success(monkeypatch):
    responses = [
        {"success": True, "response": "Ans1", "model": "m1"},
        {"success": True, "response": "Ans2", "model": "m2"},
    ]
    
    monkeypatch.setattr(pr, "call_model", lambda *a, **k: "Synthesized Answer")
    res = pr.synthesize_responses("q", responses)
    assert res == "Synthesized Answer"


def test_synthesize_responses_fallback_on_error(monkeypatch):
    responses = [
        {"success": True, "response": "Ans1", "model": "m1"},
        {"success": True, "response": "Ans2", "model": "m2"},
    ]
    
    def fake_fail(*a, **k):
        raise RuntimeError("Synthesis failed")
        
    monkeypatch.setattr(pr, "call_model", fake_fail)
    
    # Should fallback to first successful response
    res = pr.synthesize_responses("q", responses)
    assert res == "Ans1"


def test_parallel_query_defaults(monkeypatch):
    monkeypatch.setattr(pr, "get_available_models", lambda: ["m1", "m2"])
    
    # Mock both stages
    monkeypatch.setattr(
        pr, 
        "run_parallel_models", 
        lambda p, models: [{"model": m, "response": "r", "success": True} for m in models]
    )
    monkeypatch.setattr(pr, "synthesize_responses", lambda q, res: "Final")
    
    result = pr.parallel_query("q")
    assert result["question"] == "q"
    assert result["synthesized_answer"] == "Final"
    assert len(result["individual_responses"]) == 2

