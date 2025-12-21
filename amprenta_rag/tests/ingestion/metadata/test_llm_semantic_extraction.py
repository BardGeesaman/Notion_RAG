from __future__ import annotations

from types import SimpleNamespace

from amprenta_rag.ingestion.metadata import llm_semantic_extraction as lse


def test_extract_metadata_from_llm_response_handles_missing():
    resp = {"choices": [{"message": {"content": '{"disease":"flu","organism":"human"}'}}]}
    extracted = lse.extract_metadata_from_llm_response(resp)
    assert extracted["disease"] == "flu"
    assert extracted["organism"] == "human"

    # Invalid JSON path should return empty dict
    resp_bad = {"choices": [{"message": {"content": "not json"}}]}
    assert lse.extract_metadata_from_llm_response(resp_bad) == {}


def test_build_prompt_context_includes_fields():
    ctx = lse.build_prompt_context("ds1", "summary", {"key": "val"})
    assert "ds1" in ctx
    assert "summary" in ctx
    assert "key" in ctx

