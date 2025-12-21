from __future__ import annotations

from types import SimpleNamespace

from amprenta_rag.ingestion.metadata import llm_semantic_extraction as lse


def test_extract_semantic_metadata_short_text_returns_empty():
    res = lse.extract_semantic_metadata_with_llm("too short")
    assert res["diseases"] == []
    assert res["targets"] == []


def test_extract_semantic_metadata_disabled_config(monkeypatch):
    cfg = SimpleNamespace(
        pipeline=SimpleNamespace(enable_llm_semantic_extraction=False),
        openai=SimpleNamespace(model="gpt"),
    )
    monkeypatch.setattr(lse, "get_config", lambda: cfg)
    res = lse.extract_semantic_metadata_with_llm("x" * 200)
    assert res["diseases"] == []


def test_enhance_metadata_merges(monkeypatch):
    monkeypatch.setattr(
        lse,
        "extract_semantic_metadata_with_llm",
        lambda text, source_type="generic", max_length=8000: {
            "diseases": ["Flu"],
            "targets": ["TP53"],
            "signatures": ["S1"],
            "phenotype_axes": ["Axis"],
            "biomarker_roles": ["diagnostic"],
        },
    )
    base = {
        "diseases": ["Flu"],
        "targets": [],
        "lipid_signatures": [],
        "phenotype_axes": [],
        "biomarker_role": [],
    }
    out = lse.enhance_metadata_with_llm(base, "x" * 200)
    assert out["targets"] == ["TP53"]

