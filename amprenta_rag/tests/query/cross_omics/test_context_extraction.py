from __future__ import annotations

from amprenta_rag.query.cross_omics import context_extraction as ce


def test_extract_context_helpers():
    page = {"properties": {}}
    assert ce.extract_disease_context(page) == []
    assert ce.extract_matrix_context(page) == []
    assert ce.extract_model_system_context(page) == []


def test_extract_aggregated_context_returns_empty_and_logs():
    ctx = ce.extract_aggregated_context(["p1", "p2"])
    assert ctx["page_count"] == 0
    assert ctx["diseases"] == []


def test_format_and_identify_comparative_context():
    ctx = {"diseases": ["Flu"], "matrix": ["Plasma"], "model_systems": ["patient"], "page_count": 2}
    formatted = ce.format_context_for_prompt(ctx)
    assert "Flu" in formatted and "Plasma" in formatted

    comp = ce.identify_comparative_context(ctx)
    assert comp and comp["suggests_comparison"] is True

