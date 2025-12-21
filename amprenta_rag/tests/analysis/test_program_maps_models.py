from __future__ import annotations

from amprenta_rag.analysis.program_maps import models


def test_program_signature_score_defaults():
    ps = models.ProgramSignatureScore(
        program_id="p",
        signature_id="s",
        program_name="P",
        signature_name="S",
        overall_score=0.4,
    )
    assert ps.score_by_omics == {}
    assert ps.matching_datasets == []
    assert ps.coverage_fraction == 0.0


def test_program_signature_map_defaults():
    m = models.ProgramSignatureMap(
        program_id="p",
        program_name="P",
    )
    assert m.signature_scores == []
    assert m.top_signatures == []
    assert m.convergence_indicators == {}

