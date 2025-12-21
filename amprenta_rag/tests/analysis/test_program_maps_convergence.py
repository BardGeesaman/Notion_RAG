from __future__ import annotations

from amprenta_rag.analysis.program_maps import convergence
from amprenta_rag.analysis.program_maps.models import ProgramSignatureScore


def _score(overall=0.5, omics=None):
    return ProgramSignatureScore(
        program_id="p",
        signature_id="s",
        program_name="P",
        signature_name="S",
        overall_score=overall,
        score_by_omics=omics or {},
    )


def test_compute_convergence_indicators_empty():
    assert convergence.compute_convergence_indicators([]) == {}


def test_compute_convergence_indicators_counts_multi_omics():
    scores = [
        _score(omics={"rna": 0.1, "prot": 0.2}),
        _score(omics={"rna": 0.0, "prot": 0.3}),
    ]
    metrics = convergence.compute_convergence_indicators(scores)
    assert metrics["multi_omics_signature_count"] == 1
    assert metrics["convergence_fraction"] == 0.5
    assert metrics["avg_omics_per_signature"] == 1.5

