"""
Tests for automatic Program inference from dataset metadata.
"""

from __future__ import annotations

from typing import Any, Iterable, List
from uuid import uuid4


from amprenta_rag.database.models import Program
from amprenta_rag.ingestion.auto_linking import infer_program_from_metadata


class _DummySession:
    """Minimal stand-in for SQLAlchemy session with a pre-populated Program list."""

    def __init__(self, programs: List[Program]):
        self._programs = programs

    class _Query:
        def __init__(self, programs: List[Program]):
            self._programs = programs

        def all(self) -> List[Program]:
            return self._programs

    def query(self, model: Any) -> "_DummySession._Query":
        # Only Program is queried in infer_program_from_metadata
        return self._Query(self._programs)

    def close(self) -> None:  # pragma: no cover - trivial
        return None


def _fake_get_db(session: _DummySession):
    try:
        yield session
    finally:
        session.close()


def _program(name: str, diseases: Iterable[str]) -> Program:
    return Program(id=uuid4(), name=name, disease=list(diseases))


def test_infer_program_with_disease_overlap(monkeypatch):
    """
    When diseases overlap clearly with a single Program, infer_program_from_metadata
    should return that program ID with confidence >= min_confidence.
    """
    prog_als = _program("ALS-CSF Program", ["ALS"])
    prog_ad = _program("Alzheimer Plasma Program", ["Alzheimer"])
    session = _DummySession([prog_als, prog_ad])

    monkeypatch.setattr(
        "amprenta_rag.ingestion.auto_linking.get_db",
        lambda: _fake_get_db(session),
    )

    program_id, confidence = infer_program_from_metadata(
        diseases=["ALS"],
        keywords=[],
        filename=None,
        min_confidence=0.8,
    )

    assert program_id == str(prog_als.id)
    assert confidence >= 0.8


def test_infer_program_ambiguous_candidates_returns_none(monkeypatch):
    """
    If two candidates get the same top score, inference should return (None, score).
    """
    # Both programs share same disease, no distinguishing keywords
    prog1 = _program("ALS Program A", ["ALS"])
    prog2 = _program("ALS Program B", ["ALS"])
    session = _DummySession([prog1, prog2])

    monkeypatch.setattr(
        "amprenta_rag.ingestion.auto_linking.get_db",
        lambda: _fake_get_db(session),
    )

    program_id, confidence = infer_program_from_metadata(
        diseases=["ALS"],
        keywords=[],
        filename=None,
        min_confidence=0.0,
    )

    assert program_id is None
    # confidence should reflect the top score even in ambiguous case
    assert confidence > 0.0


def test_infer_program_below_confidence_threshold_rejects(monkeypatch):
    """
    When the best candidate score is below min_confidence, inference should reject.
    """
    # Single weak match: score will be positive but small
    weak_prog = _program("General Neurodegeneration Program", ["neurodegeneration"])
    session = _DummySession([weak_prog])

    monkeypatch.setattr(
        "amprenta_rag.ingestion.auto_linking.get_db",
        lambda: _fake_get_db(session),
    )

    program_id, confidence = infer_program_from_metadata(
        diseases=["ALS"],  # weak overlap based on heuristic
        keywords=[],
        filename=None,
        min_confidence=0.99,  # very high threshold
    )

    assert program_id is None
    # confidence should expose the (too-low) score
    assert 0.0 <= confidence < 0.99


def test_infer_program_keyword_matching_from_filename(monkeypatch):
    """
    Keywords derived from filename should contribute to confidence and allow matching.
    """
    prog = _program("ALS-Lipidomics Study", ["ALS"])
    other = _program("Unrelated Program", ["Cancer"])
    session = _DummySession([prog, other])

    monkeypatch.setattr(
        "amprenta_rag.ingestion.auto_linking.get_db",
        lambda: _fake_get_db(session),
    )

    filename = "ALS_Lipidomics_CSF_dataset.csv"

    program_id, confidence = infer_program_from_metadata(
        diseases=["ALS"],
        keywords=[],
        filename=filename,
        min_confidence=0.5,
    )

    assert program_id == str(prog.id)
    assert confidence >= 0.5


