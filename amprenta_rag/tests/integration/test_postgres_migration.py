"""
Postgres migration verification tests for cross-omics summaries.

These tests ensure that the new *postgres* summary functions:
- Use Postgres-backed models via `get_db`
- Do not make any Notion HTTP calls
- Are robust when Pinecone/OpenAI are fully mocked out
"""

from __future__ import annotations

from typing import Any, Dict, Generator, List, Optional, Type
from uuid import uuid4

import pytest
from unittest.mock import patch

from amprenta_rag.database.models import Dataset, Experiment, Feature, Program, Signature
from amprenta_rag.query.cross_omics.dataset_summary_postgres import (
    cross_omics_dataset_summary_postgres,
)
from amprenta_rag.query.cross_omics.feature_summary_postgres import (
    cross_omics_feature_summary_postgres,
)
from amprenta_rag.query.cross_omics.program_summary_postgres import (
    cross_omics_program_summary_postgres,
)
from amprenta_rag.query.cross_omics.signature_summary_postgres import (
    cross_omics_signature_summary_postgres,
)


class _DummyQuery:
    """Very small stand-in for SQLAlchemy query object."""

    def __init__(self, result: Any):
        self._result = result

    def filter(self, *args: Any, **kwargs: Any) -> "_DummyQuery":
        return self

    def first(self) -> Any:
        return self._result


class _DummySession:
    """Dummy DB session that returns pre-built model instances."""

    def __init__(
        self,
        program: Optional[Program] = None,
        dataset: Optional[Dataset] = None,
        signature: Optional[Signature] = None,
        feature: Optional[Feature] = None,
    ):
        self._program = program
        self._dataset = dataset
        self._signature = signature
        self._feature = feature

    def query(self, model: Type[Any]) -> _DummyQuery:
        if model is Program and self._program is not None:
            return _DummyQuery(self._program)
        if model is Dataset and self._dataset is not None:
            return _DummyQuery(self._dataset)
        if model is Signature and self._signature is not None:
            return _DummyQuery(self._signature)
        if model is Feature and self._feature is not None:
            return _DummyQuery(self._feature)
        # Default: return None so .first() → None
        return _DummyQuery(None)

    def close(self) -> None:  # pragma: no cover - trivial
        return None


def _fake_get_db(session: _DummySession) -> Generator[_DummySession, None, None]:
    """Generator matching amprenta_rag.database.base.get_db signature."""
    try:
        yield session
    finally:
        session.close()


def _make_program_with_context() -> Program:
    program = Program(
        id=uuid4(),
        name="ALS-CSF Program",
        disease=["ALS"],
    )
    # Attach one experiment and one dataset with additional context
    exp = Experiment(
        id=uuid4(),
        name="ALS-CSF-Patient-Study",
        disease=["ALS"],
        matrix=["CSF"],
        model_systems=["patient"],
    )
    ds = Dataset(
        id=uuid4(),
        name="ALS-CSF-Lipidomics",
        omics_type="lipidomics",
        disease=["ALS"],
        notion_page_id="dataset-page-1",
    )
    program.experiments.append(exp)
    program.datasets.append(ds)
    return program


def _make_dataset_with_context() -> Dataset:
    return Dataset(
        id=uuid4(),
        name="ALS-CSF-Metabolomics",
        omics_type="metabolomics",
        disease=["ALS"],
        notion_page_id="dataset-page-2",
    )


def _make_signature_with_context(ds: Dataset) -> Signature:
    sig = Signature(
        id=uuid4(),
        name="ALS-CSF-Core-6Ceramides",
        modalities=["lipid"],
        notion_page_id="signature-page-1",
    )
    sig.datasets.append(ds)
    return sig


def _make_feature_with_context(ds: Dataset) -> Feature:
    feat = Feature(
        id=uuid4(),
        name="TP53",
        feature_type="gene",
        normalized_name="TP53",
    )
    feat.datasets.append(ds)
    return feat


@pytest.mark.parametrize(
    "which_summary, build_objects",
    [
        (
            "program",
            _make_program_with_context,
        ),
        (
            "dataset",
            _make_dataset_with_context,
        ),
    ],
)
def test_program_and_dataset_summaries_use_postgres(monkeypatch, which_summary, build_objects):
    """
    cross_omics_program_summary_postgres and cross_omics_dataset_summary_postgres
    should complete successfully with:
    - get_db patched to a dummy Postgres session
    - Pinecone and OpenAI fully mocked
    - requests mocked and never called (no Notion API usage)
    """
    obj = build_objects()

    # Build dummy session for the specific object
    session = _DummySession(
        program=obj if which_summary == "program" else None,
        dataset=obj if which_summary == "dataset" else None,
    )

    # Patch get_db at the database.base level (used by all summary modules)
    monkeypatch.setattr(
        "amprenta_rag.database.base.get_db",
        lambda: _fake_get_db(session),
    )

    # Patch helpers that would hit Pinecone / LLM
    with patch(
        "amprenta_rag.query.cross_omics.helpers.retrieve_chunks_for_objects",
        return_value=[],
    ), patch(
        "amprenta_rag.query.cross_omics.synthesis.synthesize_cross_omics_summary",
        return_value="FAKE_SUMMARY",
    ), patch(
        "requests.get"
    ) as mock_get, patch(
        "requests.post"
    ) as mock_post:
        if which_summary == "program":
            result = cross_omics_program_summary_postgres(obj.id)
        else:
            result = cross_omics_dataset_summary_postgres(obj.id)

    # Function should return a non-empty summary string
    assert isinstance(result, str)
    assert result

    # No HTTP calls into Notion or other APIs should have been made
    assert mock_get.call_count == 0
    assert mock_post.call_count == 0


def test_signature_and_feature_summaries_use_postgres(monkeypatch):
    """
    cross_omics_signature_summary_postgres and cross_omics_feature_summary_postgres
    should complete successfully with:
    - get_db patched to a dummy Postgres session
    - Pinecone and OpenAI fully mocked
    - requests mocked and never called
    """
    ds = _make_dataset_with_context()
    signature = _make_signature_with_context(ds)
    feature = _make_feature_with_context(ds)

    session = _DummySession(signature=signature, feature=feature)

    monkeypatch.setattr(
        "amprenta_rag.database.base.get_db",
        lambda: _fake_get_db(session),
    )

    # Patch Pinecone + LLM stack
    with patch(
        "amprenta_rag.query.cross_omics.feature_summary_postgres.query_pinecone",
        return_value=[],
    ), patch(
        "amprenta_rag.query.cross_omics.signature_summary_postgres.query_pinecone",
        return_value=[],
    ), patch(
        "amprenta_rag.query.cross_omics.helpers.retrieve_chunks_for_objects",
        return_value=[],
    ), patch(
        "amprenta_rag.query.cross_omics.synthesis.synthesize_cross_omics_summary",
        return_value="FAKE_SUMMARY",
    ), patch(
        "requests.get"
    ) as mock_get, patch(
        "requests.post"
    ) as mock_post:
        sig_result = cross_omics_signature_summary_postgres(signature_id=signature.id)
        feat_result = cross_omics_feature_summary_postgres(feature_id=feature.id)

    assert isinstance(sig_result, str) and sig_result
    assert isinstance(feat_result, str) and feat_result

    # No HTTP calls should have been made (Postgres + Pinecone only)
    assert mock_get.call_count == 0
    assert mock_post.call_count == 0


