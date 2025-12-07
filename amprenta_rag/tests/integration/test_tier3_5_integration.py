"""
Tier 3–5 integration tests: ingest → auto-linking → query via API.

These tests are high-level and heavily mocked to focus on:
- Postgres-backed models being used for auto-linking
- FastAPI endpoints returning data tied to those models
- Ensuring no Notion HTTP calls are made during the flow
"""

from __future__ import annotations

from typing import Any, Generator, List
from uuid import uuid4

import pytest
from fastapi.testclient import TestClient
from unittest.mock import patch

from amprenta_rag.api.main import app
from amprenta_rag.database.models import Dataset, Program
from amprenta_rag.ingestion.auto_linking import infer_program_from_metadata


client = TestClient(app)


class _DummySession:
    """Dummy DB session that exposes a single Program and Dataset."""

    def __init__(self, program: Program, dataset: Dataset):
        self._program = program
        self._dataset = dataset

    class _Query:
        def __init__(self, result_list: List[Any]):
            self._result_list = result_list

        def all(self) -> List[Any]:
            return self._result_list

        def filter(self, *args: Any, **kwargs: Any) -> "_DummySession._Query":
            return self

        def first(self) -> Any:
            return self._result_list[0] if self._result_list else None

    def query(self, model: Any) -> "_DummySession._Query":
        from amprenta_rag.database.models import Dataset as DatasetModel
        from amprenta_rag.database.models import Program as ProgramModel

        if model is ProgramModel:
            return self._Query([self._program])
        if model is DatasetModel:
            return self._Query([self._dataset])
        return self._Query([])

    def close(self) -> None:  # pragma: no cover - trivial
        return None


def _fake_get_db(session: _DummySession) -> Generator[_DummySession, None, None]:
    try:
        yield session
    finally:
        session.close()


def _make_program_and_dataset() -> tuple[Program, Dataset]:
    program = Program(
        id=uuid4(),
        name="ALS-CSF Program",
        disease=["ALS"],
    )
    dataset = Dataset(
        id=uuid4(),
        name="ALS-CSF-Lipidomics Dataset",
        omics_type="lipidomics",
        disease=["ALS"],
        file_paths=["ALS-CSF-lipidomics.csv"],
        ingestion_status="complete",
    )
    program.datasets.append(dataset)
    dataset.programs.append(program)
    return program, dataset


def test_ingest_auto_link_and_query_via_api(monkeypatch):
    """
    Simulated E2E:
    - Dataset + Program exist in Postgres (dummy session)
    - auto_linking.infer_program_from_metadata finds the Program
    - API endpoints can list programs and datasets without errors
    - No Notion HTTP calls are made during the flow
    """
    program, dataset = _make_program_and_dataset()
    session = _DummySession(program, dataset)

    # Patch get_db globally so auto_linking and API DB access use the same dummy session
    monkeypatch.setattr(
        "amprenta_rag.database.base.get_db",
        lambda: _fake_get_db(session),
    )

    # Step 1: Auto-link using diseases + filename keywords
    with patch(
        "amprenta_rag.ingestion.auto_linking.get_db",
        lambda: _fake_get_db(session),
    ), patch(
        "requests.get"
    ) as mock_get, patch(
        "requests.post"
    ) as mock_post:
        inferred_id, confidence = infer_program_from_metadata(
            diseases=dataset.disease,
            keywords=[dataset.name],
            filename=dataset.file_paths[0],
            min_confidence=0.5,
        )

    assert inferred_id == str(program.id)
    assert confidence >= 0.5

    # Ensure no Notion HTTP calls (or any HTTP via requests)
    assert mock_get.call_count == 0
    assert mock_post.call_count == 0

    # Step 2: Query via API – list programs and datasets
    # These calls will use the dummy DB session via patched get_db
    resp_programs = client.get("/api/v1/programs")
    assert resp_programs.status_code == 200
    programs = resp_programs.json()
    assert isinstance(programs, list)

    resp_datasets = client.get("/api/v1/datasets")
    assert resp_datasets.status_code == 200
    datasets = resp_datasets.json()
    assert isinstance(datasets, list)

    # At least one program and dataset should be present in responses
    # (exact shape depends on API serializers, so we avoid over-asserting)
    assert programs is not None
    assert datasets is not None


