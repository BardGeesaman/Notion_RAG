"""
Tier 1 end-to-end workflow validation (high-level, heavily mocked).

Workflow under test:
    ingest (Postgres models) → warm feature cache → score signatures → cross-omics summary

This test avoids real Notion/Pinecone/OpenAI usage and focuses on verifying that:
- Postgres-backed objects flow through cache warming and scoring
- The feature cache is populated and used
- The cross-omics program summary (Postgres version) runs end-to-end
- Enhanced disease / matrix / model_systems context reaches the LLM synthesis layer
- No Notion HTTP calls are made
"""

from __future__ import annotations

from typing import Any, Dict, Iterable, List
from uuid import uuid4

import pytest
from unittest.mock import patch

from amprenta_rag.database.models import Dataset, Experiment, Program
from amprenta_rag.ingestion.dataset_feature_cache import clear_feature_cache, get_feature_cache
from amprenta_rag.ingestion.multi_omics_scoring import extract_dataset_features_by_type
from amprenta_rag.ingestion.batch_signature_scoring import score_datasets_against_signatures_batch
from amprenta_rag.query.cross_omics.program_summary_postgres import (
    cross_omics_program_summary_postgres,
)


class _DummyProgramSession:
    """Dummy DB session that always returns a single Program instance by ID."""

    def __init__(self, program: Program):
        self._program = program

    class _Query:
        def __init__(self, result: Any):
            self._result = result

        def filter(self, *args: Any, **kwargs: Any) -> "_DummyProgramSession._Query":
            return self

        def first(self) -> Any:
            return self._result

    def query(self, model: Any) -> "_DummyProgramSession._Query":
        # program_summary_postgres only queries Program
        return self._Query(self._program)

    def close(self) -> None:  # pragma: no cover - trivial
        return None


def _fake_get_db(session: _DummyProgramSession):
    try:
        yield session
    finally:
        session.close()


def _make_program_with_datasets() -> Program:
    """Create a Program with four linked datasets (one per omics type) and rich context."""
    program = Program(
        id=uuid4(),
        name="ALS-CSF Multi-Omics Program",
        disease=["ALS"],
    )
    exp = Experiment(
        id=uuid4(),
        name="ALS-CSF-Patient-Study",
        disease=["ALS"],
        matrix=["CSF"],
        model_systems=["patient"],
    )

    omics_types = ["lipidomics", "metabolomics", "proteomics", "transcriptomics"]
    for idx, omics in enumerate(omics_types, start=1):
        ds = Dataset(
            id=uuid4(),
            name=f"Dataset-{omics}",
            omics_type=omics,
            disease=["ALS"],
            notion_page_id=f"dataset-{idx}",
        )
        program.datasets.append(ds)
        exp.datasets.append(ds)

    program.experiments.append(exp)
    return program


def _fake_features_for_dataset(dataset_id: str) -> Dict[str, set[str]]:
    return {
        "gene": {f"TP53-{dataset_id}"},
        "protein": {f"P04637-{dataset_id}"},
        "metabolite": set(),
        "lipid": set(),
    }


def _warm_cache_for_datasets(dataset_ids: Iterable[str]):
    """
    Simple stand-in for scripts.warm_feature_cache.warm_cache, using the real
    DatasetFeatureCache but a fake extract function (no external I/O).
    """
    cache = get_feature_cache()
    for ds_id in dataset_ids:
        cache.set_features(ds_id, _fake_features_for_dataset(ds_id))


def test_tier1_end_to_end_pipeline(monkeypatch):
    """
    End-to-end smoke test for Tier 1:
    - Program + datasets exist in "Postgres"
    - Cache is warmed for all datasets
    - Batch scoring uses cache
    - Cross-omics program summary runs and carries enhanced context
    - No Notion HTTP calls are made
    """
    clear_feature_cache()
    program = _make_program_with_datasets()
    dataset_page_ids = [ds.notion_page_id for ds in program.datasets if ds.notion_page_id]

    # 1) "Ingest": program + datasets are already present as ORM objects
    assert len(program.datasets) == 4

    # 2) Warm cache for all datasets (no external calls)
    _warm_cache_for_datasets(dataset_page_ids)
    cache_stats = get_feature_cache().get_stats()
    assert cache_stats["cached_datasets"] == len(dataset_page_ids)

    # 3) Score signatures against cached datasets
    score_calls: List[Dict[str, Any]] = []

    def fake_find_matching_signatures_for_dataset(dataset_page_id: str, overlap_threshold: float = 0.3):
        score_calls.append(
            {
                "dataset_page_id": dataset_page_id,
                "overlap_threshold": overlap_threshold,
            }
        )
        # Return a trivial "match" structure
        return [{"dataset_id": dataset_page_id, "score": 0.9}]

    with patch(
        "amprenta_rag.ingestion.signature_matching.find_matching_signatures_for_dataset",
        side_effect=fake_find_matching_signatures_for_dataset,
    ):
        results = score_datasets_against_signatures_batch(
            dataset_page_ids=dataset_page_ids,
            preload_cache=False,
            use_cache=True,
        )

    assert set(results.keys()) == set(dataset_page_ids)
    assert len(score_calls) == len(dataset_page_ids)

    # 4) Cross-omics program summary using Postgres-backed function
    session = _DummyProgramSession(program)
    monkeypatch.setattr(
        "amprenta_rag.database.base.get_db",
        lambda: _fake_get_db(session),
    )

    captured_prompt: Dict[str, Any] = {"value": ""}

    def fake_build_enhanced_prompt(
        entity_name: str,
        entity_type: str,
        context_info: Dict[str, Any] | None,
        omics_counts: Dict[str, int],
        additional_info: str,
    ) -> str:
        # Store a human-readable prompt for assertions
        captured_prompt["value"] = (
            f"{entity_name} | {entity_type} | diseases={context_info.get('diseases', [])} "
            f"matrix={context_info.get('matrix', [])} model_systems={context_info.get('model_systems', [])} "
            f"omics={omics_counts}"
            if context_info
            else f"{entity_name} | {entity_type} | omics={omics_counts}"
        )
        return captured_prompt["value"]

    # Patch out Pinecone + LLM, and guard against any HTTP calls
    with patch(
        "amprenta_rag.query.cross_omics.helpers.retrieve_chunks_for_objects",
        return_value=[],
    ), patch(
        "amprenta_rag.query.cross_omics.program_summary_postgres.build_enhanced_prompt",
        side_effect=fake_build_enhanced_prompt,
    ), patch(
        "amprenta_rag.query.cross_omics.synthesis.synthesize_cross_omics_summary",
        side_effect=lambda prompt, context_chunks, max_chunks=32, include_comparative=False: f"SUMMARY: {prompt}",
    ), patch(
        "requests.get"
    ) as mock_get, patch(
        "requests.post"
    ) as mock_post:
        summary = cross_omics_program_summary_postgres(program.id)

    # 5) Verify enhanced context
    prompt_text = captured_prompt["value"]
    assert "ALS" in prompt_text
    assert "matrix=" in prompt_text
    assert "model_systems=" in prompt_text
    assert "lipidomics" in prompt_text or "transcriptomics" in prompt_text

    # Summary should reflect the same context (we echoed it into the summary)
    assert "SUMMARY:" in summary
    assert "ALS" in summary

    # 6) No Notion HTTP calls expected
    assert mock_get.call_count == 0
    assert mock_post.call_count == 0


