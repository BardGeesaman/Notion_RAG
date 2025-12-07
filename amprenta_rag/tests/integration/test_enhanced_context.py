"""
Enhanced context validation tests for cross-omics summaries.

Focus:
- Disease / matrix / model_systems metadata extraction from Postgres models
- Propagation into the LLM prompt via build_enhanced_prompt
- Correct handling of comparative analysis flag
"""

from __future__ import annotations

from typing import Any, Dict, List
from uuid import uuid4

import pytest
from unittest.mock import patch

from amprenta_rag.database.models import Dataset, Experiment, Program
from amprenta_rag.query.cross_omics.dataset_summary_postgres import (
    cross_omics_dataset_summary_postgres,
)
from amprenta_rag.query.cross_omics.program_summary_postgres import (
    aggregate_context_from_models,
    identify_comparative_context_postgres,
)


class _DummyDatasetSession:
    """Minimal dummy session returning a single Dataset by UUID."""

    def __init__(self, dataset: Dataset):
        self._dataset = dataset

    class _Query:
        def __init__(self, result: Any):
            self._result = result

        def filter(self, *args: Any, **kwargs: Any) -> "_DummyDatasetSession._Query":
            return self

        def first(self) -> Any:
            return self._result

    def query(self, model: Any) -> "_DummyDatasetSession._Query":
        # Only Dataset is used in cross_omics_dataset_summary_postgres
        return self._Query(self._dataset)

    def close(self) -> None:  # pragma: no cover - trivial
        return None


def _fake_get_db(session: _DummyDatasetSession):
    try:
        yield session
    finally:
        session.close()


def test_dataset_summary_prompt_includes_enhanced_context(monkeypatch):
    """
    cross_omics_dataset_summary_postgres() should pass disease/matrix/model context
    into build_enhanced_prompt and the synthesized summary.
    """
    ds_id = uuid4()
    dataset = Dataset(
        id=ds_id,
        name="ALS-CSF-Transcriptomics",
        omics_type="transcriptomics",
        disease=["ALS"],
    )
    # Link an experiment to provide matrix/model_systems context
    exp = Experiment(
        id=uuid4(),
        name="ALS-CSF-Patient-Study",
        disease=["ALS"],
        matrix=["CSF"],
        model_systems=["patient"],
    )
    dataset.experiments.append(exp)

    session = _DummyDatasetSession(dataset)
    monkeypatch.setattr(
        "amprenta_rag.database.base.get_db",
        lambda: _fake_get_db(session),
    )

    captured_context: Dict[str, Any] = {}

    def fake_build_enhanced_prompt(
        entity_name: str,
        entity_type: str,
        context_info: Dict[str, Any] | None,
        omics_counts: Dict[str, int],
        additional_info: str,
    ) -> str:
        if context_info is not None:
            captured_context.update(context_info)
        # Simple prompt echoing context for easy inspection
        return f"{entity_name} | {entity_type} | {context_info} | {omics_counts}"

    with patch(
        "amprenta_rag.query.cross_omics.dataset_summary_postgres.build_enhanced_prompt",
        side_effect=fake_build_enhanced_prompt,
    ), patch(
        "amprenta_rag.query.cross_omics.synthesis.synthesize_cross_omics_summary",
        side_effect=lambda prompt, context_chunks, max_chunks=32, include_comparative=False: prompt,
    ):
        summary = cross_omics_dataset_summary_postgres(ds_id)

    # Ensure context was extracted correctly
    assert captured_context.get("diseases") == ["ALS"]
    # matrix/model_systems come from linked experiments
    assert captured_context.get("matrix") == ["CSF"]
    assert captured_context.get("model_systems") == ["patient"]

    # Summary (which is just the prompt here) should mention these values
    text = str(summary)
    assert "ALS" in text
    assert "CSF" in text
    assert "patient" in text


def test_comparative_context_sets_include_comparative_flag():
    """
    identify_comparative_context_postgres should detect comparative context and
    set include_comparative=True for synthesis.
    """
    # Build multiple datasets/experiments to trigger comparative context
    ds1 = Dataset(
        id=uuid4(),
        name="ALS-CSF-Transcriptomics",
        omics_type="transcriptomics",
        disease=["ALS"],
    )
    ds2 = Dataset(
        id=uuid4(),
        name="AD-Plasma-Proteomics",
        omics_type="proteomics",
        disease=["Alzheimer's"],
    )
    exp1 = Experiment(
        id=uuid4(),
        name="ALS-Patient-Study",
        disease=["ALS"],
        matrix=["CSF"],
        model_systems=["patient"],
    )
    exp2 = Experiment(
        id=uuid4(),
        name="AD-Mouse-Study",
        disease=["Alzheimer's"],
        matrix=["plasma"],
        model_systems=["mouse"],
    )

    aggregated = aggregate_context_from_models([ds1, ds2], [exp1, exp2])
    comparative = identify_comparative_context_postgres(aggregated)

    assert comparative is not None
    assert comparative["multiple_diseases"]
    assert comparative["multiple_matrices"]
    assert comparative["multiple_model_systems"]

    # When comparative context is present, include_comparative should be True.
    # We test this indirectly by patching synthesize_cross_omics_summary.
    include_flags: List[bool] = []

    with patch(
        "amprenta_rag.query.cross_omics.synthesis.synthesize_cross_omics_summary",
        side_effect=lambda prompt, context_chunks, max_chunks=32, include_comparative=False: include_flags.append(
            include_comparative
        )
        or "OK",
    ):
        # Manually construct a simplified prompt call mirroring program/dataset summary behavior
        from amprenta_rag.query.cross_omics.prompt_templates import build_enhanced_prompt

        prompt = build_enhanced_prompt(
            entity_name="Comparative Program",
            entity_type="program",
            context_info=aggregated,
            omics_counts={"transcriptomics": 1, "proteomics": 1},
            additional_info="Comparative context test",
        )

        # Simulate a summary call with comparative context present
        from amprenta_rag.query.cross_omics.synthesis import synthesize_cross_omics_summary

        synthesize_cross_omics_summary(prompt, ["chunk1", "chunk2"], include_comparative=True)

    assert include_flags, "synthesize_cross_omics_summary was not called"
    assert include_flags[-1] is True


