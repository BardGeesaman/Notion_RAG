from __future__ import annotations

"""
Tests for study grouping and metadata editing behavior.

Coverage:
- _ensure_experiment_for_study creates/reuses Experiments for a study
- Multiple datasets linking to the same Experiment
- Dataset objects reflect experiment membership
- Dataset update service applies all basic metadata fields
"""

from typing import List, Dict
from uuid import uuid4

import pytest

from amprenta_rag.api.schemas import DatasetUpdate
from amprenta_rag.api.services.datasets import update_dataset
from amprenta_rag.database.models import Dataset, Experiment, Program
from scripts.dashboard.pages.repositories import _ensure_experiment_for_study


class DummyQuery:
    """Minimal stand-in for SQLAlchemy Query for a single model type."""

    def __init__(self, result):
        self._result = result

    def filter(self, *args, **kwargs):
        # For these tests we ignore the filter expressions and just propagate self.
        return self

    def first(self):
        if isinstance(self._result, list):
            return self._result[0] if self._result else None
        return self._result


class DummyDB:
    """Dummy DB session for _ensure_experiment_for_study tests."""

    def __init__(self, existing_experiment: Experiment | None = None):
        self.existing_experiment = existing_experiment
        self.added: List[Experiment] = []
        self.commits: int = 0

    def query(self, model):
        # Only Experiment is queried in _ensure_experiment_for_study
        assert model is Experiment
        return DummyQuery(self.existing_experiment)

    def add(self, obj):
        self.added.append(obj)
        # Once added, subsequent queries should see this experiment
        self.existing_experiment = obj

    def commit(self):
        self.commits += 1

    def refresh(self, obj):
        # No-op for tests
        return None


def test_ensure_experiment_for_study_creates_new_experiment():
    """
    _ensure_experiment_for_study should create and persist a new Experiment
    when no matching study exists.
    """
    db = DummyDB(existing_experiment=None)
    repo = "GEO"
    accession = "GSE12345"
    title = "ALS CSF Transcriptomics Study"

    exp = _ensure_experiment_for_study(db=db, repo=repo, accession=accession, title=title)

    assert isinstance(exp, Experiment)
    assert exp in db.added
    assert db.commits == 1
    assert exp.description == f"{repo} study {accession}"
    # External IDs should contain the repository-specific key
    assert exp.external_ids.get("geo_accession") == accession


def test_ensure_experiment_for_study_reuses_existing_experiment():
    """
    When an Experiment for the study already exists, _ensure_experiment_for_study
    should return it without creating a new one.
    """
    existing = Experiment(
        id=uuid4(),
        name="Existing GEO Study",
        description="GEO study GSE12345",
        external_ids={"geo_accession": "GSE12345"},
    )
    db = DummyDB(existing_experiment=existing)

    exp = _ensure_experiment_for_study(db=db, repo="GEO", accession="GSE12345", title="Some title")

    assert exp is existing
    # No new experiments should have been added or committed
    assert db.added == []
    assert db.commits == 0


def test_multiple_datasets_link_to_same_experiment():
    """
    Multiple Dataset objects should be able to link to the same Experiment instance.
    """
    exp = Experiment(id=uuid4(), name="Grouped Study", description="Grouped study", external_ids={})

    ds1 = Dataset(id=uuid4(), name="DS1", omics_type="lipidomics")
    ds2 = Dataset(id=uuid4(), name="DS2", omics_type="proteomics")

    ds1.experiments.append(exp)
    ds2.experiments.append(exp)

    assert exp in ds1.experiments
    assert exp in ds2.experiments
    # Both datasets should refer to the exact same Experiment object
    assert ds1.experiments[0] is ds2.experiments[0] is exp


def test_dataset_experiment_membership_reflected_in_model():
    """
    Dataset detail view can reflect experiment membership via the ORM relationship.
    """
    db = DummyDB(existing_experiment=None)
    exp = _ensure_experiment_for_study(
        db=db,
        repo="GEO",
        accession="GSE99999",
        title="Study GSE99999",
    )

    dataset = Dataset(id=uuid4(), name="DS-Study-GSE99999", omics_type="transcriptomics")
    dataset.experiments.append(exp)

    experiment_names = [e.name for e in dataset.experiments]
    assert "Study GSE99999" in experiment_names


class DummyUpdateDB:
    """Dummy DB session for update_dataset tests."""

    def __init__(self, dataset: Dataset):
        self._dataset = dataset
        self.commits = 0
        self.refreshed = False

    def query(self, model):
        assert model is Dataset
        return DummyQuery(self._dataset)

    def commit(self):
        self.commits += 1

    def refresh(self, obj):
        assert obj is self._dataset
        self.refreshed = True


def test_metadata_update_saves_basic_fields(monkeypatch):
    """
    update_dataset should apply DatasetUpdate fields correctly to the ORM model.

    This approximates the behavior of the metadata editing form, which builds a
    DatasetUpdate from user input and calls update_dataset.
    """
    # Initial dataset
    ds = Dataset(
        id=uuid4(),
        name="Old Name",
        omics_type="lipidomics",
        description="Old description",
        file_paths=[],
        file_urls=[],
        organism=[],
        sample_type=[],
        disease=[],
    )
    db = DummyUpdateDB(ds)

    # Build update payload: change most fields
    update_payload = DatasetUpdate(
        name="New Name",
        description="New description",
        file_paths=["/data/new.csv"],
        file_urls=["http://example.com/new.csv"],
        organism=["Human"],
        sample_type=["CSF"],
        disease=["ALS"],
    )

    # Call service
    updated = update_dataset(db=db, dataset_id=ds.id, dataset=update_payload)

    assert updated is ds
    assert db.commits == 1
    assert db.refreshed is True

    assert ds.name == "New Name"
    assert ds.description == "New description"
    assert ds.file_paths == ["/data/new.csv"]
    assert ds.file_urls == ["http://example.com/new.csv"]
    assert ds.organism == ["Human"]
    assert ds.sample_type == ["CSF"]
    assert ds.disease == ["ALS"]


