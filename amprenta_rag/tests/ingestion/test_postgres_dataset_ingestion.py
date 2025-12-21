from __future__ import annotations

from uuid import uuid4

from amprenta_rag.ingestion import postgres_dataset_ingestion as pdi


class _FakeDataset:
    def __init__(self):
        self.id = uuid4()
        self.name = "Dataset"
        self.description = "desc"
        self.omics_type = "metabolite"
        self.organism = ["human"]
        self.sample_type = ["plasma"]
        self.disease = ["covid"]
        self.programs = []
        self.experiments = []
        self.file_urls = ["http://example.com/file"]
        self.external_ids = {"doi": "10.1/xyz", "mw_study_id": "MW1"}
        self.methods = "methods"
        self.summary = "summary"
        self.results = "results"
        self.conclusions = "conclusions"
        self.dataset_source_type = "repo"
        self.data_origin = "external"


def test_build_dataset_text_content():
    ds = _FakeDataset()
    text = pdi.build_dataset_text_content(ds)
    assert "Dataset: Dataset" in text
    assert "Omics Type" in text
    assert "DOI" in text


def test_get_dataset_metadata_from_postgres():
    ds = _FakeDataset()
    meta = pdi.get_dataset_metadata_from_postgres(ds)
    assert meta["dataset_id"] == str(ds.id)
    assert meta["omics_type"] == "metabolite"
    assert meta["doi"] == "10.1/xyz"

