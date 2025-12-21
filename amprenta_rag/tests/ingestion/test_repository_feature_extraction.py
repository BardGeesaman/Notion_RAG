from __future__ import annotations

from types import SimpleNamespace
from pathlib import Path

import pytest

from amprenta_rag.ingestion import repository_feature_extraction as rfe


def test_extract_geo_features_handles_error(monkeypatch):
    # Force GEOparse.get_GEO to raise to ensure graceful handling
    monkeypatch.setattr(rfe, "GEOparse", SimpleNamespace(get_GEO=lambda *a, **k: (_ for _ in ()).throw(RuntimeError("fail"))))
    genes = rfe.extract_geo_features_with_geoparse("GSE1", download_dir=Path("./tmp_geo"))
    assert genes == set()


def test_extract_pride_proteins_no_files(monkeypatch):
    class FakeRepo:
        def fetch_study_data_files(self, study_id):
            return []

    monkeypatch.setattr(rfe, "PRIDERepository", FakeRepo)
    proteins = rfe.extract_pride_proteins_from_data_files("PXD1", download_dir=Path("./tmp_pride"))
    assert proteins == set()


def test_extract_metabolights_metabolites_no_files(monkeypatch, tmp_path):
    class FakeRepo:
        def fetch_study_data_files(self, study_id):
            return []

    monkeypatch.setattr(rfe, "MetaboLightsRepository", FakeRepo)
    mets = rfe.extract_metabolites_from_metabolights("MTBLS1", download_dir=tmp_path)
    assert mets == set()


def test_extract_mw_metabolites_handles_error(monkeypatch):
    class FakeRepo:
        def fetch_metabolites(self, study_id):
            raise RuntimeError("boom")

    monkeypatch.setattr(rfe, "MWRepository", FakeRepo)
    mets = rfe.extract_mw_metabolites("MW1")
    assert mets == set()

