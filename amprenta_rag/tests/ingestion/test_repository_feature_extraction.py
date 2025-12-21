from __future__ import annotations

from types import SimpleNamespace
from pathlib import Path

import pytest

import sys
from types import ModuleType

fake_geoparse = ModuleType("GEOparse")
sys.modules["GEOparse"] = fake_geoparse

from amprenta_rag.ingestion import repository_feature_extraction as rfe


def test_extract_geo_features_handles_error(monkeypatch):
    # Force GEOparse.get_GEO to raise to ensure graceful handling
    monkeypatch.setattr(rfe, "GEOparse", SimpleNamespace(get_GEO=lambda *a, **k: (_ for _ in ()).throw(RuntimeError("fail"))))
    genes = rfe.extract_geo_features_with_geoparse("GSE1", download_dir=Path("./tmp_geo"))
    assert genes == set()


def test_extract_pride_proteins_no_files(monkeypatch):
    fake_mod = ModuleType("amprenta_rag.ingestion.repositories.pride")

    class FakeRepo:
        def fetch_study_data_files(self, study_id):
            return []

    fake_mod.PRIDERepository = lambda: FakeRepo()
    monkeypatch.setitem(sys.modules, "amprenta_rag.ingestion.repositories.pride", fake_mod)
    proteins = rfe.extract_pride_proteins_from_data_files("PXD1", download_dir=Path("./tmp_pride"))
    assert proteins == set()


def test_extract_metabolights_metabolites_no_files(monkeypatch, tmp_path):
    fake_mod = ModuleType("amprenta_rag.ingestion.repositories.metabolights")

    class FakeRepo:
        def fetch_study_data_files(self, study_id):
            return []

    fake_mod.MetaboLightsRepository = lambda: FakeRepo()
    monkeypatch.setitem(sys.modules, "amprenta_rag.ingestion.repositories.metabolights", fake_mod)
    mets = rfe.extract_metabolights_metabolites_from_isa_tab("MTBLS1", download_dir=tmp_path)
    assert mets == set()


def test_extract_mw_metabolites_handles_error(monkeypatch):
    # MW extraction uses requests.get + REPOSITORY_USER_AGENT; mock both.
    fake_repo_root = ModuleType("amprenta_rag.ingestion.repositories")
    fake_repo_root.REPOSITORY_USER_AGENT = "ua"
    monkeypatch.setitem(sys.modules, "amprenta_rag.ingestion.repositories", fake_repo_root)

    monkeypatch.setattr(rfe.requests, "get", lambda *a, **k: (_ for _ in ()).throw(rfe.requests.exceptions.RequestException("boom")))
    mets = rfe.extract_mw_metabolites_from_data_endpoint("ST000001")
    assert mets == set()

