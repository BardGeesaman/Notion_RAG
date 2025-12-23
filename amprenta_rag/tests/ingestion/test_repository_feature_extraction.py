from __future__ import annotations

from types import SimpleNamespace
from pathlib import Path


import sys
from types import ModuleType

fake_geoparse = ModuleType("GEOparse")
sys.modules["GEOparse"] = fake_geoparse

from amprenta_rag.ingestion import repository_feature_extraction as rfe


def test_extract_geo_features_handles_error(monkeypatch):
    # Force GEOparse.get_GEO to raise to ensure graceful handling
    fake_geoparse.get_GEO = lambda *a, **k: (_ for _ in ()).throw(RuntimeError("fail"))
    genes = rfe.extract_geo_features_with_geoparse("GSE1", download_dir=Path("./tmp_geo"))
    assert genes == set()


def test_extract_pride_proteins_no_files(monkeypatch):
    class FakeRepo:
        def fetch_study_data_files(self, study_id):
            return []

    import amprenta_rag.ingestion.repositories.pride as pride_mod

    monkeypatch.setattr(pride_mod, "PRIDERepository", lambda: FakeRepo())
    proteins = rfe.extract_pride_proteins_from_data_files("PXD1", download_dir=Path("./tmp_pride"))
    assert proteins == set()


def test_extract_metabolights_metabolites_no_files(monkeypatch, tmp_path):
    class FakeRepo:
        def fetch_study_data_files(self, study_id):
            return []

    import amprenta_rag.ingestion.repositories.metabolights as mbl_mod

    monkeypatch.setattr(mbl_mod, "MetaboLightsRepository", lambda: FakeRepo())
    mets = rfe.extract_metabolights_metabolites_from_isa_tab("MTBLS1", download_dir=tmp_path)
    assert mets == set()


def test_extract_mw_metabolites_handles_error(monkeypatch):
    import amprenta_rag.ingestion.repositories.mw as mw_mod
    import requests

    monkeypatch.setattr(mw_mod, "REPOSITORY_USER_AGENT", "ua")
    monkeypatch.setattr(
        mw_mod.requests,
        "get",
        lambda *a, **k: (_ for _ in ()).throw(requests.exceptions.RequestException("boom")),
    )
    mets = rfe.extract_mw_metabolites_from_data_endpoint("ST000001")
    assert mets == set()


def test_extract_pride_proteins_mztab_success(monkeypatch, tmp_path):
    file_obj = SimpleNamespace(
        filename="results.mzTab",
        file_type="SEARCH",
        download_url="https://example.com/results.mzTab",
    )

    class FakeRepo:
        def fetch_study_data_files(self, study_id):
            return [file_obj]

    import amprenta_rag.ingestion.repositories.pride as pride_mod
    import requests

    monkeypatch.setattr(pride_mod, "PRIDERepository", lambda: FakeRepo())
    monkeypatch.setattr(pride_mod, "REPOSITORY_USER_AGENT", "ua")

    # Patch proteomics normalization import target
    fake_norm = ModuleType("amprenta_rag.ingestion.proteomics.normalization")
    fake_norm.normalize_protein_identifier = lambda s: s
    monkeypatch.setitem(sys.modules, "amprenta_rag.ingestion.proteomics.normalization", fake_norm)

    class Resp:
        status_code = 200
        text = "PRH\taccession\tother\nPRT\tP12345\tX\nPRT\tQ9ZZZ9\tY\n"
        content = text.encode("utf-8")

        def raise_for_status(self):
            return None

    monkeypatch.setattr(requests, "get", lambda *a, **k: Resp())

    proteins = rfe.extract_pride_proteins_from_data_files("PXD1", download_dir=tmp_path)
    assert proteins == {"P12345", "Q9ZZZ9"}


def test_extract_metabolights_metabolites_success(monkeypatch, tmp_path):
    class FakeMeta:
        def __init__(self):
            self.raw_metadata = {"mtblsStudy": {"studyHttpUrl": "https://example.com/MTBLS1"}}

    class FakeRepo:
        def fetch_study_metadata(self, study_id):
            return FakeMeta()

    import amprenta_rag.ingestion.repositories.metabolights as mbl_mod

    monkeypatch.setattr(mbl_mod, "MetaboLightsRepository", lambda: FakeRepo())
    monkeypatch.setattr(mbl_mod, "REPOSITORY_USER_AGENT", "ua")

    # metabolomics normalization
    fake_norm = ModuleType("amprenta_rag.ingestion.metabolomics.normalization")
    fake_norm.normalize_metabolite_name = lambda s: s.strip()
    monkeypatch.setitem(sys.modules, "amprenta_rag.ingestion.metabolomics.normalization", fake_norm)

    inv_url = "https://example.com/MTBLS1/i_Investigation.txt"
    maf_url = "https://example.com/MTBLS1/m_MTBLS1_maf.tsv"

    class InvResp:
        status_code = 200
        text = "some stuff m_MTBLS1_maf.tsv other"

    class MafResp:
        status_code = 200

        def raise_for_status(self):
            return None

        def iter_content(self, chunk_size=8192):
            data = (
                "metabolite_identification\tsample name\tS1\n"
                "Glucose\ts1\t1\n"
                "Lactate\ts1\t2\n"
            ).encode("utf-8")
            yield data

    def fake_get(url, *a, **k):
        if url == inv_url:
            return InvResp()
        if url == maf_url:
            return MafResp()
        raise AssertionError(f"unexpected url {url}")

    monkeypatch.setattr(mbl_mod.requests, "get", fake_get)
    # ensure HEAD never called (we found via investigation)
    monkeypatch.setattr(mbl_mod.requests, "head", lambda *a, **k: SimpleNamespace(status_code=404))

    mets = rfe.extract_metabolights_metabolites_from_isa_tab("MTBLS1", download_dir=tmp_path)
    assert mets == {"Glucose", "Lactate"}


def test_extract_mw_metabolites_success(monkeypatch):
    import amprenta_rag.ingestion.repositories.mw as mw_mod

    monkeypatch.setattr(mw_mod, "REPOSITORY_USER_AGENT", "ua")

    fake_norm = ModuleType("amprenta_rag.ingestion.metabolomics.normalization")
    fake_norm.normalize_metabolite_name = lambda s: s.strip()
    monkeypatch.setitem(sys.modules, "amprenta_rag.ingestion.metabolomics.normalization", fake_norm)

    class Resp:
        status_code = 200

        def json(self):
            return {
                "1": {"metabolite_name": "Glucose"},
                "2": {"refmet_name": "Lactate"},
            }

    monkeypatch.setattr(mw_mod.requests, "get", lambda *a, **k: Resp())
    mets = rfe.extract_mw_metabolites_from_data_endpoint("ST000001")
    assert mets == {"Glucose", "Lactate"}

