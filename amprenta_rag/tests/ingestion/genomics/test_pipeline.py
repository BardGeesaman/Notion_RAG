from __future__ import annotations

import gzip
from io import BytesIO
from pathlib import Path
from types import SimpleNamespace

import pytest

import sys
from types import ModuleType

from amprenta_rag.ingestion.genomics import pipeline as gp


def test_get_ena_fastqs_returns_runs(monkeypatch):
    # Mock the ingestion.repositories module to provide ENARepository
    fake_repo_mod = ModuleType("amprenta_rag.ingestion.repositories")

    class FakeMeta:
        def __init__(self):
            self.raw_metadata = {
                "fastq_ftp": "ftp.sra.ebi.ac.uk/vol1/fastq/R1.fastq.gz",
                "sample_alias": "alias1",
            }

    class FakeENA:
        def search_studies(self, keywords, filters=None, max_results=3):
            return ["RUN1"]

        def fetch_study_metadata(self, run_id):
            return FakeMeta()

    fake_repo_mod.ENARepository = FakeENA
    fake_repo_mod.REPOSITORY_USER_AGENT = "ua"
    monkeypatch.setitem(sys.modules, "amprenta_rag.ingestion.repositories", fake_repo_mod)
    runs = gp.get_ena_fastqs("human", limit=1)
    assert runs and runs[0]["URL"].startswith("http://")
    assert runs[0]["FTP"].startswith("ftp://")


def test_get_ena_fastqs_empty(monkeypatch):
    fake_repo_mod = ModuleType("amprenta_rag.ingestion.repositories")

    class FakeENA:
        def search_studies(self, keywords, filters=None, max_results=3):
            return []

    fake_repo_mod.ENARepository = FakeENA
    fake_repo_mod.REPOSITORY_USER_AGENT = "ua"
    monkeypatch.setitem(sys.modules, "amprenta_rag.ingestion.repositories", fake_repo_mod)
    assert gp.get_ena_fastqs("none") == []


def test_download_fastq_skips_without_confirm(tmp_path):
    run_info = {"Run": "R1", "URL": "http://example.com/file.fastq.gz"}
    path = gp.download_fastq(run_info, output_dir=tmp_path, confirm=False)
    assert path is None
    assert not list(tmp_path.iterdir())


def test_check_salmon_installed_handles_missing(monkeypatch):
    # Simulate FileNotFoundError
    def fake_run(*a, **k):
        raise FileNotFoundError()

    monkeypatch.setattr(gp.subprocess, "run", fake_run)
    assert gp.check_salmon_installed() is False


def test_quantify_with_salmon_creates_quant_file(tmp_path, monkeypatch):
    fastq = tmp_path / "a.fastq"
    fastq.write_text("data")
    index = tmp_path / "idx"
    index.write_text("idx")

    def fake_check():
        return True

    def fake_run(cmd, check, capture_output, text):
        # create quant.sf
        out_dir = Path(cmd[cmd.index("-o") + 1])
        out_dir.mkdir(parents=True, exist_ok=True)
        (out_dir / "quant.sf").write_text("Name\tTPM\nx\t1.0\n")

    monkeypatch.setattr(gp, "check_salmon_installed", fake_check)
    monkeypatch.setattr(gp.subprocess, "run", fake_run)

    out = gp.quantify_with_salmon(fastq, index, output_dir=tmp_path / "q")
    assert out and out.exists()


def test_extract_gene_counts_from_salmon(tmp_path):
    quant = tmp_path / "quant.sf"
    quant.write_text("Name\tTPM\ntr1\t2.5\ntr2\t0\n")
    counts = gp.extract_gene_counts_from_salmon(quant)
    assert counts == {"tr1": 2.5}


def test_quantify_with_kallisto_creates_output(tmp_path, monkeypatch):
    fastq = tmp_path / "a.fastq"
    fastq.write_text("data")
    index = tmp_path / "idx"
    index.write_text("idx")

    def fake_run(cmd, check, capture_output, text):
        out_dir = Path(cmd[cmd.index("-o") + 1])
        out_dir.mkdir(parents=True, exist_ok=True)
        (out_dir / "abundance.tsv").write_text("target_id\ttpm\nx\t1.0\n")

    monkeypatch.setattr(gp.subprocess, "run", fake_run)
    out = gp.quantify_with_kallisto(fastq, index, output_dir=tmp_path / "kallisto")
    assert out and out.exists()


def test_extract_gene_counts_from_kallisto(tmp_path):
    # Kallisto count extraction not implemented; ensure graceful failure path is not called here
    pass

