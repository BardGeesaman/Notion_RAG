from __future__ import annotations

from uuid import uuid4


from amprenta_rag.reporting import evidence_report as er


class _FakeQuery:
    def __init__(self, obj):
        self.obj = obj

    def filter(self, *a, **k):
        return self

    def first(self):
        return self.obj


class _FakeDB:
    def __init__(self, obj=None):
        self.obj = obj

    def query(self, model):
        return _FakeQuery(self.obj)

    def __enter__(self):
        return self

    def __exit__(self, *exc):
        return False


def test_generate_program_report_uses_program(monkeypatch):
    class FakeProgram:
        def __init__(self):
            self.name = "Prog"
            self.datasets = [type("D", (), {"id": uuid4()})()]
            self.signatures = [type("S", (), {"id": uuid4()})()]

    monkeypatch.setattr(er, "db_session", lambda: _FakeDB(FakeProgram()))
    report = er.generate_program_report(uuid4())
    assert report.entity_type == "program"
    assert report.sections[0].title == "Overview"


def test_generate_dataset_report(monkeypatch):
    class FakeDataset:
        def __init__(self):
            self.name = "DS"

    monkeypatch.setattr(er, "db_session", lambda: _FakeDB(FakeDataset()))
    report = er.generate_dataset_report(uuid4())
    assert report.entity_type == "dataset"
    assert "Dataset" in report.sections[0].summary_text


def test_generate_signature_report(monkeypatch):
    class FakeVal:
        def __init__(self):
            self.summary = "sum"
            self.matched_dataset_ids = [uuid4()]

    monkeypatch.setattr(er, "db_session", lambda: _FakeDB(type("Sig", (), {"name": "Sig"})()))
    monkeypatch.setattr(er, "validate_signature_against_all_datasets", lambda sid: FakeVal())
    report = er.generate_signature_report(uuid4())
    assert report.entity_type == "signature"
    assert any(s.title == "Validation Metrics" for s in report.sections)


def test_legacy_generators(monkeypatch):
    monkeypatch.setattr(er, "cross_omics_program_summary_postgres", lambda **k: "prog-summary")
    monkeypatch.setattr(er, "cross_omics_dataset_summary_postgres", lambda **k: "ds-summary")
    monkeypatch.setattr(er, "cross_omics_signature_summary_postgres", lambda **k: "sig-summary")
    monkeypatch.setattr(er, "cross_omics_feature_summary_postgres", lambda **k: "feat-summary")

    prog = er.generate_program_evidence_report(str(uuid4()))
    exp = er.generate_experiment_evidence_report(str(uuid4()))
    ds = er.generate_dataset_evidence_report(str(uuid4()))
    sig = er.generate_signature_evidence_report(str(uuid4()))
    feat = er.generate_feature_evidence_report("F", "gene")

    assert prog.summary == "prog-summary"
    assert exp.summary == "ds-summary"
    assert ds.summary == "ds-summary"
    assert sig.summary == "sig-summary"
    assert feat.metadata["feature_type"] == "gene"


def test_format_and_write_noop():
    rep = er.EvidenceReportLegacy(
        entity_type="dataset",
        entity_id="id",
        entity_name="Name",
        summary="summary",
        metadata={"k": "v"},
    )
    doc = er.format_evidence_report(rep, include_metadata=True)
    assert "Evidence Report" in doc and "K" in doc
    assert er.write_evidence_report_to_notion(rep) is False

