from __future__ import annotations

from uuid import uuid4

import pytest

from amprenta_rag.analysis import moa_inference as mi


class FakeResult:
    def __init__(self, target=None, ic50=None, ec50=None, kd=None, ki=None):
        self.target = target
        self.ic50 = ic50
        self.ec50 = ec50
        self.kd = kd
        self.ki = ki


class FakeFeature:
    def __init__(self, name):
        self.name = name


class FakeAssoc:
    def __init__(self, feature_id, dataset_id):
        self.c = type("C", (), {"feature_id": feature_id, "dataset_id": dataset_id})


class FakeQuery:
    def __init__(self, data):
        self._data = data

    def filter(self, *a, **k):
        return self

    def all(self):
        return list(self._data)


class FakeDB:
    def __init__(self, results=None, features=None):
        self.results = results or []
        self.features = features or []

    def query(self, model):
        if model is mi.BiochemicalResult:
            return FakeQuery(self.results)
        if model is mi.Feature:
            return FakeQuery(self.features)
        return FakeQuery([])

    def __enter__(self):
        return self

    def __exit__(self, *exc):
        return False


def test_generate_moa_candidates_collects_targets(monkeypatch):
    res = [FakeResult(target="T1"), FakeResult(target=None), FakeResult(target="T1")]
    monkeypatch.setattr(mi, "db_session", lambda: FakeDB(results=res))
    targets = mi.generate_moa_candidates(uuid4())
    assert targets == ["T1"]


def test_score_bioactivity_evidence_no_results(monkeypatch):
    monkeypatch.setattr(mi, "db_session", lambda: FakeDB(results=[]))
    assert mi.score_bioactivity_evidence(uuid4(), "T") == 0.0


def test_score_bioactivity_evidence_scoring(monkeypatch):
    res = [FakeResult(target="T", ic50=500), FakeResult(target="T", ki=2000)]
    monkeypatch.setattr(mi, "db_session", lambda: FakeDB(results=res))
    score = mi.score_bioactivity_evidence(uuid4(), "T")
    assert 0.0 < score <= 1.0


def test_collect_features_and_score_omics(monkeypatch):
    ds_id = uuid4()
    feat1 = FakeFeature("A")
    feat1.id = uuid4()
    feat2 = FakeFeature(None)
    features = [feat1, feat2]
    assoc = FakeAssoc(feat1.id, ds_id)

    def fake_query_features():
        q = FakeQuery(features)
        q.join = lambda *a, **k: q
        q.filter = lambda *a, **k: q
        return q

    monkeypatch.setattr(
        mi,
        "db_session",
        lambda: type(
            "CTX",
            (),
            {"__enter__": lambda s: s, "__exit__": lambda *a: False, "query": lambda self, model: fake_query_features()},
        )(),
    )
    feats = mi._collect_features([ds_id])
    assert "A" in feats

    monkeypatch.setattr(mi, "_collect_features", lambda ids: {"A"})
    assert mi.score_omics_concordance([ds_id], "A") == 0.6
    assert mi.score_omics_concordance([ds_id], "B") == 0.0


def test_score_pathway_enrichment_handles_errors(monkeypatch):
    monkeypatch.setattr(mi, "perform_pathway_enrichment", lambda *a, **k: (_ for _ in ()).throw(RuntimeError("x")))
    assert mi.score_pathway_enrichment(set(), "A") == 0.0
    assert mi.score_pathway_enrichment({"A"}, "A") == 0.0

    monkeypatch.setattr(mi, "perform_pathway_enrichment", lambda *a, **k: [])
    assert mi.score_pathway_enrichment({"A"}, "A") == 0.0

    class FakeRes:
        adjusted_p_value = 0.2

    monkeypatch.setattr(mi, "perform_pathway_enrichment", lambda *a, **k: [FakeRes()])
    assert mi.score_pathway_enrichment({"A"}, "B") > 0.0


def test_fuse_evidence_scores_bounds():
    contribs = [
        mi.EvidenceContribution("a", value=1.0, weight=0.5),
        mi.EvidenceContribution("b", value=0.0, weight=0.5),
    ]
    score = mi.fuse_evidence_scores(contribs)
    assert 0.0 <= score <= 1.0
    assert score == 0.5


def test_infer_moa_ranks(monkeypatch):
    monkeypatch.setattr(mi, "_collect_features", lambda ids: {"A"})
    monkeypatch.setattr(mi, "generate_moa_candidates", lambda cid: ["A", "B"])
    monkeypatch.setattr(mi, "score_bioactivity_evidence", lambda cid, cand: 0.5 if cand == "A" else 0.1)
    monkeypatch.setattr(mi, "score_omics_concordance", lambda ids, cand: 0.6 if cand == "A" else 0.0)
    monkeypatch.setattr(mi, "score_pathway_enrichment", lambda feats, cand: 0.4 if cand == "A" else 0.0)

    ranked = mi.infer_moa(uuid4(), [uuid4()])
    assert ranked[0].candidate_id == "A"
    assert ranked[0].rank == 1

