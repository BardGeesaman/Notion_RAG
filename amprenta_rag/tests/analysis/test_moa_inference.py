from __future__ import annotations

from contextlib import contextmanager
from types import SimpleNamespace
from uuid import uuid4
from unittest.mock import Mock

import pytest
from fastapi import HTTPException
from fastapi.testclient import TestClient


@pytest.mark.unit
def test_generate_moa_candidates_returns_unique_targets(monkeypatch):
    from amprenta_rag.analysis import moa_inference as moa

    compound_id = uuid4()
    rows = [
        SimpleNamespace(target="CDK2"),
        SimpleNamespace(target="EGFR"),
        SimpleNamespace(target="CDK2"),
        SimpleNamespace(target=None),
        SimpleNamespace(target=""),
    ]

    db = SimpleNamespace(
        query=lambda model: SimpleNamespace(
            filter=lambda *args, **kwargs: SimpleNamespace(all=lambda: rows)
        )
    )

    @contextmanager
    def fake_db_session():
        yield db

    monkeypatch.setattr(moa, "db_session", fake_db_session)
    out = moa.generate_moa_candidates(compound_id)
    assert set(out) == {"CDK2", "EGFR"}


@pytest.mark.unit
def test_score_bioactivity_evidence_no_results_returns_0(monkeypatch):
    from amprenta_rag.analysis import moa_inference as moa

    compound_id = uuid4()

    db = SimpleNamespace(
        query=lambda model: SimpleNamespace(
            filter=lambda *args, **kwargs: SimpleNamespace(all=lambda: [])
        )
    )

    @contextmanager
    def fake_db_session():
        yield db

    monkeypatch.setattr(moa, "db_session", fake_db_session)
    assert moa.score_bioactivity_evidence(compound_id, "CDK2") == 0.0


@pytest.mark.unit
def test_score_bioactivity_evidence_no_numeric_values_returns_0_1(monkeypatch):
    from amprenta_rag.analysis import moa_inference as moa

    compound_id = uuid4()
    rows = [
        SimpleNamespace(ic50=None, ec50=None, kd=None, ki=None),
        SimpleNamespace(ic50=0, ec50=-1, kd=None, ki=None),
    ]

    db = SimpleNamespace(
        query=lambda model: SimpleNamespace(
            filter=lambda *args, **kwargs: SimpleNamespace(all=lambda: rows)
        )
    )

    @contextmanager
    def fake_db_session():
        yield db

    monkeypatch.setattr(moa, "db_session", fake_db_session)
    assert moa.score_bioactivity_evidence(compound_id, "CDK2") == 0.1


@pytest.mark.unit
def test_score_bioactivity_evidence_uses_best_potency(monkeypatch):
    from amprenta_rag.analysis import moa_inference as moa

    compound_id = uuid4()
    # best value is 10 -> score = 1/(1+10/1000)=~0.990099
    rows = [
        SimpleNamespace(ic50=1000.0, ec50=None, kd=None, ki=None),
        SimpleNamespace(ic50=None, ec50=10.0, kd=None, ki=None),
    ]

    db = SimpleNamespace(
        query=lambda model: SimpleNamespace(
            filter=lambda *args, **kwargs: SimpleNamespace(all=lambda: rows)
        )
    )

    @contextmanager
    def fake_db_session():
        yield db

    monkeypatch.setattr(moa, "db_session", fake_db_session)
    score = moa.score_bioactivity_evidence(compound_id, "CDK2")
    assert score == pytest.approx(1.0 / (1.0 + 10.0 / 1000.0))


@pytest.mark.unit
def test_fuse_evidence_scores_weighted_average_and_bounds():
    from amprenta_rag.analysis.moa_inference import EvidenceContribution, fuse_evidence_scores

    contribs = [
        EvidenceContribution(feature_name="a", value=1.0, weight=0.5),
        EvidenceContribution(feature_name="b", value=0.0, weight=0.5),
    ]
    assert fuse_evidence_scores(contribs) == 0.5

    assert fuse_evidence_scores([]) == 0.0
    assert fuse_evidence_scores([EvidenceContribution("x", value=1.0, weight=0.0)]) == 0.0
    # bounds
    assert fuse_evidence_scores([EvidenceContribution("x", value=10.0, weight=1.0)]) == 1.0


@pytest.mark.unit
def test_infer_moa_returns_ranked_candidates(monkeypatch):
    from amprenta_rag.analysis import moa_inference as moa

    compound_id = uuid4()
    ds_ids = [uuid4()]

    monkeypatch.setattr(moa, "_collect_features", Mock(return_value={"CDK2"}))
    monkeypatch.setattr(moa, "generate_moa_candidates", Mock(return_value=["CDK2", "EGFR"]))
    # Make CDK2 clearly higher
    monkeypatch.setattr(moa, "score_bioactivity_evidence", Mock(side_effect=[0.9, 0.1]))
    monkeypatch.setattr(moa, "score_omics_concordance", Mock(side_effect=[0.6, 0.0]))
    monkeypatch.setattr(moa, "score_pathway_enrichment", Mock(side_effect=[0.8, 0.0]))

    out = moa.infer_moa(compound_id, ds_ids)
    assert [c.candidate_id for c in out] == ["CDK2", "EGFR"]
    assert [c.rank for c in out] == [1, 2]
    assert out[0].probability >= out[1].probability
    assert out[0].type == "target"
    assert len(out[0].contributions) == 3


@pytest.fixture()
def api_client():
    from amprenta_rag.api.main import app

    return TestClient(app)


@pytest.mark.api
def test_post_moa_infer_returns_schema(api_client, monkeypatch):
    from amprenta_rag.api.routers import moa as moa_router
    from amprenta_rag.analysis.moa_inference import EvidenceContribution, MOACandidate

    compound_id = uuid4()
    ds_id = uuid4()

    monkeypatch.setattr(moa_router, "_validate_compound", lambda _id: None)
    monkeypatch.setattr(
        moa_router,
        "infer_moa",
        Mock(
            return_value=[
                MOACandidate(
                    candidate_id="CDK2",
                    type="target",
                    probability=0.9,
                    rank=1,
                    contributions=[EvidenceContribution("bioactivity", 0.9, 0.4)],
                )
            ]
        ),
    )

    resp = api_client.post(
        "/api/v1/moa/infer",
        json={"compound_id": str(compound_id), "dataset_ids": [str(ds_id)]},
    )
    assert resp.status_code == 200
    body = resp.json()
    assert body["compound_id"] == str(compound_id)
    assert body["candidates"][0]["candidate_id"] == "CDK2"
    assert body["candidates"][0]["contributions"][0]["feature_name"] == "bioactivity"


@pytest.mark.api
def test_get_compound_moa_returns_schema(api_client, monkeypatch):
    from amprenta_rag.api.routers import moa as moa_router
    from amprenta_rag.analysis.moa_inference import EvidenceContribution, MOACandidate

    compound_id = uuid4()

    monkeypatch.setattr(moa_router, "_validate_compound", lambda _id: None)
    monkeypatch.setattr(
        moa_router,
        "infer_moa",
        Mock(
            return_value=[
                MOACandidate(
                    candidate_id="EGFR",
                    type="target",
                    probability=0.7,
                    rank=1,
                    contributions=[EvidenceContribution("omics_concordance", 0.6, 0.3)],
                )
            ]
        ),
    )

    resp = api_client.get(f"/api/v1/compounds/{compound_id}/moa")
    assert resp.status_code == 200
    assert resp.json()["compound_id"] == str(compound_id)
    assert resp.json()["candidates"][0]["candidate_id"] == "EGFR"


@pytest.mark.api
def test_moa_endpoints_404_when_compound_missing(api_client, monkeypatch):
    from amprenta_rag.api.routers import moa as moa_router

    compound_id = uuid4()

    def _raise(_id):
        raise HTTPException(status_code=404, detail="Compound not found")

    monkeypatch.setattr(moa_router, "_validate_compound", _raise)

    resp = api_client.get(f"/api/v1/compounds/{compound_id}/moa")
    assert resp.status_code == 404
    assert resp.json()["detail"] == "Compound not found"


