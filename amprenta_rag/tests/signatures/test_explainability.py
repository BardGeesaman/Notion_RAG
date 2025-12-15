from __future__ import annotations

from contextlib import contextmanager
from types import SimpleNamespace
from uuid import uuid4
from unittest.mock import Mock

import pytest
from fastapi.testclient import TestClient


@pytest.mark.unit
def test_componentmatch_contribution_computed_by_score_signature(monkeypatch):
    """
    Validate that ComponentMatch.contribution == weight * match_value (pre-normalization),
    for both match and conflict directions.
    """
    from amprenta_rag.signatures.signature_loader import Signature, SignatureComponent
    from amprenta_rag.signatures import signature_scoring as scoring

    # Make matching deterministic and avoid testing species-matching here.
    monkeypatch.setattr(scoring, "match_species", lambda dataset_species, signature_species, refmet_map=None: {s: s for s in signature_species if s in dataset_species})
    monkeypatch.setattr(scoring, "normalize_species_name", lambda s: s)

    sig = Signature(
        name="sig",
        components=[
            SignatureComponent(feature_name="A", direction="↑", weight=2.0),
            SignatureComponent(feature_name="B", direction="↓", weight=0.5),
        ],
    )
    dataset_species = {"A", "B"}
    dataset_directions = {"A": "↑", "B": "↑"}  # B conflicts with signature's ↓

    result = scoring.score_signature(sig, dataset_species, dataset_directions=dataset_directions)
    by_species = {m.signature_species: m for m in result.component_matches}

    assert by_species["A"].direction_match == "match"
    assert by_species["A"].contribution == pytest.approx(2.0 * 1.0)

    assert by_species["B"].direction_match == "conflict"
    assert by_species["B"].contribution == pytest.approx(0.5 * -1.0)


@pytest.mark.unit
def test_get_top_contributors_sorted_positive_and_negative():
    from amprenta_rag.signatures.signature_scoring import ComponentMatch, SignatureScoreResult

    matches = [
        ComponentMatch(signature_species="A", matched_dataset_species="A", contribution=0.2),
        ComponentMatch(signature_species="B", matched_dataset_species="B", contribution=1.5),
        ComponentMatch(signature_species="C", matched_dataset_species="C", contribution=-0.1),
        ComponentMatch(signature_species="D", matched_dataset_species="D", contribution=-2.0),
        ComponentMatch(signature_species="E", matched_dataset_species="E", contribution=0.0),
    ]
    res = SignatureScoreResult(total_score=0.0, component_matches=matches)

    top = res.get_top_contributors(n=10)
    assert [m.signature_species for m in top["positive"]] == ["B", "A"]
    # Negative list is sorted ascending (most negative first)
    assert [m.signature_species for m in top["negative"]] == ["D", "C"]


@pytest.mark.unit
def test_get_direction_concordance_percentage():
    from amprenta_rag.signatures.signature_scoring import ComponentMatch, SignatureScoreResult

    matches = [
        ComponentMatch(signature_species="A", matched_dataset_species="A", direction_match="match"),
        ComponentMatch(signature_species="B", matched_dataset_species="B", direction_match="conflict"),
        ComponentMatch(signature_species="C", matched_dataset_species=None, direction_match="missing"),
    ]
    res = SignatureScoreResult(total_score=0.0, component_matches=matches)
    assert res.get_direction_concordance() == pytest.approx(1 / 2)


@pytest.fixture()
def api_client():
    from amprenta_rag.api.main import app

    return TestClient(app)


class _Query:
    def __init__(self, first_result):
        self._first = first_result

    def filter(self, *args, **kwargs):
        return self

    def first(self):
        return self._first


@pytest.mark.api
def test_explain_endpoint_returns_expected_schema(api_client, monkeypatch):
    """
    Validate 200 response shape for:
      GET /api/v1/signatures/{id}/explain/{dataset_id}
    """
    from amprenta_rag.api.routers import explainability as exp_router
    from amprenta_rag.signatures.signature_scoring import ComponentMatch, SignatureScoreResult

    sig_id = uuid4()
    ds_id = uuid4()

    sig_model = SimpleNamespace(
        id=sig_id,
        name="TestSig",
        # Components are stored as dicts in the DB model and rehydrated into SignatureComponent
        components=[
            {"feature_name": "A", "feature_type": "gene", "direction": "↑", "weight": 2.0},
            {"feature_name": "B", "feature_type": "gene", "direction": "↓", "weight": 1.0},
        ],
    )
    dataset = SimpleNamespace(
        id=ds_id,
        name="TestDataset",
        features=[SimpleNamespace(name="A"), SimpleNamespace(name="B")],
    )

    # Fake DB session that returns signature then dataset
    def fake_db_query(model):
        if model is exp_router.Signature:
            return _Query(sig_model)
        if model is exp_router.Dataset:
            return _Query(dataset)
        raise AssertionError(f"Unexpected model queried: {model}")

    @contextmanager
    def fake_db_session():
        yield SimpleNamespace(query=fake_db_query)

    monkeypatch.setattr(exp_router, "db_session", fake_db_session)

    # Patch scoring to a deterministic result (we unit-test scoring separately)
    result = SignatureScoreResult(
        total_score=0.8,
        component_matches=[
            ComponentMatch(
                signature_species="A",
                matched_dataset_species="A",
                match_type="exact",
                direction_match="match",
                weight=2.0,
                contribution=2.0,
            ),
            ComponentMatch(
                signature_species="B",
                matched_dataset_species="B",
                match_type="exact",
                direction_match="conflict",
                weight=1.0,
                contribution=-1.0,
            ),
        ],
    )
    monkeypatch.setattr(exp_router, "score_signature", Mock(return_value=result))

    resp = api_client.get(f"/api/v1/signatures/{sig_id}/explain/{ds_id}?top_n=1")
    assert resp.status_code == 200
    body = resp.json()

    # Top-level schema fields
    assert body["signature_id"] == str(sig_id)
    assert body["signature_name"] == "TestSig"
    assert body["dataset_id"] == str(ds_id)
    assert body["dataset_name"] == "TestDataset"
    assert body["total_score"] == pytest.approx(0.8)

    # Concordance computed from result methods (1 match / 2 matched)
    assert body["direction_concordance"] == pytest.approx(0.5)

    # top_n=1 should cap lists
    assert len(body["top_positive"]) == 1
    assert len(body["top_negative"]) == 1

    # Ensure contribution objects are shaped
    pos = body["top_positive"][0]
    neg = body["top_negative"][0]
    assert set(pos.keys()) == {
        "feature_name",
        "matched_to",
        "direction_expected",
        "direction_actual",
        "weight",
        "contribution",
        "match_type",
        "direction_match",
    }
    assert pos["feature_name"] == "A"
    assert pos["contribution"] == pytest.approx(2.0)
    assert pos["direction_expected"] == "↑"
    assert neg["feature_name"] == "B"
    assert neg["contribution"] == pytest.approx(-1.0)
    assert neg["direction_expected"] == "↓"

    assert len(body["all_contributions"]) == 2


@pytest.mark.api
def test_explain_endpoint_404_for_invalid_signature(api_client, monkeypatch):
    from amprenta_rag.api.routers import explainability as exp_router

    sig_id = uuid4()
    ds_id = uuid4()

    def fake_db_query(model):
        if model is exp_router.Signature:
            return _Query(None)
        raise AssertionError("Dataset query should not be reached when signature missing")

    @contextmanager
    def fake_db_session():
        yield SimpleNamespace(query=fake_db_query)

    monkeypatch.setattr(exp_router, "db_session", fake_db_session)

    resp = api_client.get(f"/api/v1/signatures/{sig_id}/explain/{ds_id}")
    assert resp.status_code == 404
    assert resp.json()["detail"] == "Signature not found"


@pytest.mark.api
def test_explain_endpoint_404_for_invalid_dataset(api_client, monkeypatch):
    from amprenta_rag.api.routers import explainability as exp_router

    sig_id = uuid4()
    ds_id = uuid4()

    sig_model = SimpleNamespace(id=sig_id, name="Sig", components=[])

    def fake_db_query(model):
        if model is exp_router.Signature:
            return _Query(sig_model)
        if model is exp_router.Dataset:
            return _Query(None)
        raise AssertionError(f"Unexpected model queried: {model}")

    @contextmanager
    def fake_db_session():
        yield SimpleNamespace(query=fake_db_query)

    monkeypatch.setattr(exp_router, "db_session", fake_db_session)

    resp = api_client.get(f"/api/v1/signatures/{sig_id}/explain/{ds_id}")
    assert resp.status_code == 404
    assert resp.json()["detail"] == "Dataset not found"


