from __future__ import annotations

from contextlib import contextmanager
from types import SimpleNamespace
from uuid import uuid4
from unittest.mock import Mock

import pytest
from fastapi import HTTPException
from fastapi.testclient import TestClient


@pytest.mark.unit
def test_calculate_z_prime_known_controls():
    from amprenta_rag.analysis.hts_qc import calculate_z_prime

    # Perfect separation, zero variance -> Z' = 1.0
    pos = [1.0, 1.0, 1.0, 1.0]
    neg = [0.0, 0.0, 0.0, 0.0]
    assert calculate_z_prime(pos, neg) == 1.0


@pytest.mark.unit
def test_calculate_z_prime_returns_none_for_insufficient_or_zero_denom():
    from amprenta_rag.analysis.hts_qc import calculate_z_prime

    assert calculate_z_prime([1.0], [0.0, 0.0]) is None
    assert calculate_z_prime([1.0, 1.0], [0.0]) is None
    # denom == 0
    assert calculate_z_prime([1.0, 1.0], [1.0, 1.0]) is None


@pytest.mark.unit
def test_calculate_hit_rate_percentage_and_rounding():
    from amprenta_rag.analysis.hts_qc import calculate_hit_rate

    results = [
        SimpleNamespace(hit_flag=True),
        SimpleNamespace(hit_flag=False),
        SimpleNamespace(hit_flag=True),
        SimpleNamespace(hit_flag=False),
    ]
    assert calculate_hit_rate(results) == 50.0
    assert calculate_hit_rate([]) == 0.0


@pytest.mark.unit
def test_get_plate_qc_summary_returns_plateqcsummary(monkeypatch):
    from amprenta_rag.analysis import hts_qc

    cid = uuid4()
    # Include controls with normalized_value None to ensure they are excluded from control counts
    results = [
        SimpleNamespace(
            campaign_id=cid,
            well_position="A01",
            normalized_value=1.0,
            z_score=3.0,
            hit_flag=True,
            hit_category="pos",
            compound_id=uuid4(),
            result_id="R1",
        ),
        SimpleNamespace(
            campaign_id=cid,
            well_position="A02",
            normalized_value=0.0,
            z_score=-1.0,
            hit_flag=False,
            hit_category="neg",
            compound_id=uuid4(),
            result_id="R2",
        ),
        SimpleNamespace(
            campaign_id=cid,
            well_position="A03",
            normalized_value=None,
            z_score=None,
            hit_flag=False,
            hit_category="pos",
            compound_id=uuid4(),
            result_id="R3",
        ),
    ]

    db = SimpleNamespace(
        query=lambda model: SimpleNamespace(filter=lambda *args, **kwargs: SimpleNamespace(all=lambda: results))
    )

    @contextmanager
    def fake_db_session():
        yield db

    monkeypatch.setattr(hts_qc, "db_session", fake_db_session)

    summary = hts_qc.get_plate_qc_summary(cid)
    assert summary.campaign_id == cid
    assert summary.total_wells == 3
    assert summary.hits == 1
    assert summary.hit_rate == 33.33  # 1 / 3 * 100, rounded to 2 decimals
    # Only one pos control (normalized_value not None) and one neg control
    assert summary.pos_controls == 1
    assert summary.neg_controls == 1
    # insufficient controls (<2 each) => None
    assert summary.z_prime is None


@pytest.fixture()
def api_client():
    from amprenta_rag.api.main import app

    return TestClient(app)


@pytest.mark.api
def test_get_campaign_qc_endpoint_returns_schema(api_client, monkeypatch):
    from amprenta_rag.api.routers import hts as hts_router
    from amprenta_rag.analysis.hts_qc import PlateQCSummary

    cid = uuid4()
    monkeypatch.setattr(hts_router, "_campaign_exists", lambda _cid: None)
    monkeypatch.setattr(
        hts_router,
        "get_plate_qc_summary",
        Mock(return_value=PlateQCSummary(cid, 96, 1.23, 0.5, 2, 8, 8)),
    )

    resp = api_client.get(f"/api/v1/hts/campaigns/{cid}/qc")
    assert resp.status_code == 200
    body = resp.json()
    assert body["campaign_id"] == str(cid)
    assert body["total_wells"] == 96
    assert body["hit_rate"] == 1.23
    assert body["z_prime"] == 0.5


@pytest.mark.api
def test_get_campaign_plate_endpoint_returns_schema(api_client, monkeypatch):
    from amprenta_rag.api.routers import hts as hts_router
    from amprenta_rag.analysis.hts_qc import WellData

    cid = uuid4()
    monkeypatch.setattr(hts_router, "_campaign_exists", lambda _cid: None)
    monkeypatch.setattr(
        hts_router,
        "get_plate_heatmap_data",
        Mock(
            return_value=[
                WellData("A01", 0.1, 1.2, True, uuid4(), "R1"),
                WellData("A02", None, None, None, None, None),
            ]
        ),
    )

    resp = api_client.get(f"/api/v1/hts/campaigns/{cid}/plate")
    assert resp.status_code == 200
    body = resp.json()
    assert body[0]["well_position"] == "A01"
    assert body[0]["hit_flag"] is True
    assert body[1]["well_position"] == "A02"


@pytest.mark.api
def test_get_campaign_hits_endpoint_returns_schema(api_client, monkeypatch):
    from amprenta_rag.api.routers import hts as hts_router
    from amprenta_rag.analysis.hts_qc import HitCompound

    cid = uuid4()
    cmpd = uuid4()
    monkeypatch.setattr(hts_router, "_campaign_exists", lambda _cid: None)
    monkeypatch.setattr(
        hts_router,
        "get_hit_compounds",
        Mock(return_value=[HitCompound(result_id="R1", compound_id=cmpd, well_position="A01", normalized_value=0.1, z_score=3.3)]),
    )

    resp = api_client.get(f"/api/v1/hts/campaigns/{cid}/hits")
    assert resp.status_code == 200
    body = resp.json()
    assert body == [
        {
            "result_id": "R1",
            "compound_id": str(cmpd),
            "well_position": "A01",
            "normalized_value": 0.1,
            "z_score": 3.3,
        }
    ]


@pytest.mark.api
def test_campaign_missing_returns_404(api_client, monkeypatch):
    from amprenta_rag.api.routers import hts as hts_router

    cid = uuid4()

    def _raise(_cid):
        raise HTTPException(status_code=404, detail="HTS campaign not found")

    monkeypatch.setattr(hts_router, "_campaign_exists", _raise)

    resp = api_client.get(f"/api/v1/hts/campaigns/{cid}/qc")
    assert resp.status_code == 404
    assert resp.json()["detail"] == "HTS campaign not found"


