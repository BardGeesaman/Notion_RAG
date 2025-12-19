"""
Tests for chemistry-related FastAPI endpoints (compounds and screening).

These tests:
- Use a temporary SQLite database path for chemistry data
- Exercise core compounds endpoints and basic error handling
"""

from __future__ import annotations

import pytest
from fastapi import HTTPException
from fastapi.testclient import TestClient

from amprenta_rag.api.main import app


client = TestClient(app)


@pytest.fixture
def chemistry_service_stub(monkeypatch):
    """
    Stub chemistry service functions to avoid real database dependencies.
    """
    sample_compound = {
        "compound_id": "CMPD-001",
        "smiles": "CCO",
        "inchi_key": "INCHI-KEY-001",
        "canonical_smiles": "CCO",
        "molecular_formula": "C2H6O",
        "molecular_weight": 46.07,
        "logp": -0.3,
        "hbd_count": 1,
        "hba_count": 1,
        "rotatable_bonds": 1,
    }

    def fake_list_compounds():
        return [sample_compound]

    def fake_get_compound_by_id(compound_id: str):
        if compound_id == "CMPD-001":
            return sample_compound
        return None

    def fake_get_compound_programs(compound_id: str):
        if compound_id != "CMPD-001":
            return []
        return [
            {
                "program_id": "PROG-ALS",
                "status": "active",
                "notes": "ALS core compound",
                "created_at": "2024-01-01T00:00:00Z",
            }
        ]

    monkeypatch.setattr(
        "amprenta_rag.api.services.compounds.list_compounds",
        fake_list_compounds,
    )
    monkeypatch.setattr(
        "amprenta_rag.api.services.compounds.get_compound_by_id",
        fake_get_compound_by_id,
    )
    monkeypatch.setattr(
        "amprenta_rag.api.services.compounds.get_compound_programs",
        fake_get_compound_programs,
    )

    return sample_compound


class TestCompoundsAPI:
    """Tests for /api/v1/compounds endpoints."""

    def test_list_compounds(self, chemistry_service_stub):
        response = client.get("/api/v1/compounds/")
        assert response.status_code == 200
        data = response.json()
        assert isinstance(data, list)
        assert len(data) == 1
        assert data[0]["compound_id"] == "CMPD-001"

    def test_get_compound_by_id(self, chemistry_service_stub):
        response = client.get("/api/v1/compounds/CMPD-001")
        assert response.status_code == 200
        compound = response.json()
        assert compound["compound_id"] == "CMPD-001"
        assert compound["molecular_formula"] == "C2H6O"

    def test_get_compound_not_found(self, chemistry_service_stub):
        response = client.get("/api/v1/compounds/UNKNOWN")
        assert response.status_code == 404
        detail = response.json().get("detail")
        assert "not found" in detail.lower()

    def test_get_compound_programs(self, chemistry_service_stub):
        response = client.get("/api/v1/compounds/CMPD-001/programs")
        assert response.status_code == 200
        programs = response.json()
        assert isinstance(programs, list)
        assert len(programs) == 1
        assert programs[0]["program_id"] == "PROG-ALS"

    def test_database_error_handling(self, monkeypatch, chemistry_service_stub):
        """
        Simulate a database error and verify that FastAPI surfaces a 500 error.
        """

        def broken_list():
            raise HTTPException(status_code=500, detail="Simulated DB failure")

        monkeypatch.setattr(
            "amprenta_rag.api.services.compounds.list_compounds",
            broken_list,
        )

        response = client.get("/api/v1/compounds/")
        # Depending on how errors bubble, this may be 500 or 502; we assert 5xx.
        assert response.status_code >= 500


