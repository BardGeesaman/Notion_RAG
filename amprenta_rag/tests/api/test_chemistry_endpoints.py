"""
Tests for chemistry-related FastAPI endpoints (compounds and screening).

These tests:
- Use a temporary SQLite database path for chemistry data
- Exercise core compounds endpoints and basic error handling
"""

from __future__ import annotations

import sqlite3
from pathlib import Path

import pytest
from fastapi.testclient import TestClient

from amprenta_rag.api.main import app


client = TestClient(app)


def _init_temp_chemistry_db(tmp_path: Path) -> Path:
    """
    Initialize a minimal chemistry.db with just enough schema/data for tests.
    """
    db_path = tmp_path / "chemistry.db"
    conn = sqlite3.connect(str(db_path))
    cursor = conn.cursor()

    # Minimal schema subset for compounds and compound_program
    cursor.execute(
        """
        CREATE TABLE compounds (
            compound_id TEXT PRIMARY KEY,
            smiles TEXT,
            inchi_key TEXT,
            canonical_smiles TEXT,
            molecular_formula TEXT,
            molecular_weight REAL,
            logp REAL,
            hbd_count INTEGER,
            hba_count INTEGER,
            rotatable_bonds INTEGER,
            updated_at TEXT
        )
        """
    )

    cursor.execute(
        """
        CREATE TABLE compound_program (
            compound_id TEXT,
            program_id TEXT,
            status TEXT,
            notes TEXT,
            created_at TEXT
        )
        """
    )

    # Seed with a single compound and linked program
    cursor.execute(
        """
        INSERT INTO compounds (
            compound_id, smiles, inchi_key, canonical_smiles,
            molecular_formula, molecular_weight, logp,
            hbd_count, hba_count, rotatable_bonds, updated_at
        ) VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, CURRENT_TIMESTAMP)
        """,
        (
            "CMPD-001",
            "CCO",
            "INCHI-KEY-001",
            "CCO",
            "C2H6O",
            46.07,
            -0.3,
            1,
            1,
            1,
        ),
    )

    cursor.execute(
        """
        INSERT INTO compound_program (
            compound_id, program_id, status, notes, created_at
        ) VALUES (?, ?, ?, ?, CURRENT_TIMESTAMP)
        """,
        (
            "CMPD-001",
            "PROG-ALS",
            "active",
            "ALS core compound",
        ),
    )

    conn.commit()
    conn.close()
    return db_path


@pytest.fixture
def chemistry_db_env(monkeypatch, tmp_path):
    """
    Fixture that forces the chemistry DB path to point to a temporary SQLite file.
    """
    db_path = _init_temp_chemistry_db(tmp_path)

    def fake_get_chemistry_db_path():
        return db_path

    monkeypatch.setattr(
        "amprenta_rag.chemistry.database.get_chemistry_db_path",
        fake_get_chemistry_db_path,
    )
    monkeypatch.setattr(
        "amprenta_rag.api.services.compounds.get_chemistry_db_path",
        fake_get_chemistry_db_path,
    )

    return db_path


class TestCompoundsAPI:
    """Tests for /api/v1/compounds endpoints."""

    def test_list_compounds(self, chemistry_db_env):
        response = client.get("/api/v1/compounds/")
        assert response.status_code == 200
        data = response.json()
        assert isinstance(data, list)
        assert len(data) == 1
        assert data[0]["compound_id"] == "CMPD-001"

    def test_get_compound_by_id(self, chemistry_db_env):
        response = client.get("/api/v1/compounds/CMPD-001")
        assert response.status_code == 200
        compound = response.json()
        assert compound["compound_id"] == "CMPD-001"
        assert compound["molecular_formula"] == "C2H6O"

    def test_get_compound_not_found(self, chemistry_db_env):
        response = client.get("/api/v1/compounds/UNKNOWN")
        assert response.status_code == 404
        detail = response.json().get("detail")
        assert "not found" in detail.lower()

    def test_get_compound_programs(self, chemistry_db_env):
        response = client.get("/api/v1/compounds/CMPD-001/programs")
        assert response.status_code == 200
        programs = response.json()
        assert isinstance(programs, list)
        assert len(programs) == 1
        assert programs[0]["program_id"] == "PROG-ALS"

    def test_database_error_handling(self, monkeypatch, chemistry_db_env):
        """
        Simulate a database error and verify that FastAPI surfaces a 500 error.
        """

        def broken_connect(*args, **kwargs):
            raise sqlite3.OperationalError("Simulated DB failure")

        monkeypatch.setattr("sqlite3.connect", broken_connect)

        response = client.get("/api/v1/compounds/")
        # Depending on how errors bubble, this may be 500 or 502; we assert 5xx.
        assert response.status_code >= 500


