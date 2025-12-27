"""Tests for SAR grid API endpoints."""

from __future__ import annotations

from types import SimpleNamespace
from uuid import uuid4

import pytest
from fastapi.testclient import TestClient

from amprenta_rag.api.main import app


pytest.importorskip("fastapi")

# Check RDKit availability for tests that need it
import importlib.util
RDKIT_AVAILABLE = importlib.util.find_spec("rdkit") is not None

client = TestClient(app)


@pytest.mark.skipif(not RDKIT_AVAILABLE, reason="RDKit required for scaffold extraction")
def test_scaffolds_endpoint(monkeypatch) -> None:
    """Test GET /api/v1/sar/scaffolds endpoint."""
    from amprenta_rag.api.routers import sar as sar_router
    
    # Mock database compounds
    mock_compounds = [
        SimpleNamespace(smiles="c1ccccc1CCO"),  # Benzene scaffold
        SimpleNamespace(smiles="c1ccccc1CCC"),  # Benzene scaffold
        SimpleNamespace(smiles="c1ccc2ccccc2c1C"),  # Naphthalene scaffold
    ]
    
    def mock_get_db():
        db = SimpleNamespace()
        db.query = lambda model: SimpleNamespace(
            filter=lambda condition: SimpleNamespace(
                all=lambda: mock_compounds
            )
        )
        yield db
    
    app.dependency_overrides[sar_router.get_db] = mock_get_db
    
    try:
        resp = client.get("/api/v1/sar/scaffolds")
        assert resp.status_code == 200
        
        data = resp.json()
        assert isinstance(data, list)
        
        # Should have scaffold summaries
        if data:  # Only check if RDKit is available
            scaffold = data[0]
            assert "scaffold_smiles" in scaffold
            assert "compound_count" in scaffold
            assert isinstance(scaffold["compound_count"], int)
            assert scaffold["compound_count"] > 0
            
    finally:
        del app.dependency_overrides[sar_router.get_db]


@pytest.mark.skipif(RDKIT_AVAILABLE, reason="Test only applies when RDKit is not available")
def test_scaffolds_endpoint_no_rdkit() -> None:
    """Test scaffolds endpoint when RDKit is not available."""
    # This test only runs when RDKit is actually not available
    # If RDKit is available, this test is skipped
    resp = client.get("/api/v1/sar/scaffolds")
    assert resp.status_code == 503
    assert "RDKit not available" in resp.json()["detail"]


def test_grid_endpoint_with_compounds(monkeypatch) -> None:
    """Test POST /api/v1/sar/grid endpoint with valid compounds."""
    from amprenta_rag.api.routers import sar as sar_router
    
    # Mock database compounds
    compound_ids = [str(uuid4()), str(uuid4())]
    mock_compounds = [
        SimpleNamespace(compound_id=compound_ids[0], smiles="c1ccccc1CCO", molecular_weight=122.0),
        SimpleNamespace(compound_id=compound_ids[1], smiles="c1ccccc1CCC", molecular_weight=120.0),
    ]
    
    def mock_get_db():
        db = SimpleNamespace()
        db.query = lambda model: SimpleNamespace(
            filter=lambda condition: SimpleNamespace(
                all=lambda: mock_compounds
            )
        )
        yield db
    
    # Mock R-group functions
    def mock_find_common_core(smiles_list):
        return "c1ccccc1"  # Benzene core
    
    def mock_decompose_rgroups(smiles_list, core_smiles):
        return [
            {"smiles": "c1ccccc1CCO", "R1": "CCO", "R2": "H"},
            {"smiles": "c1ccccc1CCC", "R1": "CCC", "R2": "H"},
        ]
    
    monkeypatch.setattr(sar_router, "find_common_core", mock_find_common_core)
    monkeypatch.setattr(sar_router, "decompose_rgroups", mock_decompose_rgroups)
    
    app.dependency_overrides[sar_router.get_db] = mock_get_db
    
    try:
        resp = client.post(
            "/api/v1/sar/grid",
            json={
                "compound_ids": compound_ids,
                "x_axis": "R1",
                "y_axis": "R2"
            }
        )
        assert resp.status_code == 200
        
        data = resp.json()
        assert "matrix" in data
        assert "row_labels" in data
        assert "col_labels" in data
        assert "core_smiles" in data
        assert "total_compounds" in data
        
        assert data["core_smiles"] == "c1ccccc1"
        assert data["total_compounds"] == 2
        assert isinstance(data["matrix"], list)
        assert isinstance(data["row_labels"], list)
        assert isinstance(data["col_labels"], list)
        
    finally:
        del app.dependency_overrides[sar_router.get_db]


def test_grid_endpoint_auto_core(monkeypatch) -> None:
    """Test grid endpoint with automatic core detection."""
    from amprenta_rag.api.routers import sar as sar_router
    
    compound_ids = [str(uuid4())]
    mock_compounds = [
        SimpleNamespace(compound_id=compound_ids[0], smiles="c1ccccc1CCO", molecular_weight=122.0),
    ]
    
    def mock_get_db():
        db = SimpleNamespace()
        db.query = lambda model: SimpleNamespace(
            filter=lambda condition: SimpleNamespace(
                all=lambda: mock_compounds
            )
        )
        yield db
    
    # Track if find_common_core was called
    core_detection_called = False
    
    def mock_find_common_core(smiles_list):
        nonlocal core_detection_called
        core_detection_called = True
        return "c1ccccc1"
    
    def mock_decompose_rgroups(smiles_list, core_smiles):
        return [{"smiles": "c1ccccc1CCO", "R1": "CCO", "R2": "H"}]
    
    monkeypatch.setattr(sar_router, "find_common_core", mock_find_common_core)
    monkeypatch.setattr(sar_router, "decompose_rgroups", mock_decompose_rgroups)
    
    app.dependency_overrides[sar_router.get_db] = mock_get_db
    
    try:
        resp = client.post(
            "/api/v1/sar/grid",
            json={
                "compound_ids": compound_ids,
                # No core_smarts provided - should auto-detect
                "x_axis": "R1",
                "y_axis": "R2"
            }
        )
        assert resp.status_code == 200
        assert core_detection_called, "find_common_core should have been called for auto-detection"
        
    finally:
        del app.dependency_overrides[sar_router.get_db]


def test_grid_endpoint_manual_core(monkeypatch) -> None:
    """Test grid endpoint with manual core specification."""
    from amprenta_rag.api.routers import sar as sar_router
    
    compound_ids = [str(uuid4())]
    mock_compounds = [
        SimpleNamespace(compound_id=compound_ids[0], smiles="c1ccccc1CCO", molecular_weight=122.0),
    ]
    
    def mock_get_db():
        db = SimpleNamespace()
        db.query = lambda model: SimpleNamespace(
            filter=lambda condition: SimpleNamespace(
                all=lambda: mock_compounds
            )
        )
        yield db
    
    # Track if find_common_core was called (should NOT be called)
    core_detection_called = False
    
    def mock_find_common_core(smiles_list):
        nonlocal core_detection_called
        core_detection_called = True
        return "c1ccccc1"
    
    def mock_decompose_rgroups(smiles_list, core_smiles):
        assert core_smiles == "c1ccncc1", "Manual core should be passed to decompose_rgroups"
        return [{"smiles": "c1ccccc1CCO", "R1": "CCO", "R2": "H"}]
    
    monkeypatch.setattr(sar_router, "find_common_core", mock_find_common_core)
    monkeypatch.setattr(sar_router, "decompose_rgroups", mock_decompose_rgroups)
    
    app.dependency_overrides[sar_router.get_db] = mock_get_db
    
    try:
        resp = client.post(
            "/api/v1/sar/grid",
            json={
                "compound_ids": compound_ids,
                "core_smarts": "c1ccncc1",  # Manual core provided
                "x_axis": "R1",
                "y_axis": "R2"
            }
        )
        assert resp.status_code == 200
        assert not core_detection_called, "find_common_core should NOT be called when manual core provided"
        
    finally:
        del app.dependency_overrides[sar_router.get_db]


def test_grid_endpoint_invalid_compounds() -> None:
    """Test grid endpoint with non-existent compound IDs."""
    from amprenta_rag.api.routers import sar as sar_router
    
    def mock_get_db():
        db = SimpleNamespace()
        db.query = lambda model: SimpleNamespace(
            filter=lambda condition: SimpleNamespace(
                all=lambda: []  # No compounds found
            )
        )
        yield db
    
    app.dependency_overrides[sar_router.get_db] = mock_get_db
    
    try:
        resp = client.post(
            "/api/v1/sar/grid",
            json={
                "compound_ids": ["non-existent-id-1", "non-existent-id-2"],
                "x_axis": "R1",
                "y_axis": "R2"
            }
        )
        assert resp.status_code == 404
        assert "No compounds found" in resp.json()["detail"]
        
    finally:
        del app.dependency_overrides[sar_router.get_db]


def test_grid_endpoint_no_valid_smiles() -> None:
    """Test grid endpoint with compounds that have no valid SMILES."""
    from amprenta_rag.api.routers import sar as sar_router
    
    compound_ids = [str(uuid4())]
    mock_compounds = [
        SimpleNamespace(compound_id=compound_ids[0], smiles=None, molecular_weight=122.0),  # No SMILES
    ]
    
    def mock_get_db():
        db = SimpleNamespace()
        db.query = lambda model: SimpleNamespace(
            filter=lambda condition: SimpleNamespace(
                all=lambda: mock_compounds
            )
        )
        yield db
    
    app.dependency_overrides[sar_router.get_db] = mock_get_db
    
    try:
        resp = client.post(
            "/api/v1/sar/grid",
            json={
                "compound_ids": compound_ids,
                "x_axis": "R1",
                "y_axis": "R2"
            }
        )
        assert resp.status_code == 400
        assert "No valid SMILES found" in resp.json()["detail"]
        
    finally:
        del app.dependency_overrides[sar_router.get_db]


def test_grid_endpoint_core_detection_fails(monkeypatch) -> None:
    """Test grid endpoint when core detection fails."""
    from amprenta_rag.api.routers import sar as sar_router
    
    compound_ids = [str(uuid4())]
    mock_compounds = [
        SimpleNamespace(compound_id=compound_ids[0], smiles="CCO", molecular_weight=46.0),  # Simple alcohol
    ]
    
    def mock_get_db():
        db = SimpleNamespace()
        db.query = lambda model: SimpleNamespace(
            filter=lambda condition: SimpleNamespace(
                all=lambda: mock_compounds
            )
        )
        yield db
    
    def mock_find_common_core(smiles_list):
        return None  # Core detection fails
    
    monkeypatch.setattr(sar_router, "find_common_core", mock_find_common_core)
    
    app.dependency_overrides[sar_router.get_db] = mock_get_db
    
    try:
        resp = client.post(
            "/api/v1/sar/grid",
            json={
                "compound_ids": compound_ids,
                "x_axis": "R1",
                "y_axis": "R2"
            }
        )
        assert resp.status_code == 400
        assert "Could not find common core" in resp.json()["detail"]
        
    finally:
        del app.dependency_overrides[sar_router.get_db]


def test_grid_endpoint_decomposition_fails(monkeypatch) -> None:
    """Test grid endpoint when R-group decomposition fails."""
    from amprenta_rag.api.routers import sar as sar_router
    
    compound_ids = [str(uuid4())]
    mock_compounds = [
        SimpleNamespace(compound_id=compound_ids[0], smiles="c1ccccc1CCO", molecular_weight=122.0),
    ]
    
    def mock_get_db():
        db = SimpleNamespace()
        db.query = lambda model: SimpleNamespace(
            filter=lambda condition: SimpleNamespace(
                all=lambda: mock_compounds
            )
        )
        yield db
    
    def mock_find_common_core(smiles_list):
        return "c1ccccc1"
    
    def mock_decompose_rgroups(smiles_list, core_smiles):
        return []  # Decomposition fails
    
    monkeypatch.setattr(sar_router, "find_common_core", mock_find_common_core)
    monkeypatch.setattr(sar_router, "decompose_rgroups", mock_decompose_rgroups)
    
    app.dependency_overrides[sar_router.get_db] = mock_get_db
    
    try:
        resp = client.post(
            "/api/v1/sar/grid",
            json={
                "compound_ids": compound_ids,
                "x_axis": "R1",
                "y_axis": "R2"
            }
        )
        assert resp.status_code == 400
        assert "R-group decomposition failed" in resp.json()["detail"]
        
    finally:
        del app.dependency_overrides[sar_router.get_db]


def test_grid_endpoint_matrix_build_fails(monkeypatch) -> None:
    """Test grid endpoint when SAR matrix building fails."""
    from amprenta_rag.api.routers import sar as sar_router
    
    compound_ids = [str(uuid4())]
    mock_compounds = [
        SimpleNamespace(compound_id=compound_ids[0], smiles="c1ccccc1CCO", molecular_weight=122.0),
    ]
    
    def mock_get_db():
        db = SimpleNamespace()
        db.query = lambda model: SimpleNamespace(
            filter=lambda condition: SimpleNamespace(
                all=lambda: mock_compounds
            )
        )
        yield db
    
    def mock_find_common_core(smiles_list):
        return "c1ccccc1"
    
    def mock_decompose_rgroups(smiles_list, core_smiles):
        # Return data that will cause matrix building to fail
        return [{"smiles": "c1ccccc1CCO"}]  # Missing R-group columns
    
    monkeypatch.setattr(sar_router, "find_common_core", mock_find_common_core)
    monkeypatch.setattr(sar_router, "decompose_rgroups", mock_decompose_rgroups)
    
    app.dependency_overrides[sar_router.get_db] = mock_get_db
    
    try:
        resp = client.post(
            "/api/v1/sar/grid",
            json={
                "compound_ids": compound_ids,
                "x_axis": "R1",
                "y_axis": "R2"
            }
        )
        assert resp.status_code == 400
        assert "SAR grid generation failed" in resp.json()["detail"]
        
    finally:
        del app.dependency_overrides[sar_router.get_db]
