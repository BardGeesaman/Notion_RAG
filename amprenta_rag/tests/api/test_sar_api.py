"""Integration tests for SAR API endpoints."""

import pytest
from fastapi.testclient import TestClient
from unittest.mock import MagicMock, patch
from uuid import uuid4

from amprenta_rag.api.main import app
from amprenta_rag.api.dependencies import get_current_user


@pytest.fixture
def client():
    """Test client with mocked auth only."""
    def mock_user():
        user = MagicMock()
        user.id = uuid4()
        user.email = "test@example.com"
        return user
    
    app.dependency_overrides[get_current_user] = mock_user
    try:
        yield TestClient(app)
    finally:
        app.dependency_overrides.clear()


class TestSARAPI:
    """Integration tests for SAR API endpoints."""

    def test_list_targets_success(self, client):
        """GET /sar/targets returns target list."""
        with patch("amprenta_rag.api.routers.sar.service.list_targets") as mock_service:
            mock_service.return_value = [
                {"target_id": "EGFR", "target_name": "Epidermal Growth Factor Receptor", "compound_count": 150}
            ]
            
            response = client.get("/api/v1/sar/targets")
            assert response.status_code == 200
            data = response.json()
            assert len(data) == 1
            assert data[0]["target_id"] == "EGFR"

    def test_list_targets_with_limit(self, client):
        """GET /sar/targets respects limit parameter."""
        with patch("amprenta_rag.api.routers.sar.service.list_targets") as mock_service:
            mock_service.return_value = []
            
            response = client.get("/api/v1/sar/targets?limit=50")
            assert response.status_code == 200
            mock_service.assert_called_once_with(limit=50)

    def test_get_compounds_by_target_success(self, client):
        """GET /sar/targets/{target}/compounds returns compound list."""
        with patch("amprenta_rag.api.routers.sar.service.get_compounds_by_target") as mock_service:
            mock_service.return_value = [
                {"compound_id": "test123", "smiles": "CCO", "ic50": 10.5}
            ]
            
            response = client.get("/api/v1/sar/targets/EGFR/compounds")
            assert response.status_code == 200
            data = response.json()
            assert len(data) == 1
            assert data[0]["compound_id"] == "test123"

    def test_get_activity_cliffs_success(self, client):
        """GET /sar/targets/{target}/cliffs returns activity cliffs."""
        with patch("amprenta_rag.api.routers.sar.service.get_activity_cliffs_for_target") as mock_service:
            mock_service.return_value = [
                {
                    "compound_1": "CCO", 
                    "compound_2": "CCCO",
                    "similarity": 0.8,
                    "fold_change": 15.0
                }
            ]
            
            response = client.get("/api/v1/sar/targets/EGFR/cliffs")
            assert response.status_code == 200
            data = response.json()
            assert len(data) == 1
            assert data[0]["similarity"] == 0.8

    def test_get_activity_cliffs_with_params(self, client):
        """GET /sar/targets/{target}/cliffs respects query parameters."""
        with patch("amprenta_rag.api.routers.sar.service.get_activity_cliffs_for_target") as mock_service:
            mock_service.return_value = []
            
            response = client.get("/api/v1/sar/targets/EGFR/cliffs?similarity_threshold=0.7&fold_change=20.0&limit=25")
            assert response.status_code == 200
            mock_service.assert_called_once_with(
                target="EGFR",
                similarity_threshold=0.7,
                fold_change=20.0,
                limit=25
            )

    def test_validate_smiles_success(self, client):
        """POST /sar/validate validates SMILES."""
        with patch("amprenta_rag.api.routers.sar.validate_smiles") as mock_validate:
            mock_validate.return_value = True
            
            response = client.post("/api/v1/sar/validate", json={"smiles": "CCO"})
            assert response.status_code == 200
            data = response.json()
            assert data["smiles"] == "CCO"
            assert data["valid"] is True

    def test_validate_smiles_invalid(self, client):
        """POST /sar/validate handles invalid SMILES."""
        with patch("amprenta_rag.api.routers.sar.validate_smiles") as mock_validate:
            mock_validate.return_value = False
            
            response = client.post("/api/v1/sar/validate", json={"smiles": "INVALID"})
            assert response.status_code == 200
            data = response.json()
            assert data["valid"] is False

    def test_validate_smiles_import_error(self, client):
        """POST /sar/validate handles RDKit import error."""
        with patch("amprenta_rag.api.routers.sar.validate_smiles") as mock_validate:
            mock_validate.side_effect = ImportError("RDKit not available")
            
            response = client.post("/api/v1/sar/validate", json={"smiles": "CCO"})
            assert response.status_code == 503

    def test_predict_properties_success(self, client):
        """POST /sar/predict returns property predictions."""
        with patch("amprenta_rag.api.routers.sar.compare_properties") as mock_compare:
            mock_df = MagicMock()
            mock_df.to_dict.return_value = [
                {"smiles": "CCO", "mw": 46.07, "logp": -0.31}
            ]
            mock_compare.return_value = mock_df
            
            response = client.post("/api/v1/sar/predict", json={"smiles_list": ["CCO"]})
            assert response.status_code == 200
            data = response.json()
            assert len(data) == 1
            assert data[0]["smiles"] == "CCO"

    def test_predict_properties_import_error(self, client):
        """POST /sar/predict handles import error."""
        with patch("amprenta_rag.api.routers.sar.compare_properties") as mock_compare:
            mock_compare.side_effect = ImportError("RDKit not available")
            
            response = client.post("/api/v1/sar/predict", json={"smiles_list": ["CCO"]})
            assert response.status_code == 503

    def test_list_transformations_success(self, client):
        """GET /sar/transformations returns transformation list."""
        with patch("amprenta_rag.api.routers.sar.TRANSFORMATIONS", {"test_transform": {"name": "Test", "description": "Test transformation"}}):
            response = client.get("/api/v1/sar/transformations")
            assert response.status_code == 200
            data = response.json()
            assert len(data) == 1
            assert data[0]["id"] == "test_transform"

    def test_scaffold_hop_success(self, client):
        """POST /sar/scaffold-hop returns hopped products."""
        with patch("amprenta_rag.api.routers.sar.scaffold_hop") as mock_hop:
            mock_hop.return_value = ["CCCO", "CCCCO"]
            
            response = client.post("/api/v1/sar/scaffold-hop", json={
                "smiles": "CCO",
                "transformation": "test_transform"
            })
            assert response.status_code == 200
            data = response.json()
            assert data["smiles"] == "CCO"
            assert len(data["products"]) == 2

    def test_scaffold_hop_import_error(self, client):
        """POST /sar/scaffold-hop handles import error."""
        with patch("amprenta_rag.api.routers.sar.scaffold_hop") as mock_hop:
            mock_hop.side_effect = ImportError("RDKit not available")
            
            response = client.post("/api/v1/sar/scaffold-hop", json={
                "smiles": "CCO",
                "transformation": "test_transform"
            })
            assert response.status_code == 503

    @patch("amprenta_rag.api.routers.sar.get_db")
    def test_get_scaffolds_success(self, mock_get_db, client):
        """GET /sar/scaffolds returns scaffold summaries."""
        # Mock database session and compounds
        mock_db = MagicMock()
        mock_get_db.return_value = mock_db
        
        mock_compound = MagicMock()
        mock_compound.smiles = "CCO"
        mock_db.query.return_value.filter.return_value.all.return_value = [mock_compound]
        
        with patch("amprenta_rag.api.routers.sar.Chem") as mock_chem, \
             patch("amprenta_rag.api.routers.sar.Scaffolds") as mock_scaffolds:
            
            mock_mol = MagicMock()
            mock_chem.MolFromSmiles.return_value = mock_mol
            mock_scaffolds.MurckoScaffold.GetScaffoldForMol.return_value = mock_mol
            mock_chem.MolToSmiles.return_value = "c1ccccc1"
            
            response = client.get("/api/v1/sar/scaffolds")
            assert response.status_code == 200
            data = response.json()
            assert len(data) == 1
            assert data[0]["scaffold_smiles"] == "c1ccccc1"

    @patch("amprenta_rag.api.routers.sar.get_db")
    def test_get_scaffolds_rdkit_unavailable(self, mock_get_db, client):
        """GET /sar/scaffolds handles RDKit unavailable."""
        with patch("amprenta_rag.api.routers.sar.Chem", side_effect=ImportError):
            response = client.get("/api/v1/sar/scaffolds")
            assert response.status_code == 503

    @patch("amprenta_rag.api.routers.sar.get_db")
    def test_build_sar_grid_success(self, mock_get_db, client):
        """POST /sar/grid builds SAR matrix."""
        # Mock database session and compounds
        mock_db = MagicMock()
        mock_get_db.return_value = mock_db
        
        mock_compound = MagicMock()
        mock_compound.compound_id = "test123"
        mock_compound.smiles = "CCO"
        mock_compound.molecular_weight = 46.07
        mock_db.query.return_value.filter.return_value.all.return_value = [mock_compound]
        
        with patch("amprenta_rag.api.routers.sar.find_common_core") as mock_core, \
             patch("amprenta_rag.api.routers.sar.decompose_rgroups") as mock_decompose, \
             patch("amprenta_rag.api.routers.sar.build_sar_matrix") as mock_matrix:
            
            mock_core.return_value = "c1ccccc1"
            mock_decompose.return_value = [{"smiles": "CCO", "R1": "H", "R2": "H"}]
            
            mock_df = MagicMock()
            mock_df.empty = False
            mock_df.values.tolist.return_value = [[46.07]]
            mock_df.index.tolist.return_value = ["H"]
            mock_df.columns.tolist.return_value = ["H"]
            mock_matrix.return_value = mock_df
            
            response = client.post("/api/v1/sar/grid", json={
                "compound_ids": ["test123"],
                "core_smarts": "",
                "x_axis": "R1",
                "y_axis": "R2"
            })
            assert response.status_code == 200
            data = response.json()
            assert data["core_smiles"] == "c1ccccc1"
            assert data["total_compounds"] == 1

    @patch("amprenta_rag.api.routers.sar.get_db")
    def test_build_sar_grid_no_compounds(self, mock_get_db, client):
        """POST /sar/grid handles no compounds found."""
        mock_db = MagicMock()
        mock_get_db.return_value = mock_db
        mock_db.query.return_value.filter.return_value.all.return_value = []
        
        response = client.post("/api/v1/sar/grid", json={
            "compound_ids": ["nonexistent"],
            "core_smarts": "",
            "x_axis": "R1",
            "y_axis": "R2"
        })
        assert response.status_code == 404

    def test_sar_endpoints_require_auth(self):
        """SAR endpoints require authentication."""
        # Use client without auth override
        client = TestClient(app)
        
        # Test a few key endpoints
        response = client.get("/api/v1/sar/targets")
        assert response.status_code == 401
        
        response = client.post("/api/v1/sar/validate", json={"smiles": "CCO"})
        assert response.status_code == 401
        
        response = client.get("/api/v1/sar/scaffolds")
        assert response.status_code == 401
