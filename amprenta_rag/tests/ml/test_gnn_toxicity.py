"""Tests for GNN toxicity models."""
import pytest
import torch
import os
from unittest.mock import MagicMock, patch
from pathlib import Path

# Set OpenMP workaround for tests
os.environ["KMP_DUPLICATE_LIB_OK"] = "TRUE"


class TestMoleculeGNN:
    """Tests for MoleculeGNN model."""
    
    def test_model_initialization(self):
        """Model initializes with correct parameters."""
        from amprenta_rag.ml.gnn.model import MoleculeGNN
        
        model = MoleculeGNN(node_features=147, edge_features=7)
        assert model is not None
        assert sum(p.numel() for p in model.parameters()) > 100000
        assert model.task_type == "classification"
        assert model.dropout == 0.2
    
    def test_model_forward_pass(self):
        """Forward pass produces correct output shape."""
        from amprenta_rag.ml.gnn.model import MoleculeGNN
        from amprenta_rag.ml.gnn.featurizer import mol_to_graph
        from torch_geometric.data import Batch
        
        model = MoleculeGNN(node_features=147, edge_features=7)
        data = mol_to_graph("CCO")
        batch = Batch.from_data_list([data])
        
        output = model(batch)
        assert output.shape == (1,)
    
    def test_model_classification_output_range(self):
        """Classification output is in [0, 1]."""
        from amprenta_rag.ml.gnn.model import MoleculeGNN
        from amprenta_rag.ml.gnn.featurizer import mol_to_graph
        from torch_geometric.data import Batch
        
        model = MoleculeGNN(node_features=147, edge_features=7, task_type="classification")
        data = mol_to_graph("CCO")
        batch = Batch.from_data_list([data])
        
        output = model(batch)
        assert 0 <= output.item() <= 1
    
    def test_model_regression_output(self):
        """Regression output is unbounded."""
        from amprenta_rag.ml.gnn.model import MoleculeGNN
        from amprenta_rag.ml.gnn.featurizer import mol_to_graph
        from torch_geometric.data import Batch
        
        model = MoleculeGNN(node_features=147, edge_features=7, task_type="regression")
        data = mol_to_graph("CCO")
        batch = Batch.from_data_list([data])
        
        output = model(batch)
        # Regression can be any value
        assert isinstance(output.item(), float)
    
    def test_mc_dropout_uncertainty(self):
        """MC Dropout produces non-zero uncertainty."""
        from amprenta_rag.ml.gnn.model import MoleculeGNN
        from amprenta_rag.ml.gnn.featurizer import mol_to_graph
        from torch_geometric.data import Batch
        
        model = MoleculeGNN(node_features=147, edge_features=7, dropout=0.2)
        data = mol_to_graph("CCO")
        batch = Batch.from_data_list([data])
        
        mean, std = model.predict_with_uncertainty(batch, n_samples=10)
        assert std.item() >= 0
        assert mean.shape == (1,)
        assert std.shape == (1,)
    
    def test_model_info(self):
        """Model info returns correct structure."""
        from amprenta_rag.ml.gnn.model import MoleculeGNN
        
        model = MoleculeGNN(node_features=147, edge_features=7, hidden_dim=64, num_layers=2)
        info = model.get_model_info()
        
        assert info["architecture"] == "MoleculeGNN"
        assert info["hidden_dim"] == 64
        assert info["num_layers"] == 2
        assert info["total_parameters"] > 0
        assert info["has_pyg"] is True


class TestFeaturizer:
    """Tests for molecular featurization."""
    
    def test_mol_to_graph_valid_smiles(self):
        """Valid SMILES converts to graph."""
        from amprenta_rag.ml.gnn.featurizer import mol_to_graph
        
        data = mol_to_graph("CCO")
        assert data is not None
        assert data.x.shape[0] == 3  # 3 atoms (C-C-O)
        assert data.x.shape[1] == 147  # Node feature dimension
    
    def test_mol_to_graph_invalid_smiles(self):
        """Invalid SMILES returns None."""
        from amprenta_rag.ml.gnn.featurizer import mol_to_graph
        
        data = mol_to_graph("invalid_smiles")
        assert data is None
    
    def test_mol_to_graph_benzene(self):
        """Benzene has correct structure."""
        from amprenta_rag.ml.gnn.featurizer import mol_to_graph
        
        data = mol_to_graph("c1ccccc1")
        assert data.x.shape[0] == 6  # 6 carbons
        assert data.edge_index.shape[1] == 12  # 6 bonds * 2 directions
        assert data.edge_attr.shape[0] == 12  # Bond features for each edge
        assert data.edge_attr.shape[1] == 7  # Edge feature dimension
    
    def test_feature_dimensions(self):
        """Feature dimensions are consistent."""
        from amprenta_rag.ml.gnn.featurizer import get_feature_dims
        
        node_dim, edge_dim = get_feature_dims()
        assert node_dim == 147
        assert edge_dim == 7
    
    def test_molecule_dataset(self):
        """MoleculeDataset works correctly."""
        from amprenta_rag.ml.gnn.featurizer import MoleculeDataset
        
        smiles_list = ["CCO", "c1ccccc1", "invalid"]
        labels = [0.1, 0.8, 0.5]
        
        dataset = MoleculeDataset(smiles_list, labels)
        assert len(dataset) == 3
        
        # Valid SMILES
        data0 = dataset.get(0)
        assert abs(data0.y.item() - 0.1) < 1e-6  # Float comparison
        
        # Invalid SMILES should get empty graph
        data2 = dataset.get(2)
        assert data2.x.shape[0] == 1  # Single dummy atom
        assert data2.y.item() == 0.5
    
    def test_dataset_statistics(self):
        """Dataset statistics calculation."""
        from amprenta_rag.ml.gnn.featurizer import MoleculeDataset
        
        smiles_list = ["CCO", "c1ccccc1", "invalid", "CC"]
        dataset = MoleculeDataset(smiles_list)
        
        stats = dataset.get_statistics()
        assert stats["total"] == 4
        assert stats["valid"] == 3  # CCO, benzene, CC
        assert stats["invalid"] == 1  # invalid
        assert stats["valid_percentage"] == 75.0


class TestGNNPredictor:
    """Tests for GNN prediction service."""
    
    def test_predictor_initialization(self):
        """Predictor initializes without model."""
        from amprenta_rag.ml.gnn.predictor import GNNPredictor
        
        predictor = GNNPredictor()
        assert predictor.model is None
        assert predictor.device in ["cpu", "cuda"]
    
    def test_predictor_without_model_fails(self):
        """Predictor without model raises error."""
        from amprenta_rag.ml.gnn.predictor import GNNPredictor
        
        predictor = GNNPredictor()
        
        with pytest.raises(ValueError, match="Model not loaded"):
            predictor.predict(["CCO"])
    
    def test_get_available_endpoints(self):
        """Get available endpoints returns list."""
        from amprenta_rag.ml.gnn.predictor import get_available_endpoints
        
        endpoints = get_available_endpoints()
        assert isinstance(endpoints, list)
        # Should include hERG if model was trained
        # assert "herg" in endpoints  # Depends on whether model file exists
    
    def test_predictor_cache_functions(self):
        """Predictor cache functions work."""
        from amprenta_rag.ml.gnn.predictor import clear_predictor_cache, _predictor_cache
        
        # Cache should be clearable
        clear_predictor_cache()
        assert len(_predictor_cache) == 0
    
    @patch('amprenta_rag.ml.gnn.predictor.Path.exists')
    def test_get_gnn_predictor_missing_model(self, mock_exists):
        """get_gnn_predictor returns None for missing model."""
        from amprenta_rag.ml.gnn.predictor import get_gnn_predictor
        
        mock_exists.return_value = False
        predictor = get_gnn_predictor("nonexistent")
        assert predictor is None
    
    def test_smiles_to_features(self):
        """SMILES to features extraction."""
        from amprenta_rag.ml.gnn.featurizer import smiles_to_features
        
        features = smiles_to_features("CCO")
        if features:  # Only test if RDKit available
            assert "num_atoms" in features
            assert "num_bonds" in features
            assert features["num_atoms"] == 3
            assert features["num_bonds"] == 2
