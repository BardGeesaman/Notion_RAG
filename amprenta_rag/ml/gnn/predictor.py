"""GNN predictor wrapper for production use."""
import torch
import logging
from typing import List, Optional, Dict, Any
from pathlib import Path

from amprenta_rag.ml.gnn.model import MoleculeGNN
from amprenta_rag.ml.gnn.featurizer import mol_to_graph, MoleculeDataset, get_feature_dimensions

logger = logging.getLogger(__name__)

try:
    from torch_geometric.data import Batch
    HAS_PYG = True
except ImportError:
    HAS_PYG = False
    Batch = None


class GNNPredictor:
    """Production wrapper for GNN molecular property prediction."""
    
    def __init__(
        self,
        model_path: Optional[str] = None,
        task_type: str = "classification",
        device: str = "cpu"
    ):
        """Initialize GNN predictor.
        
        Args:
            model_path: Path to saved model weights
            task_type: "classification" or "regression"
            device: "cpu" or "cuda"
        """
        if not HAS_PYG:
            raise ImportError("torch_geometric required for GNNPredictor")
        
        self.task_type = task_type
        self.device = device
        self.model = None
        self.feature_dims = get_feature_dimensions()
        
        if model_path and Path(model_path).exists():
            self.load_model(model_path)
        else:
            # Initialize with default architecture
            self.model = MoleculeGNN(
                node_features=self.feature_dims["node_features"],
                edge_features=self.feature_dims["edge_features"],
                task_type=task_type
            ).to(device)
    
    def load_model(self, model_path: str) -> None:
        """Load model weights from file."""
        try:
            checkpoint = torch.load(model_path, map_location=self.device)
            
            # Extract model config if available
            config = checkpoint.get("config", {})
            
            self.model = MoleculeGNN(
                node_features=config.get("node_features", self.feature_dims["node_features"]),
                edge_features=config.get("edge_features", self.feature_dims["edge_features"]),
                hidden_dim=config.get("hidden_dim", 128),
                num_layers=config.get("num_layers", 3),
                dropout=config.get("dropout", 0.2),
                task_type=config.get("task_type", self.task_type)
            ).to(self.device)
            
            self.model.load_state_dict(checkpoint["model_state_dict"])
            self.model.eval()
            
            logger.info(f"Loaded GNN model from {model_path}")
            
        except Exception as e:
            logger.error(f"Failed to load model from {model_path}: {e}")
            raise
    
    def save_model(self, model_path: str, metadata: Optional[dict] = None) -> None:
        """Save model weights and config to file."""
        if self.model is None:
            raise ValueError("No model to save")
        
        checkpoint = {
            "model_state_dict": self.model.state_dict(),
            "config": {
                "node_features": self.feature_dims["node_features"],
                "edge_features": self.feature_dims["edge_features"],
                "hidden_dim": self.model.hidden_dim,
                "num_layers": self.model.num_layers,
                "dropout": self.model.dropout,
                "task_type": self.model.task_type
            },
            "metadata": metadata or {}
        }
        
        torch.save(checkpoint, model_path)
        logger.info(f"Saved GNN model to {model_path}")
    
    def predict(self, smiles_list: List[str]) -> List[float]:
        """Predict properties for a list of SMILES.
        
        Args:
            smiles_list: List of SMILES strings
            
        Returns:
            List of predictions
        """
        if self.model is None:
            raise ValueError("Model not initialized")
        
        # Convert to graphs
        graphs = []
        for smiles in smiles_list:
            graph = mol_to_graph(smiles)
            if graph is not None:
                graphs.append(graph)
            else:
                # Add dummy graph for invalid SMILES
                graphs.append(Data(
                    x=torch.zeros((1, self.feature_dims["node_features"]), dtype=torch.float),
                    edge_index=torch.zeros((2, 0), dtype=torch.long),
                    edge_attr=torch.zeros((0, self.feature_dims["edge_features"]), dtype=torch.float)
                ))
        
        if not graphs:
            return []
        
        # Create batch
        batch = Batch.from_data_list(graphs).to(self.device)
        
        # Predict
        self.model.eval()
        with torch.no_grad():
            predictions = self.model(batch)
        
        return predictions.cpu().numpy().tolist()
    
    def predict_with_uncertainty(
        self, 
        smiles_list: List[str], 
        n_samples: int = 10
    ) -> List[Dict[str, float]]:
        """Predict with uncertainty estimation using MC Dropout.
        
        Args:
            smiles_list: List of SMILES strings
            n_samples: Number of MC Dropout samples
            
        Returns:
            List of dicts with "prediction" and "uncertainty" keys
        """
        if self.model is None:
            raise ValueError("Model not initialized")
        
        # Convert to graphs
        graphs = []
        for smiles in smiles_list:
            graph = mol_to_graph(smiles)
            if graph is not None:
                graphs.append(graph)
            else:
                # Add dummy graph for invalid SMILES
                graphs.append(Data(
                    x=torch.zeros((1, self.feature_dims["node_features"]), dtype=torch.float),
                    edge_index=torch.zeros((2, 0), dtype=torch.long),
                    edge_attr=torch.zeros((0, self.feature_dims["edge_features"]), dtype=torch.float)
                ))
        
        if not graphs:
            return []
        
        # Create batch
        batch = Batch.from_data_list(graphs).to(self.device)
        
        # MC Dropout prediction
        mean_pred, std_pred = self.model.predict_with_uncertainty(batch, n_samples)
        
        results = []
        for i in range(len(smiles_list)):
            results.append({
                "prediction": float(mean_pred[i]),
                "uncertainty": float(std_pred[i])
            })
        
        return results
    
    def get_model_info(self) -> dict:
        """Get model information."""
        if self.model is None:
            return {"status": "not_initialized"}
        
        info = self.model.get_model_info()
        info.update({
            "device": self.device,
            "feature_dimensions": self.feature_dims
        })
        return info
