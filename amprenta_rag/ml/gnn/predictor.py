"""GNN Predictor for toxicity endpoint inference."""
import torch
from pathlib import Path
from typing import List, Dict, Any, Optional
import logging

logger = logging.getLogger(__name__)

try:
    from torch_geometric.data import Batch
    HAS_PYG = True
except ImportError:
    HAS_PYG = False
    Batch = None


class GNNPredictor:
    """Production predictor for GNN toxicity models.
    
    Supports:
    - Loading models from file or MLModelRegistry
    - Batch inference with uncertainty (MC Dropout)
    - Error handling for invalid SMILES
    """
    
    def __init__(self, model_path: Optional[str] = None, device: str = "auto"):
        if not HAS_PYG:
            raise ImportError("torch_geometric required for GNNPredictor")
        
        self.device = "cuda" if device == "auto" and torch.cuda.is_available() else "cpu"
        self.model = None
        self.config = None
        self.endpoint = None
        
        if model_path:
            self.load(model_path)
    
    def load(self, model_path: str) -> None:
        """Load model from checkpoint."""
        from amprenta_rag.ml.gnn.model import MoleculeGNN
        
        try:
            checkpoint = torch.load(model_path, map_location=self.device)
            self.config = checkpoint["config"]
            self.endpoint = checkpoint.get("endpoint")
            
            self.model = MoleculeGNN(
                node_features=self.config["node_features"],
                edge_features=self.config["edge_features"],
                task_type=self.config["task_type"]
            ).to(self.device)
            
            self.model.load_state_dict(checkpoint["state_dict"])
            self.model.eval()
            
            logger.info(f"Loaded GNN model for {self.endpoint} on {self.device}")
        except Exception as e:
            logger.error(f"Failed to load GNN model from {model_path}: {e}")
            raise
    
    def predict(
        self,
        smiles_list: List[str],
        with_uncertainty: bool = True,
        n_samples: int = 10
    ) -> List[Dict[str, Any]]:
        """Predict toxicity for list of SMILES.
        
        Args:
            smiles_list: List of SMILES strings
            with_uncertainty: Use MC Dropout for uncertainty
            n_samples: Number of MC Dropout samples
            
        Returns:
            List of prediction dicts with keys:
            - smiles: input SMILES
            - prediction: model output (probability for classification)
            - uncertainty: std from MC Dropout (if with_uncertainty)
            - error: error message if prediction failed
        """
        if self.model is None:
            raise ValueError("Model not loaded")
        
        from amprenta_rag.ml.gnn.featurizer import mol_to_graph
        
        results = []
        
        # Process each SMILES
        valid_indices = []
        data_list = []
        
        for i, smiles in enumerate(smiles_list):
            data = mol_to_graph(smiles)
            if data is None:
                results.append({
                    "smiles": smiles,
                    "error": "Failed to parse SMILES"
                })
            else:
                valid_indices.append(i)
                data_list.append(data)
                results.append({"smiles": smiles})
        
        if not data_list:
            return results
        
        # Batch inference
        batch = Batch.from_data_list(data_list).to(self.device)
        
        try:
            if with_uncertainty:
                mean, std = self.model.predict_with_uncertainty(batch, n_samples)
                predictions = mean.cpu().numpy()
                uncertainties = std.cpu().numpy()
            else:
                self.model.eval()
                with torch.no_grad():
                    predictions = self.model(batch).cpu().numpy()
                uncertainties = [None] * len(predictions)
            
            # Fill in results
            for idx, (pred, unc) in zip(valid_indices, zip(predictions, uncertainties)):
                results[idx]["prediction"] = float(pred)
                if unc is not None:
                    results[idx]["uncertainty"] = float(unc)
                
                # Add interpretation for classification
                if self.config["task_type"] == "classification":
                    results[idx]["class"] = 1 if pred > 0.5 else 0
                    results[idx]["label"] = "toxic" if pred > 0.5 else "non-toxic"
                    results[idx]["confidence"] = abs(pred - 0.5) * 2  # Distance from boundary
        
        except Exception as e:
            logger.error(f"GNN prediction error: {e}")
            for idx in valid_indices:
                results[idx]["error"] = f"Prediction failed: {str(e)}"
        
        return results
    
    def get_model_info(self) -> Dict[str, Any]:
        """Get model information."""
        if self.model is None:
            return {"status": "not_loaded"}
        
        info = self.model.get_model_info()
        info.update({
            "endpoint": self.endpoint,
            "device": self.device,
            "config": self.config
        })
        return info


# Global predictor cache
_predictor_cache: Dict[str, GNNPredictor] = {}


def get_gnn_predictor(endpoint: str) -> Optional[GNNPredictor]:
    """Get cached GNN predictor for endpoint."""
    if endpoint not in _predictor_cache:
        model_path = Path(f"models/gnn/gnn_{endpoint}.pt")
        if not model_path.exists():
            logger.warning(f"GNN model not found: {model_path}")
            return None
        
        try:
            _predictor_cache[endpoint] = GNNPredictor(str(model_path))
        except Exception as e:
            logger.error(f"Failed to load GNN predictor for {endpoint}: {e}")
            return None
    
    return _predictor_cache[endpoint]


def clear_predictor_cache() -> None:
    """Clear the predictor cache."""
    global _predictor_cache
    _predictor_cache.clear()
    logger.info("GNN predictor cache cleared")


def get_available_endpoints() -> List[str]:
    """Get list of endpoints with trained models."""
    available = []
    for endpoint in ["herg", "ames", "dili", "ld50", "clintox"]:
        model_path = Path(f"models/gnn/gnn_{endpoint}.pt")
        if model_path.exists():
            available.append(endpoint)
    return available