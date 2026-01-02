"""Graph Neural Network architecture for molecules."""
import torch
import torch.nn as nn
import torch.nn.functional as F
from typing import Optional
import logging

logger = logging.getLogger(__name__)

try:
    from torch_geometric.nn import GINEConv, global_mean_pool
    from torch_geometric.data import Data
    HAS_PYG = True
except ImportError:
    HAS_PYG = False
    Data = None


class MoleculeGNN(nn.Module):
    """Message Passing Neural Network for molecular property prediction.
    
    Architecture:
    - Node encoder (atom features → hidden)
    - N GINEConv layers with batch norm
    - Global mean pooling
    - MLP prediction head
    
    Args:
        node_features: Number of atom features (default 78)
        edge_features: Number of bond features (default 10)
        hidden_dim: Hidden layer dimension (default 128)
        num_layers: Number of message passing layers (default 3)
        dropout: Dropout rate for MC Dropout uncertainty (default 0.2)
        task_type: "classification" or "regression"
    """
    
    def __init__(
        self,
        node_features: int = 78,
        edge_features: int = 10,
        hidden_dim: int = 128,
        num_layers: int = 3,
        dropout: float = 0.2,
        task_type: str = "classification"
    ):
        super().__init__()
        
        if not HAS_PYG:
            raise ImportError("torch_geometric required for MoleculeGNN")
        
        self.task_type = task_type
        self.dropout = dropout
        self.num_layers = num_layers
        self.hidden_dim = hidden_dim
        
        # Node embedding
        self.node_encoder = nn.Linear(node_features, hidden_dim)
        
        # Message passing layers
        self.conv_layers = nn.ModuleList()
        self.batch_norms = nn.ModuleList()
        
        for _ in range(num_layers):
            mlp = nn.Sequential(
                nn.Linear(hidden_dim, hidden_dim),
                nn.ReLU(),
                nn.Linear(hidden_dim, hidden_dim)
            )
            self.conv_layers.append(GINEConv(mlp, edge_dim=edge_features))
            self.batch_norms.append(nn.BatchNorm1d(hidden_dim))
        
        # Prediction head
        self.head = nn.Sequential(
            nn.Linear(hidden_dim, hidden_dim // 2),
            nn.ReLU(),
            nn.Dropout(dropout),
            nn.Linear(hidden_dim // 2, 1)
        )
    
    def forward(self, data) -> torch.Tensor:
        """Forward pass.
        
        Args:
            data: PyG Data object with x, edge_index, edge_attr, batch
            
        Returns:
            Predictions tensor (batch_size,)
        """
        x = data.x
        edge_index = data.edge_index
        edge_attr = data.edge_attr
        batch = data.batch
        
        # Encode nodes
        x = self.node_encoder(x)
        
        # Message passing
        for conv, bn in zip(self.conv_layers, self.batch_norms):
            x = conv(x, edge_index, edge_attr)
            x = bn(x)
            x = F.relu(x)
            x = F.dropout(x, p=self.dropout, training=self.training)
        
        # Graph-level readout
        x = global_mean_pool(x, batch)
        
        # Prediction
        out = self.head(x).squeeze(-1)
        
        if self.task_type == "classification":
            out = torch.sigmoid(out)
        
        return out
    
    def predict_with_uncertainty(
        self, 
        data, 
        n_samples: int = 10
    ) -> tuple[torch.Tensor, torch.Tensor]:
        """Predict with MC Dropout uncertainty estimation.
        
        Args:
            data: PyG Data object
            n_samples: Number of forward passes for MC Dropout
            
        Returns:
            (mean_prediction, std_prediction) tensors
        """
        self.train()  # Enable dropout
        
        predictions = []
        with torch.no_grad():
            for _ in range(n_samples):
                pred = self.forward(data)
                predictions.append(pred)
        
        predictions = torch.stack(predictions)
        mean = predictions.mean(dim=0)
        std = predictions.std(dim=0)
        
        self.eval()
        return mean, std
    
    def get_model_info(self) -> dict:
        """Get model architecture information."""
        total_params = sum(p.numel() for p in self.parameters())
        trainable_params = sum(p.numel() for p in self.parameters() if p.requires_grad)
        
        return {
            "architecture": "MoleculeGNN",
            "task_type": self.task_type,
            "hidden_dim": self.hidden_dim,
            "num_layers": self.num_layers,
            "dropout": self.dropout,
            "total_parameters": total_params,
            "trainable_parameters": trainable_params,
            "has_pyg": HAS_PYG
        }
