#!/usr/bin/env python
"""Train GNN models for ADMET toxicity endpoints.

Usage:
    python scripts/train_gnn_models.py --endpoint herg --epochs 50
    python scripts/train_gnn_models.py --endpoint all --dry-run
"""
import argparse
import logging
from pathlib import Path
from typing import Dict, Any, Optional
from datetime import datetime
import json
import os
import sys

# Add project root to path
sys.path.insert(0, str(Path(__file__).parent.parent))

import torch
import torch.nn as nn
from torch.optim import Adam
from torch.optim.lr_scheduler import ReduceLROnPlateau

logger = logging.getLogger(__name__)

# TDC endpoint configurations
GNN_ENDPOINTS = {
    "herg": {
        "tdc_class": "Tox",
        "tdc_name": "hERG",
        "type": "classification",
        "description": "hERG channel inhibition (cardiotoxicity)"
    },
    "ames": {
        "tdc_class": "Tox", 
        "tdc_name": "AMES",
        "type": "classification",
        "description": "Mutagenicity (Ames test)"
    },
    "dili": {
        "tdc_class": "Tox",
        "tdc_name": "DILI",
        "type": "classification", 
        "description": "Drug-induced liver injury"
    },
    "ld50": {
        "tdc_class": "Tox",
        "tdc_name": "LD50_Zhu",
        "type": "regression",
        "description": "Acute oral toxicity (LD50)"
    },
    "clintox": {
        "tdc_class": "Tox",
        "tdc_name": "ClinTox",
        "type": "classification",
        "description": "Clinical trial toxicity"
    }
}


def load_tdc_dataset(endpoint: str) -> tuple:
    """Load TDC dataset for endpoint."""
    from tdc.single_pred import Tox, ADME
    
    config = GNN_ENDPOINTS[endpoint]
    tdc_class = Tox if config["tdc_class"] == "Tox" else ADME
    dataset = tdc_class(name=config["tdc_name"])
    splits = dataset.get_split()
    return splits['train'], splits['valid'], splits['test']


def create_dataloaders(train_df, val_df, test_df, batch_size: int = 32):
    """Create PyG DataLoaders from dataframes."""
    from torch_geometric.loader import DataLoader
    from amprenta_rag.ml.gnn.featurizer import MoleculeDataset
    
    train_dataset = MoleculeDataset(
        train_df["Drug"].tolist(),
        train_df["Y"].tolist()
    )
    val_dataset = MoleculeDataset(
        val_df["Drug"].tolist(),
        val_df["Y"].tolist()
    )
    test_dataset = MoleculeDataset(
        test_df["Drug"].tolist(),
        test_df["Y"].tolist()
    )
    
    train_loader = DataLoader(train_dataset, batch_size=batch_size, shuffle=True)
    val_loader = DataLoader(val_dataset, batch_size=batch_size)
    test_loader = DataLoader(test_dataset, batch_size=batch_size)
    
    return train_loader, val_loader, test_loader


def train_epoch(model, loader, optimizer, criterion, device):
    """Train for one epoch."""
    model.train()
    total_loss = 0
    
    for batch in loader:
        batch = batch.to(device)
        optimizer.zero_grad()
        
        out = model(batch)
        loss = criterion(out, batch.y.float())
        
        loss.backward()
        optimizer.step()
        total_loss += loss.item() * batch.num_graphs
    
    return total_loss / len(loader.dataset)


def evaluate(model, loader, criterion, device, is_classification: bool):
    """Evaluate model on loader."""
    model.eval()
    total_loss = 0
    all_preds = []
    all_labels = []
    
    with torch.no_grad():
        for batch in loader:
            batch = batch.to(device)
            out = model(batch)
            loss = criterion(out, batch.y.float())
            total_loss += loss.item() * batch.num_graphs
            
            all_preds.extend(out.cpu().numpy())
            all_labels.extend(batch.y.cpu().numpy())
    
    avg_loss = total_loss / len(loader.dataset)
    
    # Compute metrics
    import numpy as np
    preds = np.array(all_preds)
    labels = np.array(all_labels)
    
    if is_classification:
        from sklearn.metrics import roc_auc_score, accuracy_score
        auc = roc_auc_score(labels, preds) if len(np.unique(labels)) > 1 else 0.5
        acc = accuracy_score(labels, (preds > 0.5).astype(int))
        return {"loss": avg_loss, "auc": auc, "accuracy": acc}
    else:
        from sklearn.metrics import mean_squared_error, r2_score
        rmse = np.sqrt(mean_squared_error(labels, preds))
        r2 = r2_score(labels, preds)
        return {"loss": avg_loss, "rmse": rmse, "r2": r2}


def train_endpoint(
    endpoint: str,
    epochs: int = 50,
    batch_size: int = 32,
    lr: float = 0.001,
    patience: int = 10,
    device: str = "auto",
    save_dir: Optional[str] = None,
    dry_run: bool = False
) -> Dict[str, Any]:
    """Train GNN model for a single endpoint."""
    from amprenta_rag.ml.gnn.model import MoleculeGNN
    from amprenta_rag.ml.gnn.featurizer import get_feature_dims
    
    config = GNN_ENDPOINTS[endpoint]
    is_classification = config["type"] == "classification"
    
    print(f"\n{'='*60}")
    print(f"Training {endpoint} ({config['type']})")
    print(f"Description: {config['description']}")
    print(f"{'='*60}")
    
    if dry_run:
        print("DRY RUN - would train model")
        return {"endpoint": endpoint, "dry_run": True}
    
    # Device
    if device == "auto":
        device = "cuda" if torch.cuda.is_available() else "cpu"
    print(f"Device: {device}")
    
    # Load data
    print("Loading TDC dataset...")
    train_df, val_df, test_df = load_tdc_dataset(endpoint)
    print(f"Splits: train={len(train_df)}, val={len(val_df)}, test={len(test_df)}")
    
    # Create loaders
    train_loader, val_loader, test_loader = create_dataloaders(
        train_df, val_df, test_df, batch_size
    )
    
    # Create model
    node_dim, edge_dim = get_feature_dims()
    model = MoleculeGNN(
        node_features=node_dim,
        edge_features=edge_dim,
        task_type=config["type"]
    ).to(device)
    
    print(f"Model params: {sum(p.numel() for p in model.parameters()):,}")
    
    # Training setup
    optimizer = Adam(model.parameters(), lr=lr)
    scheduler = ReduceLROnPlateau(optimizer, mode='min', patience=patience//2)
    
    if is_classification:
        criterion = nn.BCELoss()
    else:
        criterion = nn.MSELoss()
    
    # Training loop
    best_val_loss = float('inf')
    best_model_state = None
    no_improve = 0
    
    print(f"\nTraining for up to {epochs} epochs...")
    for epoch in range(epochs):
        train_loss = train_epoch(model, train_loader, optimizer, criterion, device)
        val_metrics = evaluate(model, val_loader, criterion, device, is_classification)
        
        scheduler.step(val_metrics['loss'])
        
        # Early stopping
        if val_metrics['loss'] < best_val_loss:
            best_val_loss = val_metrics['loss']
            best_model_state = model.state_dict().copy()
            no_improve = 0
        else:
            no_improve += 1
        
        if epoch % 5 == 0 or no_improve == 0:
            metric_str = f"AUC={val_metrics.get('auc', 0):.3f}" if is_classification else f"RMSE={val_metrics.get('rmse', 0):.3f}"
            print(f"Epoch {epoch:3d}: train_loss={train_loss:.4f}, val_loss={val_metrics['loss']:.4f}, {metric_str}")
        
        if no_improve >= patience:
            print(f"Early stopping at epoch {epoch}")
            break
    
    # Load best model and evaluate on test
    if best_model_state:
        model.load_state_dict(best_model_state)
    test_metrics = evaluate(model, test_loader, criterion, device, is_classification)
    
    print(f"\nTest Results:")
    for k, v in test_metrics.items():
        print(f"  {k}: {v:.4f}")
    
    # Save model
    if save_dir:
        save_path = Path(save_dir)
        save_path.mkdir(parents=True, exist_ok=True)
        
        model_path = save_path / f"gnn_{endpoint}.pt"
        torch.save({
            "state_dict": best_model_state or model.state_dict(),
            "config": {
                "node_features": node_dim,
                "edge_features": edge_dim,
                "task_type": config["type"],
            },
            "metrics": test_metrics,
            "endpoint": endpoint,
            "trained_at": datetime.now().isoformat()
        }, model_path)
        print(f"Model saved to {model_path}")
    
    return {
        "endpoint": endpoint,
        "test_metrics": test_metrics,
        "epochs_trained": epoch + 1,
        "best_val_loss": best_val_loss
    }


def main():
    parser = argparse.ArgumentParser(description="Train GNN toxicity models")
    parser.add_argument("--endpoint", default="herg", help="Endpoint or 'all'")
    parser.add_argument("--epochs", type=int, default=50)
    parser.add_argument("--batch-size", type=int, default=32)
    parser.add_argument("--lr", type=float, default=0.001)
    parser.add_argument("--patience", type=int, default=10)
    parser.add_argument("--device", default="auto")
    parser.add_argument("--save-dir", default="models/gnn")
    parser.add_argument("--dry-run", action="store_true")
    
    args = parser.parse_args()
    
    # Set OpenMP workaround
    os.environ["KMP_DUPLICATE_LIB_OK"] = "TRUE"
    
    endpoints = list(GNN_ENDPOINTS.keys()) if args.endpoint == "all" else [args.endpoint]
    
    results = []
    for ep in endpoints:
        if ep not in GNN_ENDPOINTS:
            print(f"Unknown endpoint: {ep}")
            continue
        
        result = train_endpoint(
            ep,
            epochs=args.epochs,
            batch_size=args.batch_size,
            lr=args.lr,
            patience=args.patience,
            device=args.device,
            save_dir=args.save_dir,
            dry_run=args.dry_run
        )
        results.append(result)
    
    print("\n" + "="*60)
    print("Training Summary")
    print("="*60)
    for r in results:
        if "test_metrics" in r:
            metrics = r["test_metrics"]
            print(f"{r['endpoint']}: {metrics}")
        else:
            print(f"{r['endpoint']}: {r}")


if __name__ == "__main__":
    main()
