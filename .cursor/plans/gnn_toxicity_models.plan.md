# Neural Network Toxicity Models (GNN)

## Overview
Graph Neural Networks for molecular property prediction, extending beyond fingerprint-based XGBoost models. Deep learning ADMET using message passing on molecular graphs.

**Timeline:** 4-5 weeks (includes buffer for PyG installation)
**Dependencies:** PyTorch, PyTorch Geometric, TDC datasets
**GPU:** Recommended (2-4 hrs/model on GPU, 20+ hrs on CPU)

**Reviewer Status:** ✅ APPROVED WITH P1 FIXES (applied below)

### Key Decisions
- **Uncertainty:** MC Dropout (N=10 forward passes). Ensemble deferred to P2.
- **Hyperparameters:** Use TDC defaults. Tuning deferred to P2.
- **Tox21 multi-label:** Deferred to P3 (too complex for MVP)

---

## Architecture

```
┌─────────────────┐    ┌──────────────────┐    ┌─────────────────┐
│  SMILES Input   │───▶│  Mol → Graph     │───▶│  GNN Encoder    │
│                 │    │  (atoms, bonds)  │    │  (MPNN/GAT)     │
└─────────────────┘    └──────────────────┘    └─────────────────┘
                                                       │
                                                       ▼
┌─────────────────┐    ┌──────────────────┐    ┌─────────────────┐
│  Prediction     │◀───│  Readout Layer   │◀───│  Graph Pooling  │
│  (toxicity)     │    │  (MLP head)      │    │  (mean/attention)│
└─────────────────┘    └──────────────────┘    └─────────────────┘
```

---

## GNN Model Design

### Message Passing Neural Network (MPNN)

```python
class MoleculeGNN(nn.Module):
    """Graph Neural Network for molecular property prediction."""
    
    def __init__(
        self,
        node_features: int = 78,  # Atom features
        edge_features: int = 10,  # Bond features
        hidden_dim: int = 128,
        num_layers: int = 3,
        dropout: float = 0.2,
        task_type: str = "classification"  # or "regression"
    ):
        super().__init__()
        
        # Node embedding
        self.node_encoder = nn.Linear(node_features, hidden_dim)
        
        # Message passing layers
        self.conv_layers = nn.ModuleList([
            GINEConv(
                nn.Sequential(
                    nn.Linear(hidden_dim, hidden_dim),
                    nn.ReLU(),
                    nn.Linear(hidden_dim, hidden_dim)
                ),
                edge_dim=edge_features
            )
            for _ in range(num_layers)
        ])
        
        # Batch normalization
        self.batch_norms = nn.ModuleList([
            nn.BatchNorm1d(hidden_dim) for _ in range(num_layers)
        ])
        
        # Readout
        self.pool = global_mean_pool
        
        # Prediction head
        self.head = nn.Sequential(
            nn.Linear(hidden_dim, hidden_dim // 2),
            nn.ReLU(),
            nn.Dropout(dropout),
            nn.Linear(hidden_dim // 2, 1)
        )
        
        self.task_type = task_type
    
    def forward(self, data):
        x, edge_index, edge_attr, batch = (
            data.x, data.edge_index, data.edge_attr, data.batch
        )
        
        # Encode nodes
        x = self.node_encoder(x)
        
        # Message passing
        for conv, bn in zip(self.conv_layers, self.batch_norms):
            x = conv(x, edge_index, edge_attr)
            x = bn(x)
            x = F.relu(x)
        
        # Graph-level readout
        x = self.pool(x, batch)
        
        # Prediction
        out = self.head(x)
        
        if self.task_type == "classification":
            out = torch.sigmoid(out)
        
        return out.squeeze(-1)
```

### Molecular Graph Featurization

```python
def mol_to_graph(smiles: str) -> Data:
    """Convert SMILES to PyTorch Geometric Data object."""
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return None
    
    # Atom features (78 dims)
    atom_features = []
    for atom in mol.GetAtoms():
        features = [
            # Atom type one-hot (common elements)
            *one_hot(atom.GetAtomicNum(), [6, 7, 8, 9, 15, 16, 17, 35, 53]),
            # Degree
            *one_hot(atom.GetDegree(), [0, 1, 2, 3, 4, 5]),
            # Formal charge
            atom.GetFormalCharge(),
            # Hybridization
            *one_hot(atom.GetHybridization(), [SP, SP2, SP3, SP3D, SP3D2]),
            # Aromaticity
            atom.GetIsAromatic(),
            # Ring membership
            atom.IsInRing(),
            # Hydrogen count
            atom.GetTotalNumHs(),
        ]
        atom_features.append(features)
    
    x = torch.tensor(atom_features, dtype=torch.float)
    
    # Edge features (bond features)
    edge_index = []
    edge_attr = []
    for bond in mol.GetBonds():
        i, j = bond.GetBeginAtomIdx(), bond.GetEndAtomIdx()
        edge_index.extend([[i, j], [j, i]])
        
        bond_features = [
            *one_hot(bond.GetBondType(), [SINGLE, DOUBLE, TRIPLE, AROMATIC]),
            bond.GetIsConjugated(),
            bond.IsInRing(),
        ]
        edge_attr.extend([bond_features, bond_features])
    
    edge_index = torch.tensor(edge_index, dtype=torch.long).t().contiguous()
    edge_attr = torch.tensor(edge_attr, dtype=torch.float)
    
    return Data(x=x, edge_index=edge_index, edge_attr=edge_attr)
```

---

## Endpoints & Models

### Toxicity Endpoints (TDC)

| Endpoint | TDC Dataset | Type | Description |
|----------|-------------|------|-------------|
| herg_gnn | hERG | Classification | hERG channel inhibition (cardiotoxicity) |
| ames_gnn | AMES | Classification | Mutagenicity (Ames test) |
| dili_gnn | DILI | Classification | Drug-induced liver injury |
| ld50_gnn | LD50_Zhu | Regression | Acute oral toxicity |
| clintox_gnn | ClinTox | Classification | Clinical trial toxicity |

**Deferred to P3:** Tox21 multi-label (12 tasks - too complex for MVP)

### Integration with Existing ADMET

```python
# In amprenta_rag/ml/admet/predictor.py

GNN_MODELS = {
    "herg_gnn": "admet_herg_gnn",
    "ames_gnn": "admet_ames_gnn",
    "dili_gnn": "admet_dili_gnn",
    "ld50_gnn": "admet_ld50_gnn",
    "clintox_gnn": "admet_clintox_gnn",
}

class ADMETPredictor:
    def predict_gnn(
        self,
        smiles_list: List[str],
        endpoints: Optional[List[str]] = None
    ) -> List[Dict[str, Any]]:
        """Predict using GNN models."""
        # Load GNN model from registry
        # Convert SMILES to graphs
        # Batch inference
        # Return predictions with uncertainty
```

---

## Service Layer

### `amprenta_rag/ml/gnn/model.py`
- MoleculeGNN class
- Graph featurization utilities
- Model loading/saving

### `amprenta_rag/ml/gnn/trainer.py`
- Training loop with early stopping
- TDC dataset loading
- Metrics computation (AUC, RMSE)
- Model checkpointing

### `amprenta_rag/ml/gnn/predictor.py`
- GNNPredictor class
- Batch inference with error handling (return error dict for invalid SMILES)
- MC Dropout uncertainty (N=10 forward passes)
- Model loading from MLModelRegistry (model_type="gnn_classification" or "gnn_regression")

---

## API Endpoints

Add to existing `amprenta_rag/api/routers/admet.py`:

| Method | Endpoint | Description |
|--------|----------|-------------|
| POST | /admet/predict-gnn | GNN predictions for SMILES list |
| GET | /admet/gnn-models | List available GNN models |
| GET | /admet/gnn-models/{name}/info | Model metadata and metrics |
| POST | /admet/compare | Compare XGBoost vs GNN predictions |

---

## Dashboard

### New Tab: "Deep Learning ADMET"

Add to existing ADMET dashboard or create new page:

1. **Prediction Interface**
   - SMILES input (single or batch)
   - Endpoint selection (GNN models)
   - Side-by-side: XGBoost vs GNN predictions
   - Uncertainty visualization

2. **Model Comparison**
   - Performance metrics table (AUC, RMSE per endpoint)
   - Calibration plots
   - Applicability domain comparison

3. **Attention Visualization** (future)
   - Highlight important atoms/bonds
   - Interpretability via attention weights

---

## Implementation Batches

### Batch 1: Core GNN Infrastructure (Days 1-4)
- [ ] Create `amprenta_rag/ml/gnn/` module
- [ ] Implement MoleculeGNN model class
- [ ] Implement mol_to_graph featurization
- [ ] Add PyTorch Geometric to requirements

### Batch 2: Training Pipeline (Days 5-8)
- [ ] Create `scripts/train_gnn_models.py`
- [ ] TDC dataset integration
- [ ] Training loop with validation
- [ ] Model checkpointing to registry

### Batch 3: Predictor + API (Days 9-12)
- [ ] GNNPredictor class with uncertainty
- [ ] API endpoints (4)
- [ ] Integration with ADMETPredictor
- [ ] Pydantic schemas

### Batch 4: Dashboard (Days 13-15)
- [ ] New "Deep Learning ADMET" tab
- [ ] Prediction interface
- [ ] Model comparison visualization

### Batch 5: Tests + Documentation (Days 16-18)
- [ ] Service tests (12)
- [ ] API tests (10)
- [ ] E2E tests (4)
- [ ] ROADMAP update

---

## Success Criteria

- [ ] MoleculeGNN training on at least 3 TDC endpoints
- [ ] GNN models registered in MLModelRegistry
- [ ] API predictions working via `/admet/predict-gnn`
- [ ] Dashboard shows GNN vs XGBoost comparison
- [ ] 26+ tests passing
- [ ] Documentation complete

---

## Dependencies

Add to `requirements.txt`:
```
torch>=2.0.0
torch-geometric>=2.4.0
torch-scatter
torch-sparse
```

Or via conda:
```bash
conda install pytorch pytorch-cuda=12.1 -c pytorch -c nvidia
conda install pyg -c pyg
```

---

## TODOs

- [ ] Batch 1: Core GNN Infrastructure
- [ ] Batch 2: Training Pipeline
- [ ] Batch 3: Predictor + API
- [ ] Batch 4: Dashboard
- [ ] Batch 5: Tests + Documentation


