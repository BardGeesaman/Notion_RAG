"""Molecular graph featurization for GNN input."""
import torch
from typing import Optional, List
import logging

logger = logging.getLogger(__name__)

try:
    from torch_geometric.data import Data, Dataset
    from rdkit import Chem
    from rdkit.Chem import rdchem, Descriptors
    HAS_DEPS = True
except ImportError:
    HAS_DEPS = False
    Data = None
    Dataset = None


# Atom feature dimensions
ATOM_FEATURES = {
    'atomic_num': list(range(1, 119)),  # All elements
    'degree': [0, 1, 2, 3, 4, 5, 6],
    'formal_charge': [-2, -1, 0, 1, 2],
    'hybridization': [
        rdchem.HybridizationType.SP,
        rdchem.HybridizationType.SP2,
        rdchem.HybridizationType.SP3,
        rdchem.HybridizationType.SP3D,
        rdchem.HybridizationType.SP3D2
    ] if HAS_DEPS else [],
    'num_hs': [0, 1, 2, 3, 4],
}

# Bond feature dimensions
BOND_FEATURES = {
    'bond_type': [
        rdchem.BondType.SINGLE,
        rdchem.BondType.DOUBLE,
        rdchem.BondType.TRIPLE,
        rdchem.BondType.AROMATIC
    ] if HAS_DEPS else [],
}


def one_hot(value, choices: list) -> List[int]:
    """One-hot encode a value given choices."""
    encoding = [0] * (len(choices) + 1)  # +1 for unknown
    try:
        idx = choices.index(value)
        encoding[idx] = 1
    except ValueError:
        encoding[-1] = 1  # Unknown
    return encoding


def mol_to_graph(smiles: str) -> Optional[Data]:
    """Convert SMILES to PyTorch Geometric Data object.
    
    Args:
        smiles: SMILES string
        
    Returns:
        PyG Data object or None if conversion fails
    """
    if not HAS_DEPS:
        raise ImportError("RDKit and torch_geometric required")
    
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        logger.warning(f"Failed to parse SMILES: {smiles[:50]}...")
        return None
    
    # Atom features
    atom_features = []
    for atom in mol.GetAtoms():
        features = (
            one_hot(atom.GetAtomicNum(), ATOM_FEATURES['atomic_num']) +
            one_hot(atom.GetDegree(), ATOM_FEATURES['degree']) +
            one_hot(atom.GetFormalCharge(), ATOM_FEATURES['formal_charge']) +
            one_hot(atom.GetHybridization(), ATOM_FEATURES['hybridization']) +
            one_hot(atom.GetTotalNumHs(), ATOM_FEATURES['num_hs']) +
            [1 if atom.GetIsAromatic() else 0] +
            [1 if atom.IsInRing() else 0]
        )
        atom_features.append(features)
    
    x = torch.tensor(atom_features, dtype=torch.float)
    
    # Edge features (bonds)
    edge_index = []
    edge_attr = []
    
    for bond in mol.GetBonds():
        i = bond.GetBeginAtomIdx()
        j = bond.GetEndAtomIdx()
        
        # Undirected: add both directions
        edge_index.extend([[i, j], [j, i]])
        
        bond_features = (
            one_hot(bond.GetBondType(), BOND_FEATURES['bond_type']) +
            [1 if bond.GetIsConjugated() else 0] +
            [1 if bond.IsInRing() else 0]
        )
        edge_attr.extend([bond_features, bond_features])
    
    if len(edge_index) == 0:
        # Single atom molecule
        edge_index = torch.zeros((2, 0), dtype=torch.long)
        edge_attr = torch.zeros((0, 6), dtype=torch.float)
    else:
        edge_index = torch.tensor(edge_index, dtype=torch.long).t().contiguous()
        edge_attr = torch.tensor(edge_attr, dtype=torch.float)
    
    return Data(x=x, edge_index=edge_index, edge_attr=edge_attr)


def get_feature_dimensions() -> dict:
    """Get the dimensions of atom and bond features."""
    if not HAS_DEPS:
        return {"node_features": 0, "edge_features": 0}
    
    # Calculate feature dimensions
    node_dim = (
        len(ATOM_FEATURES['atomic_num']) + 1 +  # atomic_num + unknown
        len(ATOM_FEATURES['degree']) + 1 +       # degree + unknown
        len(ATOM_FEATURES['formal_charge']) + 1 + # formal_charge + unknown
        len(ATOM_FEATURES['hybridization']) + 1 + # hybridization + unknown
        len(ATOM_FEATURES['num_hs']) + 1 +        # num_hs + unknown
        2  # is_aromatic, is_in_ring
    )
    
    edge_dim = (
        len(BOND_FEATURES['bond_type']) + 1 +  # bond_type + unknown
        2  # is_conjugated, is_in_ring
    )
    
    return {"node_features": node_dim, "edge_features": edge_dim}


def get_feature_dims() -> tuple[int, int]:
    """Get feature dimensions for model initialization.
    
    Returns:
        (node_features, edge_features) tuple
    """
    dims = get_feature_dimensions()
    return dims["node_features"], dims["edge_features"]


class MoleculeDataset(Dataset):
    """PyTorch Geometric Dataset for molecules."""
    
    def __init__(self, smiles_list: List[str], labels: Optional[List[float]] = None):
        super().__init__()
        if not HAS_DEPS:
            raise ImportError("torch_geometric required for MoleculeDataset")
        
        self.smiles_list = smiles_list
        self.labels = labels
        self._data_list = None
        
        # Get feature dimensions for empty graph fallback
        self.feature_dims = get_feature_dimensions()
    
    def len(self) -> int:
        return len(self.smiles_list)
    
    def get(self, idx: int) -> Data:
        if self._data_list is None:
            self._data_list = [None] * len(self.smiles_list)
        
        if self._data_list[idx] is None:
            data = mol_to_graph(self.smiles_list[idx])
            if data is None:
                # Return empty graph for invalid SMILES
                data = Data(
                    x=torch.zeros((1, self.feature_dims["node_features"]), dtype=torch.float),
                    edge_index=torch.zeros((2, 0), dtype=torch.long),
                    edge_attr=torch.zeros((0, self.feature_dims["edge_features"]), dtype=torch.float)
                )
            
            if self.labels is not None:
                data.y = torch.tensor([self.labels[idx]], dtype=torch.float)
            
            self._data_list[idx] = data
        
        return self._data_list[idx]
    
    def get_statistics(self) -> dict:
        """Get dataset statistics."""
        if not self.smiles_list:
            return {"total": 0, "valid": 0, "invalid": 0}
        
        valid_count = 0
        invalid_count = 0
        
        for smiles in self.smiles_list:
            if mol_to_graph(smiles) is not None:
                valid_count += 1
            else:
                invalid_count += 1
        
        return {
            "total": len(self.smiles_list),
            "valid": valid_count,
            "invalid": invalid_count,
            "valid_percentage": (valid_count / len(self.smiles_list)) * 100
        }


def smiles_to_features(smiles: str) -> Optional[dict]:
    """Extract basic molecular features from SMILES for analysis.
    
    Args:
        smiles: SMILES string
        
    Returns:
        Dictionary with molecular features or None if parsing fails
    """
    if not HAS_DEPS:
        return None
    
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return None
    
    return {
        "num_atoms": mol.GetNumAtoms(),
        "num_bonds": mol.GetNumBonds(),
        "num_rings": mol.GetRingInfo().NumRings(),
        "molecular_weight": Chem.Descriptors.ExactMolWt(mol),
        "num_aromatic_atoms": sum(1 for atom in mol.GetAtoms() if atom.GetIsAromatic()),
        "num_heteroatoms": sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() != 6)
    }
