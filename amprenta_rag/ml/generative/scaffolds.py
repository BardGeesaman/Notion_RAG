"""Murcko scaffold extraction for scaffold hopping and analysis."""

from __future__ import annotations

from collections import defaultdict
from typing import Any, Dict, List, Optional, Set, Tuple, Union

from amprenta_rag.logging_utils import get_logger

logger = get_logger(__name__)


def _require_rdkit():
    """Import RDKit with helpful error message."""
    try:
        from rdkit import Chem
        from rdkit.Chem.Scaffolds import MurckoScaffold
        return Chem, MurckoScaffold
    except ImportError as e:
        raise ImportError("RDKit required for scaffold extraction") from e


class ScaffoldExtractor:
    """Extract Murcko scaffolds for scaffold hopping and molecular analysis."""
    
    @staticmethod
    def get_scaffold(smiles: str, generic: bool = False) -> Optional[str]:
        """Extract Murcko scaffold from SMILES.
        
        Args:
            smiles: Input molecule SMILES
            generic: If True, return generic scaffold (atoms->C, bonds->single)
            
        Returns:
            Scaffold SMILES or None if invalid input
        """
        Chem, MurckoScaffold = _require_rdkit()
        
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            return None
        
        try:
            # Extract Murcko scaffold
            core = MurckoScaffold.GetScaffoldForMol(mol)
            
            if core is None:
                return None
            
            # Make generic if requested
            if generic:
                core = MurckoScaffold.MakeScaffoldGeneric(core)
            
            return Chem.MolToSmiles(core)
        
        except Exception as e:
            logger.warning(f"Failed to extract scaffold from {smiles}: {e}")
            return None
    
    @staticmethod
    def get_framework(smiles: str) -> Optional[str]:
        """Extract molecular framework (scaffold + linkers).
        
        Args:
            smiles: Input molecule SMILES
            
        Returns:
            Framework SMILES or None if invalid
        """
        Chem, MurckoScaffold = _require_rdkit()
        
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            return None
        
        try:
            framework = MurckoScaffold.GetFrameworkForMol(mol)
            if framework is None:
                return None
            
            return Chem.MolToSmiles(framework)
        
        except Exception as e:
            logger.warning(f"Failed to extract framework from {smiles}: {e}")
            return None
    
    @staticmethod
    def get_sidechains(smiles: str) -> List[str]:
        """Extract R-groups/sidechains from molecule.
        
        This is a simplified approach that removes the scaffold and 
        returns the remaining fragments.
        
        Args:
            smiles: Input molecule SMILES
            
        Returns:
            List of sidechain SMILES
        """
        Chem, MurckoScaffold = _require_rdkit()
        
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            return []
        
        try:
            # Get scaffold
            core = MurckoScaffold.GetScaffoldForMol(mol)
            if core is None:
                return []
            
            # This is a simplified approach - in practice, R-group decomposition
            # is more complex and would require specialized algorithms
            # For now, we'll return an empty list as a placeholder
            
            # NOTE: R-group decomposition enhancement tracked in ROADMAP
            # Could use RDKit's RGroupDecomposition or similar approaches
            
            return []
        
        except Exception as e:
            logger.warning(f"Failed to extract sidechains from {smiles}: {e}")
            return []
    
    @staticmethod
    def group_by_scaffold(smiles_list: List[str], generic: bool = False) -> Dict[str, List[str]]:
        """Group molecules by their Murcko scaffold.
        
        Args:
            smiles_list: List of molecule SMILES
            generic: Use generic scaffolds
            
        Returns:
            Dictionary mapping scaffold SMILES -> list of molecule SMILES
        """
        scaffold_groups = defaultdict(list)
        
        for smiles in smiles_list:
            scaffold = ScaffoldExtractor.get_scaffold(smiles, generic=generic)
            
            if scaffold is not None:
                scaffold_groups[scaffold].append(smiles)
            else:
                # Group invalid molecules separately
                scaffold_groups["<invalid>"].append(smiles)
        
        # Convert to regular dict
        return dict(scaffold_groups)
    
    @staticmethod
    def get_scaffold_statistics(smiles_list: List[str]) -> Dict[str, Union[int, float, List[Tuple[str, int]]]]:
        """Get statistics about scaffold diversity in a molecule set.
        
        Args:
            smiles_list: List of molecule SMILES
            
        Returns:
            Dictionary with scaffold statistics
        """
        scaffold_groups = ScaffoldExtractor.group_by_scaffold(smiles_list)
        
        total_molecules = len(smiles_list)
        unique_scaffolds = len([k for k in scaffold_groups.keys() if k != "<invalid>"])
        invalid_molecules = len(scaffold_groups.get("<invalid>", []))
        
        # Find most common scaffolds
        scaffold_counts = {k: len(v) for k, v in scaffold_groups.items() if k != "<invalid>"}
        sorted_scaffolds = sorted(scaffold_counts.items(), key=lambda x: x[1], reverse=True)
        
        return {
            "total_molecules": total_molecules,
            "unique_scaffolds": unique_scaffolds,
            "invalid_molecules": invalid_molecules,
            "valid_molecules": total_molecules - invalid_molecules,
            "scaffold_diversity": unique_scaffolds / max(1, total_molecules - invalid_molecules),
            "most_common_scaffolds": sorted_scaffolds[:10],  # Top 10
        }
    
    @staticmethod
    def find_scaffold_analogs(
        query_smiles: str, 
        database_smiles: List[str], 
        max_analogs: int = 20
    ) -> List[Tuple[str, str]]:
        """Find molecules with the same scaffold as query.
        
        Args:
            query_smiles: Query molecule SMILES
            database_smiles: Database of molecules to search
            max_analogs: Maximum number of analogs to return
            
        Returns:
            List of (analog_smiles, scaffold_smiles) tuples
        """
        query_scaffold = ScaffoldExtractor.get_scaffold(query_smiles)
        if query_scaffold is None:
            logger.warning(f"Could not extract scaffold from query: {query_smiles}")
            return []
        
        analogs = []
        
        for smiles in database_smiles:
            if smiles == query_smiles:
                continue  # Skip the query molecule itself
            
            scaffold = ScaffoldExtractor.get_scaffold(smiles)
            if scaffold == query_scaffold:
                analogs.append((smiles, scaffold))
                
                if len(analogs) >= max_analogs:
                    break
        
        logger.debug(f"Found {len(analogs)} scaffold analogs for {query_smiles}")
        return analogs
    
    @staticmethod
    def scaffold_similarity(smiles1: str, smiles2: str) -> float:
        """Compute scaffold similarity between two molecules.
        
        This is a simple binary similarity: 1.0 if same scaffold, 0.0 otherwise.
        More sophisticated scaffold similarities could be implemented using
        graph-based methods.
        
        Args:
            smiles1: First molecule SMILES
            smiles2: Second molecule SMILES
            
        Returns:
            Scaffold similarity (0.0 or 1.0)
        """
        scaffold1 = ScaffoldExtractor.get_scaffold(smiles1)
        scaffold2 = ScaffoldExtractor.get_scaffold(smiles2)
        
        if scaffold1 is None or scaffold2 is None:
            return 0.0
        
        return 1.0 if scaffold1 == scaffold2 else 0.0
    
    @staticmethod
    def get_ring_systems(smiles: str) -> List[str]:
        """Extract individual ring systems from molecule.
        
        Args:
            smiles: Input molecule SMILES
            
        Returns:
            List of ring system SMILES
        """
        Chem, MurckoScaffold = _require_rdkit()
        
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            return []
        
        try:
            # Get ring info
            ring_info = mol.GetRingInfo()
            ring_systems = []
            
            # Group rings into systems (connected rings)
            ring_groups: List[Set[int]] = []
            for ring in ring_info.AtomRings():
                ring_set = set(ring)
                
                # Find which group this ring belongs to
                merged = False
                for i, group in enumerate(ring_groups):
                    if ring_set & group:  # Overlapping rings
                        ring_groups[i] = group | ring_set
                        merged = True
                        break
                
                if not merged:
                    ring_groups.append(ring_set)
            
            # Extract SMILES for each ring system
            for ring_group in ring_groups:
                # Create submolecule with only ring atoms
                atom_indices = list(ring_group)
                if len(atom_indices) >= 3:  # Minimum for a ring
                    try:
                        # This is a simplified extraction - proper implementation
                        # would need to handle bond connectivity carefully
                        ring_systems.append(f"ring_system_{len(atom_indices)}_atoms")
                    except Exception:
                        continue
            
            return ring_systems
        
        except Exception as e:
            logger.warning(f"Failed to extract ring systems from {smiles}: {e}")
            return []
