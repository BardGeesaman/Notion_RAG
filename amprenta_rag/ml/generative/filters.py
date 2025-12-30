"""Molecular filters for generative chemistry (NoveltyChecker, DiversityFilter)."""

from __future__ import annotations

from typing import List, Optional, Set, Tuple

from amprenta_rag.logging_utils import get_logger

logger = get_logger(__name__)


def _require_rdkit():
    """Import RDKit with helpful error message."""
    try:
        from rdkit import Chem
        from rdkit.Chem import AllChem, DataStructs
        return Chem, AllChem, DataStructs
    except ImportError as e:
        raise ImportError("RDKit required for molecular filters") from e


class NoveltyChecker:
    """Check if molecules are novel vs reference set using Tanimoto similarity."""
    
    def __init__(self, reference_smiles: List[str], similarity_threshold: float = 0.9):
        """Initialize novelty checker.
        
        Args:
            reference_smiles: Known molecules (training set, existing compounds)
            similarity_threshold: Max Tanimoto similarity to be considered "novel"
        """
        Chem, AllChem, DataStructs = _require_rdkit()
        
        self.threshold = similarity_threshold
        
        # Precompute fingerprints for reference set
        self.reference_fps = []
        valid_count = 0
        
        for smiles in reference_smiles:
            fp = self._get_fingerprint(smiles)
            if fp is not None:
                self.reference_fps.append(fp)
                valid_count += 1
        
        logger.info(f"NoveltyChecker initialized: {valid_count}/{len(reference_smiles)} valid reference molecules")
        
        if not self.reference_fps:
            logger.warning("No valid reference molecules - all molecules will be considered novel")
    
    def _get_fingerprint(self, smiles: str):
        """Get Morgan fingerprint for SMILES string.
        
        Args:
            smiles: SMILES string
            
        Returns:
            RDKit fingerprint or None if invalid
        """
        Chem, AllChem, DataStructs = _require_rdkit()
        
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            return None
        
        return AllChem.GetMorganFingerprintAsBitVect(mol, radius=2, nBits=2048)
    
    def get_max_similarity(self, smiles: str) -> float:
        """Get maximum Tanimoto similarity to reference set.
        
        Args:
            smiles: Query molecule SMILES
            
        Returns:
            Maximum Tanimoto similarity (0.0 if invalid SMILES or no references)
        """
        Chem, AllChem, DataStructs = _require_rdkit()
        
        if not self.reference_fps:
            return 0.0
        
        query_fp = self._get_fingerprint(smiles)
        if query_fp is None:
            return 0.0
        
        max_sim = 0.0
        for ref_fp in self.reference_fps:
            sim = DataStructs.TanimotoSimilarity(query_fp, ref_fp)
            max_sim = max(max_sim, sim)
        
        return max_sim
    
    def is_novel(self, smiles: str) -> bool:
        """Check if molecule is novel (not similar to any reference).
        
        Args:
            smiles: Query molecule SMILES
            
        Returns:
            True if molecule is novel (max similarity < threshold)
        """
        max_sim = self.get_max_similarity(smiles)
        return max_sim < self.threshold
    
    def filter_novel(self, smiles_list: List[str]) -> List[str]:
        """Filter list to return only novel molecules.
        
        Args:
            smiles_list: List of SMILES to filter
            
        Returns:
            List of novel SMILES
        """
        novel_molecules = []
        
        for smiles in smiles_list:
            if self.is_novel(smiles):
                novel_molecules.append(smiles)
        
        logger.debug(f"NoveltyChecker: {len(novel_molecules)}/{len(smiles_list)} molecules are novel")
        return novel_molecules
    
    def get_similarity_report(self, smiles_list: List[str]) -> List[Tuple[str, float, bool]]:
        """Get detailed similarity report for molecules.
        
        Args:
            smiles_list: List of SMILES to analyze
            
        Returns:
            List of (smiles, max_similarity, is_novel) tuples
        """
        report = []
        
        for smiles in smiles_list:
            max_sim = self.get_max_similarity(smiles)
            is_novel = max_sim < self.threshold
            report.append((smiles, max_sim, is_novel))
        
        return report


class DiversityFilter:
    """Ensure diversity in generated molecule set using greedy selection."""
    
    def __init__(self, similarity_threshold: float = 0.7):
        """Initialize diversity filter.
        
        Args:
            similarity_threshold: Max Tanimoto similarity between selected molecules
        """
        self.threshold = similarity_threshold
        logger.info(f"DiversityFilter initialized: similarity_threshold={similarity_threshold}")
    
    def _get_fingerprint(self, smiles: str):
        """Get Morgan fingerprint for SMILES string."""
        Chem, AllChem, DataStructs = _require_rdkit()
        
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            return None
        
        return AllChem.GetMorganFingerprintAsBitVect(mol, radius=2, nBits=2048)
    
    def _get_similarity(self, smiles1: str, smiles2: str) -> float:
        """Get Tanimoto similarity between two molecules.
        
        Args:
            smiles1: First molecule SMILES
            smiles2: Second molecule SMILES
            
        Returns:
            Tanimoto similarity (0.0 if either SMILES is invalid)
        """
        Chem, AllChem, DataStructs = _require_rdkit()
        
        fp1 = self._get_fingerprint(smiles1)
        fp2 = self._get_fingerprint(smiles2)
        
        if fp1 is None or fp2 is None:
            return 0.0
        
        return DataStructs.TanimotoSimilarity(fp1, fp2)
    
    def filter_diverse(self, smiles_list: List[str], max_count: int = 10) -> List[str]:
        """Select diverse subset using greedy algorithm.
        
        Greedy selection: keep molecule if Tanimoto < threshold to all already-selected.
        
        Args:
            smiles_list: Input molecules to filter
            max_count: Maximum number of molecules to select
            
        Returns:
            List of diverse SMILES (up to max_count)
        """
        if not smiles_list:
            return []
        
        # Filter valid SMILES first
        valid_smiles = []
        for smiles in smiles_list:
            if self._get_fingerprint(smiles) is not None:
                valid_smiles.append(smiles)
        
        if not valid_smiles:
            logger.warning("No valid SMILES found for diversity filtering")
            return []
        
        # Greedy selection
        selected = []
        
        for smiles in valid_smiles:
            if len(selected) >= max_count:
                break
            
            # Check if diverse from all already selected
            is_diverse = True
            for selected_smiles in selected:
                similarity = self._get_similarity(smiles, selected_smiles)
                if similarity >= self.threshold:
                    is_diverse = False
                    break
            
            if is_diverse:
                selected.append(smiles)
        
        logger.debug(f"DiversityFilter: selected {len(selected)}/{len(valid_smiles)} diverse molecules")
        return selected
    
    def get_diversity_matrix(self, smiles_list: List[str]) -> List[List[float]]:
        """Compute pairwise similarity matrix.
        
        Args:
            smiles_list: List of SMILES
            
        Returns:
            N x N similarity matrix
        """
        n = len(smiles_list)
        matrix = [[0.0] * n for _ in range(n)]
        
        for i in range(n):
            for j in range(i, n):
                if i == j:
                    similarity = 1.0
                else:
                    similarity = self._get_similarity(smiles_list[i], smiles_list[j])
                
                matrix[i][j] = similarity
                matrix[j][i] = similarity  # Symmetric
        
        return matrix
    
    def compute_set_diversity(self, smiles_list: List[str]) -> float:
        """Compute average pairwise diversity (1 - similarity) for a set.
        
        Args:
            smiles_list: List of SMILES
            
        Returns:
            Average diversity score (0-1, higher is more diverse)
        """
        if len(smiles_list) < 2:
            return 1.0  # Single molecule is maximally diverse
        
        matrix = self.get_diversity_matrix(smiles_list)
        n = len(smiles_list)
        
        total_similarity = 0.0
        count = 0
        
        for i in range(n):
            for j in range(i + 1, n):
                total_similarity += matrix[i][j]
                count += 1
        
        avg_similarity = total_similarity / count if count > 0 else 0.0
        return 1.0 - avg_similarity  # Convert similarity to diversity
    
    def filter_diverse_with_scores(self, smiles_list: List[str], max_count: int = 10) -> List[Tuple[str, float]]:
        """Select diverse molecules with diversity scores.
        
        Args:
            smiles_list: Input molecules
            max_count: Maximum number to select
            
        Returns:
            List of (smiles, diversity_score) tuples
        """
        selected_smiles = self.filter_diverse(smiles_list, max_count)
        
        results = []
        for smiles in selected_smiles:
            # Compute diversity score as min distance to other selected molecules
            min_similarity = 1.0
            for other_smiles in selected_smiles:
                if other_smiles != smiles:
                    similarity = self._get_similarity(smiles, other_smiles)
                    min_similarity = min(min_similarity, similarity)
            
            diversity_score = 1.0 - min_similarity
            results.append((smiles, diversity_score))
        
        return results
