"""R-group decomposition utilities for SAR analysis."""
from __future__ import annotations

from typing import List, Dict, Optional

from amprenta_rag.logging_utils import get_logger

logger = get_logger(__name__)

try:
    from rdkit import Chem, RDLogger
    from rdkit.Chem import rdFMCS
    RDLogger.DisableLog("rdApp.*")
    RDKIT_AVAILABLE = True
except ImportError:
    RDKIT_AVAILABLE = False
    Chem = None  # type: ignore
    rdFMCS = None  # type: ignore


def _mol_from_smiles_best_effort(smiles: str):
    """Parse SMILES with fallback for kekulization issues."""
    if not smiles:
        return None
    m = Chem.MolFromSmiles(smiles)
    if m is not None:
        return m
    try:
        m = Chem.MolFromSmiles(smiles, sanitize=False)
        if m is None:
            return None
        Chem.SanitizeMol(
            m,
            sanitizeOps=Chem.SanitizeFlags.SANITIZE_ALL ^ Chem.SanitizeFlags.SANITIZE_KEKULIZE,
        )
        return m
    except Exception:
        return None


def find_common_core(smiles_list: List[str]) -> Optional[str]:
    """
    Find maximum common substructure (MCS) from a list of SMILES.

    Args:
        smiles_list: List of SMILES strings

    Returns:
        SMARTS string of the common core, or None if MCS fails
    """
    if not RDKIT_AVAILABLE:
        logger.error("[RGROUP] RDKit not available")
        return None

    if not smiles_list:
        logger.warning("[RGROUP] Empty SMILES list")
        return None

    try:
        # Convert SMILES to molecules
        mols = []
        for smi in smiles_list:
            mol = _mol_from_smiles_best_effort(smi)
            if mol is None:
                logger.debug("[RGROUP] Invalid SMILES: %s", smi)
                continue
            mols.append(mol)

        if len(mols) < 2:
            logger.warning("[RGROUP] Need at least 2 valid molecules")
            return None

        # Find MCS
        logger.info("[RGROUP] Finding MCS for %d molecules", len(mols))
        mcs = rdFMCS.FindMCS(mols, timeout=10)

        if mcs.numAtoms == 0:
            logger.warning("[RGROUP] No common substructure found")
            return None

        core_smarts = mcs.smartsString
        logger.info("[RGROUP] Found common core with %d atoms: %s", mcs.numAtoms, core_smarts[:50])
        return core_smarts

    except Exception as e:
        logger.error("[RGROUP] Error finding common core: %r", e)
        return None


def decompose_rgroups(smiles_list: List[str], core_smarts: str) -> List[Dict[str, any]]:
    """
    Decompose compounds into R-groups based on a common core.

    Args:
        smiles_list: List of SMILES strings
        core_smarts: SMARTS pattern for the common core

    Returns:
        List of dicts with keys: smiles, R1, R2, ... (R-group positions)
    """
    if not RDKIT_AVAILABLE:
        logger.error("[RGROUP] RDKit not available")
        return []

    if not core_smarts:
        logger.warning("[RGROUP] No core SMARTS provided")
        return []

    try:
        # Create core pattern
        core_pattern = Chem.MolFromSmarts(core_smarts)
        if core_pattern is None:
            logger.error("[RGROUP] Invalid core SMARTS: %s", core_smarts)
            return []

        results = []

        for smi in smiles_list:
            mol = _mol_from_smiles_best_effort(smi)
            if mol is None:
                logger.debug("[RGROUP] Invalid SMILES: %s", smi)
                continue

            # Check if molecule matches core
            if not mol.HasSubstructMatch(core_pattern):
                logger.debug("[RGROUP] Molecule %s doesn't match core", smi)
                results.append({"smiles": smi, "error": "No match"})
                continue

            # Get matching atoms
            matches = mol.GetSubstructMatches(core_pattern)
            if not matches:
                results.append({"smiles": smi, "error": "No match"})
                continue

            # Use first match
            match = matches[0]
            match_set = set(match)

            # Find R-groups by identifying atoms not in the core match
            # This is a simplified approach - full R-group decomposition is more complex
            rgroups = {}
            rgroup_idx = 1

            # Get atoms and bonds
            atoms = [mol.GetAtomWithIdx(i) for i in range(mol.GetNumAtoms())]

            # Find attachment points (atoms in core that have bonds to non-core atoms)
            attachment_points = []
            for atom_idx in match:
                atom = mol.GetAtomWithIdx(atom_idx)
                for neighbor in atom.GetNeighbors():
                    if neighbor.GetIdx() not in match_set:
                        attachment_points.append((atom_idx, neighbor.GetIdx()))

            # Extract R-groups from attachment points
            for core_atom_idx, r_atom_idx in attachment_points:
                # Extract the R-group fragment
                # This is simplified - full implementation would need to handle multiple R-groups per position
                rgroup_mol = Chem.PathToSubmol(mol, [r_atom_idx], useQuery=True)
                if rgroup_mol:
                    rgroup_smiles = Chem.MolToSmiles(rgroup_mol)
                    rgroups[f"R{rgroup_idx}"] = rgroup_smiles
                    rgroup_idx += 1

            result = {"smiles": smi, **rgroups}
            results.append(result)

        logger.info("[RGROUP] Decomposed %d compounds", len(results))
        return results

    except Exception as e:
        logger.error("[RGROUP] Error decomposing R-groups: %r", e)
        return []


def get_rgroup_statistics(decomposition: List[Dict[str, any]]) -> Dict[str, Dict[str, int]]:
    """
    Calculate frequency statistics for R-groups at each position.

    Args:
        decomposition: List of R-group decomposition dicts from decompose_rgroups

    Returns:
        Dict mapping R-group position (R1, R2, ...) to dict of {smiles: count}
    """
    if not decomposition:
        logger.warning("[RGROUP] Empty decomposition list")
        return {}

    try:
        stats = {}

        # Collect all R-group positions
        all_positions = set()
        for item in decomposition:
            for key in item.keys():
                if key.startswith("R") and key[1:].isdigit():
                    all_positions.add(key)

        # Count frequencies for each position
        for position in sorted(all_positions):
            counts: Dict[str, int] = {}
            for item in decomposition:
                if position in item and item[position]:
                    rgroup_smiles = item[position]
                    counts[rgroup_smiles] = counts.get(rgroup_smiles, 0) + 1
            stats[position] = counts

        logger.info("[RGROUP] Calculated statistics for %d R-group positions", len(stats))
        return stats

    except Exception as e:
        logger.error("[RGROUP] Error calculating statistics: %r", e)
        return {}
