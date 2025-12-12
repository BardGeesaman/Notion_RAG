"""Global search service for searching across multiple entity types."""
from __future__ import annotations

from typing import Dict, List, Any

from amprenta_rag.database.models import Experiment, Compound, Signature, Dataset
from amprenta_rag.logging_utils import get_logger

logger = get_logger(__name__)


def global_search(query: str, db, limit: int = 5) -> Dict[str, List[Dict[str, Any]]]:
    """
    Search across multiple entity types in the database.
    
    Args:
        query: Search query string
        db: Database session
        limit: Maximum number of results per entity type
        
    Returns:
        Dict with results grouped by entity type:
        {
            "experiments": [{id, name, type}],
            "compounds": [{id, compound_id, smiles}],
            "signatures": [{id, name}],
            "datasets": [{id, name}]
        }
    """
    if not query or not query.strip():
        return {
            "experiments": [],
            "compounds": [],
            "signatures": [],
            "datasets": [],
        }
    
    search_term = f"%{query.strip()}%"
    results = {
        "experiments": [],
        "compounds": [],
        "signatures": [],
        "datasets": [],
    }
    
    # Search Experiments
    experiments = (
        db.query(Experiment)
        .filter(
            (Experiment.name.ilike(search_term)) | (Experiment.description.ilike(search_term))
        )
        .limit(limit)
        .all()
    )
    for exp in experiments:
        results["experiments"].append({
            "id": str(exp.id),
            "name": exp.name,
            "type": exp.type,
        })
    
    # Search Compounds
    compounds = (
        db.query(Compound)
        .filter(
            (Compound.compound_id.ilike(search_term)) | (Compound.canonical_smiles.ilike(search_term))
        )
        .limit(limit)
        .all()
    )
    for comp in compounds:
        results["compounds"].append({
            "id": str(comp.id),
            "compound_id": comp.compound_id,
            "smiles": comp.smiles[:100] + "..." if comp.smiles and len(comp.smiles) > 100 else comp.smiles,
        })
    
    # Search Signatures
    signatures = (
        db.query(Signature)
        .filter(Signature.name.ilike(search_term))
        .limit(limit)
        .all()
    )
    for sig in signatures:
        results["signatures"].append({
            "id": str(sig.id),
            "name": sig.name,
        })
    
    # Search Datasets
    datasets = (
        db.query(Dataset)
        .filter(Dataset.name.ilike(search_term))
        .limit(limit)
        .all()
    )
    for ds in datasets:
        results["datasets"].append({
            "id": str(ds.id),
            "name": ds.name,
        })
    
    logger.debug("[GLOBAL_SEARCH] Query '%s' returned: %d experiments, %d compounds, %d signatures, %d datasets",
                 query, len(results["experiments"]), len(results["compounds"]), 
                 len(results["signatures"]), len(results["datasets"]))
    
    return results
