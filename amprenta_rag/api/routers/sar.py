from __future__ import annotations

from typing import List

from fastapi import APIRouter, HTTPException, Query, Depends
from pydantic import BaseModel, Field
from sqlalchemy.orm import Session

from amprenta_rag.api import schemas
from amprenta_rag.api.services import sar_data as service
from amprenta_rag.analysis.sar_whatif import TRANSFORMATIONS, compare_properties, scaffold_hop, validate_smiles
from amprenta_rag.chemistry.rgroup import find_common_core, decompose_rgroups
from amprenta_rag.chemistry.sar_analysis import build_sar_matrix
from amprenta_rag.database.base import get_db
from amprenta_rag.database.models import Compound

router = APIRouter()


@router.get(
    "/targets",
    summary="List SAR targets",
    description="Retrieve a list of biological targets with associated compound count for SAR analysis.",
    response_model=List[schemas.TargetResponse],
)
def list_targets(limit: int = Query(200, ge=1, le=2000, description="Maximum number of targets to return")):
    """Get list of SAR targets with compound counts.
    
    Returns biological targets that have associated compounds and activity data
    available for structure-activity relationship analysis.
    
    Args:
        limit: Maximum number of targets to return (1-2000)
        
    Returns:
        List of targets with their names and compound counts
    """
    return service.list_targets(limit=limit)


@router.get(
    "/targets/{target}/compounds",
    summary="List compounds with activity for target",
    description="Retrieve compounds and their biological activity data for a specific target protein.",
    response_model=List[schemas.CompoundActivityResponse],
)
def get_compounds_by_target(
    target: str = Field(..., description="Target protein name or identifier"),
    limit: int = Query(2000, ge=1, le=20000, description="Maximum number of compounds to return")
):
    """Get compounds with activity data for a specific target.
    
    Returns compounds that have been tested against the specified biological target,
    including their chemical structure (SMILES) and activity measurements (IC50, etc.).
    
    Args:
        target: Target protein name or identifier
        limit: Maximum number of compounds to return (1-20000)
        
    Returns:
        List of compounds with their structures and activity measurements
    """
    return service.get_compounds_by_target(target=target, limit=limit)


@router.get(
    "/targets/{target}/cliffs",
    summary="Detect activity cliffs for target",
    description="Find pairs of structurally similar compounds with large differences in biological activity.",
    response_model=List[schemas.ActivityCliffResponse],
)
def get_activity_cliffs(
    target: str = Field(..., description="Target protein name or identifier"),
    similarity_threshold: float = Query(0.6, ge=0.0, le=1.0, description="Minimum Tanimoto similarity threshold"),
    fold_change: float = Query(10.0, ge=1.0, le=1e6, description="Minimum fold change in activity"),
    limit: int = Query(50, ge=1, le=500, description="Maximum number of cliffs to return"),
):
    """Detect activity cliffs for a specific target.
    
    Activity cliffs are pairs of compounds that are structurally very similar
    but have dramatically different biological activity. These are valuable for
    understanding structure-activity relationships and optimizing compounds.
    
    Args:
        target: Target protein name or identifier
        similarity_threshold: Minimum Tanimoto similarity (0.0-1.0)
        fold_change: Minimum fold change in activity (≥1.0)
        limit: Maximum number of cliffs to return (1-500)
        
    Returns:
        List of activity cliff pairs with similarity and activity data
    """
    return service.get_activity_cliffs_for_target(
        target=target,
        similarity_threshold=similarity_threshold,
        fold_change=fold_change,
        limit=limit,
    )


class ValidateSmilesRequest(BaseModel):
    smiles: str = Field(..., min_length=1)


class PredictRequest(BaseModel):
    smiles_list: List[str]


class ScaffoldHopRequest(BaseModel):
    smiles: str = Field(..., min_length=1)
    transformation: str = Field(..., min_length=1)


@router.post(
    "/validate",
    summary="Validate SMILES string",
    description="Check if a SMILES string represents a valid chemical structure."
)
def sar_validate(payload: ValidateSmilesRequest) -> dict:
    """Validate a SMILES string for chemical correctness.
    
    Uses RDKit to parse and validate the SMILES string to ensure it represents
    a valid chemical structure that can be used in SAR analysis.
    
    Args:
        payload: Request containing the SMILES string to validate
        
    Returns:
        Validation result with the SMILES string and validity status
        
    Raises:
        HTTPException: 503 if RDKit is not available
    """
    try:
        ok = validate_smiles(payload.smiles)
        return {"smiles": payload.smiles, "valid": bool(ok)}
    except ImportError as e:
        raise HTTPException(status_code=503, detail=str(e))


@router.post(
    "/predict",
    summary="Predict molecular properties",
    description="Calculate molecular descriptors and properties for a list of SMILES strings."
)
def sar_predict(payload: PredictRequest) -> List[dict]:
    """Predict molecular properties for multiple compounds.
    
    Calculates various molecular descriptors including molecular weight, LogP,
    hydrogen bond donors/acceptors, and other drug-like properties for SAR analysis.
    
    Args:
        payload: Request containing list of SMILES strings to analyze
        
    Returns:
        List of dictionaries containing SMILES and calculated properties
        
    Raises:
        HTTPException: 503 if RDKit is not available, 400 for invalid input
    """
    try:
        df = compare_properties(payload.smiles_list)
        return df.to_dict(orient="records")
    except ImportError as e:
        raise HTTPException(status_code=503, detail=str(e))
    except Exception as e:  # noqa: BLE001
        raise HTTPException(status_code=400, detail=str(e))


@router.get(
    "/transformations",
    summary="List chemical transformations",
    description="Get available chemical transformations for scaffold hopping and lead optimization."
)
def list_transformations() -> List[dict]:
    """Get list of available chemical transformations.
    
    Returns predefined chemical transformations that can be applied to compounds
    for scaffold hopping, lead optimization, and exploring chemical space.
    
    Returns:
        List of transformation objects with IDs, names, and descriptions
    """
    out: List[dict] = []
    for k, v in TRANSFORMATIONS.items():
        out.append({"id": k, **v})
    return out


@router.post(
    "/scaffold-hop",
    summary="Perform scaffold hopping",
    description="Apply chemical transformations to generate new compounds with different scaffolds."
)
def sar_scaffold_hop(payload: ScaffoldHopRequest) -> dict:
    """Perform scaffold hopping on a compound.
    
    Applies a specified chemical transformation to generate new compounds with
    different core scaffolds while maintaining similar properties and activity.
    Useful for lead optimization and exploring chemical space.
    
    Args:
        payload: Request containing SMILES string and transformation ID
        
    Returns:
        Original SMILES, transformation applied, and list of product SMILES
        
    Raises:
        HTTPException: 503 if RDKit is not available, 400 for invalid input
    """
    try:
        products = scaffold_hop(payload.smiles, payload.transformation)
        return {"smiles": payload.smiles, "transformation": payload.transformation, "products": products}
    except ImportError as e:
        raise HTTPException(status_code=503, detail=str(e))
    except Exception as e:  # noqa: BLE001
        raise HTTPException(status_code=400, detail=str(e))


@router.get(
    "/scaffolds",
    summary="List molecular scaffolds",
    description="Get Murcko scaffolds extracted from compounds with frequency counts.",
    response_model=List[schemas.ScaffoldSummary],
)
def get_scaffolds(db: Session = Depends(get_db)) -> List[schemas.ScaffoldSummary]:
    """Get Murcko scaffolds and their compound frequencies.
    
    Extracts Murcko scaffolds (core ring systems) from all compounds in the database
    and returns them ranked by frequency. Useful for understanding scaffold diversity
    and identifying privileged structures for SAR analysis.
    
    Returns:
        List of scaffolds with SMILES strings and compound counts, sorted by frequency
        
    Raises:
        HTTPException: 503 if RDKit is not available
    """
    try:
        from rdkit import Chem
        from rdkit.Chem import Scaffolds
    except ImportError:
        raise HTTPException(status_code=503, detail="RDKit not available")
    
    # Get all compounds with valid SMILES
    compounds = db.query(Compound).filter(Compound.smiles.isnot(None)).all()
    
    scaffold_counts = {}
    for compound in compounds:
        try:
            mol = Chem.MolFromSmiles(compound.smiles)
            if mol:
                scaffold = Scaffolds.MurckoScaffold.GetScaffoldForMol(mol)
                scaffold_smiles = Chem.MolToSmiles(scaffold)
                scaffold_counts[scaffold_smiles] = scaffold_counts.get(scaffold_smiles, 0) + 1
        except Exception:
            continue  # Skip invalid molecules
    
    # Convert to response format
    scaffolds = [
        schemas.ScaffoldSummary(scaffold_smiles=smiles, compound_count=count)
        for smiles, count in sorted(scaffold_counts.items(), key=lambda x: x[1], reverse=True)
    ]
    
    return scaffolds


@router.post(
    "/grid",
    summary="Generate SAR matrix grid",
    description="Create a structure-activity relationship matrix from R-group decomposition of compounds.",
    response_model=schemas.SARGridResponse,
)
def build_sar_grid(
    request: schemas.SARGridRequest,
    db: Session = Depends(get_db)
) -> schemas.SARGridResponse:
    """Build SAR matrix grid from R-group decomposition.
    
    Takes a set of compounds, decomposes them into R-groups around a common core,
    and creates a 2D matrix showing activity relationships across different
    R-group substitution patterns. Essential for medicinal chemistry SAR analysis.
    
    Args:
        request: Grid configuration including compound IDs, core SMARTS, and axis definitions
        db: Database session for compound lookup
        
    Returns:
        SAR matrix with activity values, row/column labels, and metadata
        
    Raises:
        HTTPException: 404 if compounds not found, 400 for analysis failures
    """
    
    # Fetch compounds by IDs
    compounds = db.query(Compound).filter(
        Compound.compound_id.in_(request.compound_ids)
    ).all()
    
    if not compounds:
        raise HTTPException(status_code=404, detail="No compounds found")
    
    # Extract SMILES
    smiles_list = [c.smiles for c in compounds if c.smiles]
    compound_data = {c.smiles: c for c in compounds if c.smiles}
    
    if not smiles_list:
        raise HTTPException(status_code=400, detail="No valid SMILES found")
    
    try:
        # Find common core if not provided
        core_smiles = request.core_smarts
        if not core_smiles:
            core_smiles = find_common_core(smiles_list)
            if not core_smiles:
                raise HTTPException(status_code=400, detail="Could not find common core")
        
        # Decompose R-groups
        decomposed = decompose_rgroups(smiles_list, core_smiles)
        
        if not decomposed:
            raise HTTPException(status_code=400, detail="R-group decomposition failed")
        
        # Add activity data (using molecular weight as placeholder for now)
        for item in decomposed:
            smiles = item.get("smiles", "")
            if smiles in compound_data:
                compound = compound_data[smiles]
                item["ic50"] = compound.molecular_weight or 0.0  # Placeholder activity
        
        # Build SAR matrix
        matrix_df = build_sar_matrix(
            decomposed,
            x_axis=request.x_axis,
            y_axis=request.y_axis,
            activity_col="ic50"
        )
        
        if matrix_df.empty:
            raise HTTPException(status_code=400, detail="Could not build SAR matrix")
        
        # Convert to response format
        matrix = matrix_df.values.tolist()
        row_labels = matrix_df.index.tolist()
        col_labels = matrix_df.columns.tolist()
        
        return schemas.SARGridResponse(
            matrix=matrix,
            row_labels=[str(label) for label in row_labels],
            col_labels=[str(label) for label in col_labels],
            core_smiles=core_smiles,
            total_compounds=len(decomposed)
        )
        
    except Exception as e:
        raise HTTPException(status_code=400, detail=f"SAR grid generation failed: {str(e)}")


