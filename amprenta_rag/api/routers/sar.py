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
    response_model=List[schemas.TargetResponse],
)
def list_targets(limit: int = Query(200, ge=1, le=2000)):
    return service.list_targets(limit=limit)


@router.get(
    "/targets/{target}/compounds",
    summary="List compounds (with activity) for a target",
    response_model=List[schemas.CompoundActivityResponse],
)
def get_compounds_by_target(target: str, limit: int = Query(2000, ge=1, le=20000)):
    return service.get_compounds_by_target(target=target, limit=limit)


@router.get(
    "/targets/{target}/cliffs",
    summary="Detect activity cliffs for a target",
    response_model=List[schemas.ActivityCliffResponse],
)
def get_activity_cliffs(
    target: str,
    similarity_threshold: float = Query(0.6, ge=0.0, le=1.0),
    fold_change: float = Query(10.0, ge=1.0, le=1e6),
    limit: int = Query(50, ge=1, le=500),
):
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


@router.post("/validate")
def sar_validate(payload: ValidateSmilesRequest) -> dict:
    try:
        ok = validate_smiles(payload.smiles)
        return {"smiles": payload.smiles, "valid": bool(ok)}
    except ImportError as e:
        raise HTTPException(status_code=503, detail=str(e))


@router.post("/predict")
def sar_predict(payload: PredictRequest) -> List[dict]:
    try:
        df = compare_properties(payload.smiles_list)
        return df.to_dict(orient="records")
    except ImportError as e:
        raise HTTPException(status_code=503, detail=str(e))
    except Exception as e:  # noqa: BLE001
        raise HTTPException(status_code=400, detail=str(e))


@router.get("/transformations")
def list_transformations() -> List[dict]:
    out: List[dict] = []
    for k, v in TRANSFORMATIONS.items():
        out.append({"id": k, **v})
    return out


@router.post("/scaffold-hop")
def sar_scaffold_hop(payload: ScaffoldHopRequest) -> dict:
    try:
        products = scaffold_hop(payload.smiles, payload.transformation)
        return {"smiles": payload.smiles, "transformation": payload.transformation, "products": products}
    except ImportError as e:
        raise HTTPException(status_code=503, detail=str(e))
    except Exception as e:  # noqa: BLE001
        raise HTTPException(status_code=400, detail=str(e))


@router.get(
    "/scaffolds",
    summary="List scaffold summaries",
    response_model=List[schemas.ScaffoldSummary],
)
def get_scaffolds(db: Session = Depends(get_db)) -> List[schemas.ScaffoldSummary]:
    """Get Murcko scaffolds and compound counts."""
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
    summary="Build SAR matrix grid",
    response_model=schemas.SARGridResponse,
)
def build_sar_grid(
    request: schemas.SARGridRequest,
    db: Session = Depends(get_db)
) -> schemas.SARGridResponse:
    """Build SAR matrix from R-group decomposition."""
    
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


