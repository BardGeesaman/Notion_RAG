"""Chemistry endpoints (ADMET prediction, etc.)."""

from __future__ import annotations

from typing import List, Optional
from uuid import UUID

from fastapi import APIRouter, HTTPException
from pydantic import BaseModel

from amprenta_rag.database.models import Compound
from amprenta_rag.database.session import db_session
from amprenta_rag.ml.admet import get_admet_predictor

router = APIRouter()


class ADMETPredictRequest(BaseModel):
    smiles_list: List[str]
    endpoints: Optional[List[str]] = None  # ["herg", "logs", "logp"]
    include_shap: bool = False


@router.post("/chemistry/predict-admet")
def predict_admet(request: ADMETPredictRequest):
    """Predict ADMET properties for compounds."""
    try:
        predictor = get_admet_predictor()
        return predictor.predict(
            smiles_list=request.smiles_list,
            endpoints=request.endpoints,
            include_shap=request.include_shap,
        )
    except ImportError as e:
        raise HTTPException(status_code=503, detail=str(e))
    except Exception as e:  # noqa: BLE001
        raise HTTPException(status_code=500, detail=f"ADMET prediction failed: {e}")


@router.post("/chemistry/compounds/{compound_id}/predict-admet")
def predict_admet_for_compound(
    compound_id: UUID,
    request: ADMETPredictRequest,
):
    """Predict ADMET properties for a specific compound by UUID (fetches SMILES from DB)."""
    try:
        with db_session() as db:
            compound = db.query(Compound).filter(Compound.id == compound_id).first()
            if compound is None:
                raise HTTPException(status_code=404, detail="Compound not found")
            smiles = getattr(compound, "smiles", None)
            if not smiles:
                raise HTTPException(status_code=422, detail="Compound missing SMILES")

        predictor = get_admet_predictor()
        out = predictor.predict(
            smiles_list=[smiles],
            endpoints=request.endpoints,
            include_shap=request.include_shap,
        )
        return out[0] if out else {"smiles": smiles, "error": "No prediction returned"}
    except HTTPException:
        raise
    except ImportError as e:
        raise HTTPException(status_code=503, detail=str(e))
    except Exception as e:  # noqa: BLE001
        raise HTTPException(status_code=500, detail=f"ADMET prediction failed: {e}")


__all__ = ["router"]


