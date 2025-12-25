"""Biomarker discovery API endpoints."""

from __future__ import annotations

from typing import Any, Dict, List

from fastapi import APIRouter, HTTPException

from amprenta_rag.api.schemas import BiomarkerDiscoverRequest, BiomarkerDiscoverResponse, BiomarkerFeature
from amprenta_rag.ml.biomarker.discovery import BiomarkerDiscoveryService


router = APIRouter(prefix="/biomarker", tags=["biomarker"])

AVAILABLE_METHODS = ["statistical", "stability", "importance"]


@router.get("/methods", response_model=List[str])
def list_methods() -> List[str]:
    return list(AVAILABLE_METHODS)


@router.post("/discover", response_model=BiomarkerDiscoverResponse)
def discover_biomarkers(request: BiomarkerDiscoverRequest) -> BiomarkerDiscoverResponse:
    methods = [str(m).lower() for m in (request.methods or [])]
    if not methods:
        raise HTTPException(status_code=400, detail="At least one method is required")
    bad = [m for m in methods if m not in AVAILABLE_METHODS]
    if bad:
        raise HTTPException(status_code=400, detail=f"Unknown methods: {bad}")

    svc = BiomarkerDiscoveryService()
    try:
        out = svc.discover(
            experiment_id=request.experiment_id,
            group1=request.group1_samples,
            group2=request.group2_samples,
            methods=methods,
        )
    except ImportError as e:
        raise HTTPException(status_code=400, detail=str(e))
    except ValueError as e:
        raise HTTPException(status_code=400, detail=str(e))
    except Exception as e:  # noqa: BLE001
        raise HTTPException(status_code=500, detail=str(e))

    # Optional: apply FDR threshold to statistical results and consensus.
    fdr_thr = float(request.fdr_threshold) if request.fdr_threshold is not None else None
    method_results: Dict[str, Any] = out.get("method_results") if isinstance(out, dict) else {}
    consensus_raw = out.get("consensus_ranking") if isinstance(out, dict) else []

    if fdr_thr is not None and "statistical" in methods and isinstance(method_results, dict):
        stats_rows = method_results.get("statistical")
        if isinstance(stats_rows, list):
            filtered = []
            allowed = set()
            for r in stats_rows:
                if not isinstance(r, dict):
                    continue
                p_adj = r.get("p_adj", r.get("p_value"))
                try:
                    if p_adj is not None and float(p_adj) <= fdr_thr:
                        filtered.append(r)
                        allowed.add(str(r.get("feature") or ""))
                except Exception:
                    continue
            method_results["statistical"] = filtered
            if isinstance(consensus_raw, list) and allowed:
                consensus_raw = [c for c in consensus_raw if isinstance(c, dict) and str(c.get("feature") or "") in allowed]

    consensus: List[BiomarkerFeature] = []
    if isinstance(consensus_raw, list):
        for r in consensus_raw:
            if not isinstance(r, dict):
                continue
            consensus.append(BiomarkerFeature.model_validate(r))

    return BiomarkerDiscoverResponse(consensus_ranking=consensus, method_results=method_results or {})


__all__ = ["router"]


