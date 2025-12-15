from __future__ import annotations

from fastapi import APIRouter, HTTPException

from amprenta_rag.api import schemas
from amprenta_rag.reports.generator import generate_report

router = APIRouter()


@router.post(
    "/reports/generate",
    summary="Generate narrative report",
    response_model=schemas.ReportResponse,
)
def generate_narrative_report(request: schemas.ReportRequest) -> schemas.ReportResponse:
    """Generate a narrative report for an entity using the notebook pipeline."""
    try:
        report_path = generate_report(
            template_name="narrative_report.ipynb",
            params={
                "entity_type": request.entity_type,
                "entity_id": request.entity_id,
                "format": request.format,
            },
            format=request.format,
        )
    except FileNotFoundError as exc:
        raise HTTPException(status_code=404, detail=str(exc)) from exc
    except ValueError as exc:
        raise HTTPException(status_code=400, detail=str(exc)) from exc
    except Exception as exc:  # pragma: no cover - fallback
        raise HTTPException(
            status_code=500, detail="Failed to generate report"
        ) from exc

    # Placeholder: could map to download URL if served via static hosting
    return schemas.ReportResponse(file_path=report_path, download_url=None)
