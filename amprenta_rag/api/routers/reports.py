"""
API router for Report Generation.
"""

import os
from pathlib import Path
from typing import Dict, Any

from fastapi import APIRouter, HTTPException, Response
from fastapi.responses import FileResponse
from pydantic import BaseModel

from amprenta_rag.reports import generate_report

router = APIRouter()


class ReportGenerateRequest(BaseModel):
    """Request schema for report generation."""
    template_name: str
    parameters: Dict[str, Any] = {}
    format: str = "html"  # "html" or "pdf"


@router.post("/generate")
async def generate_report_endpoint(request: ReportGenerateRequest):
    """Generate a report from a notebook template.
    
    Args:
        request: Report generation request with template_name, parameters, and format
        
    Returns:
        Generated report file as download
    """
    try:
        # Generate the report
        report_path = generate_report(
            template_name=request.template_name,
            params=request.parameters,
            format=request.format
        )
        
        if not os.path.exists(report_path):
            raise HTTPException(status_code=500, detail="Report generation failed")
        
        # Determine media type
        if request.format == "pdf":
            media_type = "application/pdf"
        else:
            media_type = "text/html"
        
        # Return file as download
        return FileResponse(
            report_path,
            media_type=media_type,
            filename=Path(report_path).name,
            headers={"Content-Disposition": f'attachment; filename="{Path(report_path).name}"'}
        )
    
    except FileNotFoundError as e:
        raise HTTPException(status_code=404, detail=str(e))
    except ValueError as e:
        raise HTTPException(status_code=400, detail=str(e))
    except Exception as e:
        raise HTTPException(status_code=500, detail=f"Report generation error: {str(e)}")

