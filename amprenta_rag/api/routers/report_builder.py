"""Report Builder API endpoints."""

from __future__ import annotations

import base64
import logging
from typing import List, Optional
from uuid import UUID

from fastapi import APIRouter, Depends, HTTPException, status
from sqlalchemy.orm import Session

from amprenta_rag.api import schemas
from amprenta_rag.api.dependencies import get_current_user, get_db
from amprenta_rag.database.models import User
from amprenta_rag.services import report_builder as service

logger = logging.getLogger(__name__)

router = APIRouter()


# ============================================================================
# SECTION METADATA
# ============================================================================

@router.get(
    "/sections",
    response_model=List[schemas.SectionMetadataSchema],
    summary="List available section types",
)
def list_available_sections() -> List[schemas.SectionMetadataSchema]:
    """List all available section types with metadata."""
    return [
        schemas.SectionMetadataSchema(
            type=s["type"],
            name=s["name"],
            description=s["description"],
            requires_entity=s["requires_entity"],
            entity_type=s.get("entity_type"),
            icon=s["icon"],
        )
        for s in service.SECTION_REGISTRY
    ]


# ============================================================================
# TEMPLATE CRUD
# ============================================================================

@router.post(
    "/templates",
    response_model=schemas.ReportTemplateResponse,
    status_code=status.HTTP_201_CREATED,
    summary="Create report template",
)
def create_template(
    data: schemas.ReportTemplateCreate,
    db: Session = Depends(get_db),
    current_user: User = Depends(get_current_user),
) -> schemas.ReportTemplateResponse:
    """Create a new report template."""
    try:
        sections = [{"type": s.type, "config": s.config, "order": s.order} for s in data.sections]
        template = service.create_template(
            db=db,
            name=data.name,
            description=data.description,
            sections=sections,
            is_public=data.is_public,
            program_id=data.program_id,
            created_by_id=current_user.id,
        )
        return schemas.ReportTemplateResponse.model_validate(template)
    except ValueError as e:
        raise HTTPException(status_code=status.HTTP_400_BAD_REQUEST, detail=str(e))


@router.get(
    "/templates",
    response_model=List[schemas.ReportTemplateResponse],
    summary="List report templates",
)
def list_templates(
    program_id: Optional[UUID] = None,
    include_public: bool = True,
    limit: int = 100,
    db: Session = Depends(get_db),
    current_user: User = Depends(get_current_user),
) -> List[schemas.ReportTemplateResponse]:
    """List user's templates and optionally public templates."""
    templates = service.list_templates(
        db=db,
        created_by_id=current_user.id,
        program_id=program_id,
        include_public=include_public,
        limit=limit,
    )
    return [schemas.ReportTemplateResponse.model_validate(t) for t in templates]


@router.get(
    "/templates/{template_id}",
    response_model=schemas.ReportTemplateResponse,
    summary="Get template details",
)
def get_template(
    template_id: UUID,
    db: Session = Depends(get_db),
    current_user: User = Depends(get_current_user),
) -> schemas.ReportTemplateResponse:
    """Get a template by ID."""
    template = service.get_template(db, template_id)
    if not template:
        raise HTTPException(status_code=status.HTTP_404_NOT_FOUND, detail="Template not found")
    return schemas.ReportTemplateResponse.model_validate(template)


@router.put(
    "/templates/{template_id}",
    response_model=schemas.ReportTemplateResponse,
    summary="Update template",
)
def update_template(
    template_id: UUID,
    data: schemas.ReportTemplateUpdate,
    db: Session = Depends(get_db),
    current_user: User = Depends(get_current_user),
) -> schemas.ReportTemplateResponse:
    """Update a template."""
    try:
        sections = None
        if data.sections is not None:
            sections = [{"type": s.type, "config": s.config, "order": s.order} for s in data.sections]
        
        template = service.update_template(
            db=db,
            template_id=template_id,
            name=data.name,
            description=data.description,
            sections=sections,
            is_public=data.is_public,
        )
        return schemas.ReportTemplateResponse.model_validate(template)
    except ValueError as e:
        raise HTTPException(status_code=status.HTTP_400_BAD_REQUEST, detail=str(e))


@router.delete(
    "/templates/{template_id}",
    status_code=status.HTTP_204_NO_CONTENT,
    summary="Delete template",
)
def delete_template(
    template_id: UUID,
    db: Session = Depends(get_db),
    current_user: User = Depends(get_current_user),
):
    """Delete a template."""
    success = service.delete_template(db, template_id)
    if not success:
        raise HTTPException(status_code=status.HTTP_404_NOT_FOUND, detail="Template not found")


@router.post(
    "/templates/{template_id}/clone",
    response_model=schemas.ReportTemplateResponse,
    status_code=status.HTTP_201_CREATED,
    summary="Clone template",
)
def clone_template(
    template_id: UUID,
    data: schemas.TemplateCloneRequest,
    db: Session = Depends(get_db),
    current_user: User = Depends(get_current_user),
) -> schemas.ReportTemplateResponse:
    """P1 FIX: Clone an existing template."""
    try:
        template = service.clone_template(
            db=db,
            template_id=template_id,
            new_name=data.new_name,
            created_by_id=current_user.id,
        )
        return schemas.ReportTemplateResponse.model_validate(template)
    except ValueError as e:
        raise HTTPException(status_code=status.HTTP_400_BAD_REQUEST, detail=str(e))


# ============================================================================
# REPORT GENERATION
# ============================================================================

@router.post(
    "/generate",
    response_model=schemas.ReportGenerateResponse,
    summary="Generate report",
)
def generate_report(
    data: schemas.ReportGenerateRequest,
    db: Session = Depends(get_db),
    current_user: User = Depends(get_current_user),
) -> schemas.ReportGenerateResponse:
    """Generate a report from template or ad-hoc sections."""
    # Get sections from template or request
    if data.template_id:
        template = service.get_template(db, data.template_id)
        if not template:
            raise HTTPException(status_code=status.HTTP_404_NOT_FOUND, detail="Template not found")
        sections = template.sections
    elif data.sections:
        sections = [{"type": s.type, "config": s.config, "order": s.order} for s in data.sections]
    else:
        raise HTTPException(
            status_code=status.HTTP_400_BAD_REQUEST,
            detail="Either template_id or sections must be provided"
        )
    
    # Build HTML
    html = service.build_report(sections, db, title=data.title or "Report")
    
    # Return based on format
    if data.format == "pdf":
        try:
            pdf_bytes = service.export_to_pdf(html)
            pdf_b64 = base64.b64encode(pdf_bytes).decode("utf-8")
            return schemas.ReportGenerateResponse(
                format="pdf",
                content_base64=pdf_b64,
            )
        except Exception as e:
            logger.error(f"PDF export failed: {e}")
            raise HTTPException(
                status_code=status.HTTP_500_INTERNAL_SERVER_ERROR,
                detail="PDF export failed"
            )
    else:
        return schemas.ReportGenerateResponse(
            format="html",
            content=html,
        )


@router.post(
    "/preview",
    response_model=schemas.SectionPreviewResponse,
    summary="Preview section",
)
def preview_section(
    data: schemas.SectionPreviewRequest,
    db: Session = Depends(get_db),
    current_user: User = Depends(get_current_user),
) -> schemas.SectionPreviewResponse:
    """Preview a single section (for live preview in builder)."""
    html = service.render_section(data.type, data.config, db)
    return schemas.SectionPreviewResponse(html=html)
