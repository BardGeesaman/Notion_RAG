"""Report Builder service layer."""

from __future__ import annotations

import logging
from datetime import datetime, timezone
from typing import Any, Dict, List, Optional
from uuid import UUID

from sqlalchemy.orm import Session

from amprenta_rag.database.models import ReportTemplate, User

logger = logging.getLogger(__name__)


# ============================================================================
# SECTION REGISTRY - Available section types
# ============================================================================

SECTION_REGISTRY: List[Dict[str, Any]] = [
    {
        "type": "title_page",
        "name": "Title Page",
        "description": "Report title, date, and author information",
        "requires_entity": False,
        "icon": "📄",
    },
    {
        "type": "executive_summary",
        "name": "Executive Summary",
        "description": "Auto-generated or custom summary",
        "requires_entity": False,
        "icon": "📋",
    },
    {
        "type": "compound_profile",
        "name": "Compound Profile",
        "description": "Compound properties and structure image",
        "requires_entity": True,
        "entity_type": "compound",
        "icon": "🧪",
    },
    {
        "type": "compound_table",
        "name": "Compound Table",
        "description": "Table of multiple compounds with properties",
        "requires_entity": True,
        "entity_type": "compounds",
        "icon": "📊",
    },
    {
        "type": "experiment_summary",
        "name": "Experiment Summary",
        "description": "Experiment overview and metadata",
        "requires_entity": True,
        "entity_type": "experiment",
        "icon": "🔬",
    },
    {
        "type": "dataset_stats",
        "name": "Dataset Statistics",
        "description": "Dataset feature counts, sample statistics",
        "requires_entity": True,
        "entity_type": "dataset",
        "icon": "📈",
    },
    {
        "type": "activity_chart",
        "name": "Activity Chart",
        "description": "Activity data visualization",
        "requires_entity": True,
        "entity_type": "compound",
        "icon": "📉",
    },
    {
        "type": "admet_radar",
        "name": "ADMET Radar",
        "description": "ADMET properties radar chart",
        "requires_entity": True,
        "entity_type": "compound",
        "icon": "🎯",
    },
    {
        "type": "signature_heatmap",
        "name": "Signature Heatmap",
        "description": "Top features heatmap visualization",
        "requires_entity": True,
        "entity_type": "signature",
        "icon": "🔥",
    },
    {
        "type": "pathway_enrichment",
        "name": "Pathway Enrichment",
        "description": "Pathway analysis results",
        "requires_entity": True,
        "entity_type": "signature",
        "icon": "🛤️",
    },
    {
        "type": "free_text",
        "name": "Free Text",
        "description": "Custom markdown content",
        "requires_entity": False,
        "icon": "✏️",
    },
    {
        "type": "image",
        "name": "Image",
        "description": "Custom image with caption",
        "requires_entity": False,
        "icon": "🖼️",
    },
    {
        "type": "table_of_contents",
        "name": "Table of Contents",
        "description": "Auto-generated table of contents",
        "requires_entity": False,
        "icon": "📑",
    },
    {
        "type": "appendix",
        "name": "Appendix",
        "description": "References and additional notes",
        "requires_entity": False,
        "icon": "📎",
    },
]

# Map type -> registry entry for quick lookup
SECTION_TYPE_MAP = {s["type"]: s for s in SECTION_REGISTRY}


# ============================================================================
# SECTION SCHEMAS - P1 FIX: Validation schemas for section configs
# ============================================================================

SECTION_SCHEMAS: Dict[str, Dict[str, Any]] = {
    "title_page": {
        "type": "object",
        "properties": {
            "title": {"type": "string"},
            "subtitle": {"type": "string"},
            "show_date": {"type": "boolean"},
            "show_author": {"type": "boolean"},
        },
    },
    "executive_summary": {
        "type": "object",
        "properties": {
            "content": {"type": "string"},
            "max_length": {"type": "integer"},
        },
    },
    "compound_profile": {
        "type": "object",
        "required": ["compound_id"],
        "properties": {
            "compound_id": {"type": "string", "format": "uuid"},
            "show_structure": {"type": "boolean"},
            "show_admet": {"type": "boolean"},
            "show_properties": {"type": "boolean"},
        },
    },
    "compound_table": {
        "type": "object",
        "required": ["compound_ids"],
        "properties": {
            "compound_ids": {"type": "array", "items": {"type": "string"}},
            "columns": {"type": "array", "items": {"type": "string"}},
        },
    },
    "experiment_summary": {
        "type": "object",
        "required": ["experiment_id"],
        "properties": {
            "experiment_id": {"type": "string", "format": "uuid"},
            "show_datasets": {"type": "boolean"},
        },
    },
    "dataset_stats": {
        "type": "object",
        "required": ["dataset_id"],
        "properties": {
            "dataset_id": {"type": "string", "format": "uuid"},
            "show_chart": {"type": "boolean"},
        },
    },
    "activity_chart": {
        "type": "object",
        "required": ["compound_id"],
        "properties": {
            "compound_id": {"type": "string", "format": "uuid"},
            "experiment_id": {"type": "string", "format": "uuid"},
            "chart_type": {"type": "string", "enum": ["bar", "line", "scatter"]},
        },
    },
    "admet_radar": {
        "type": "object",
        "required": ["compound_id"],
        "properties": {
            "compound_id": {"type": "string", "format": "uuid"},
            "endpoints": {"type": "array", "items": {"type": "string"}},
        },
    },
    "signature_heatmap": {
        "type": "object",
        "required": ["signature_id"],
        "properties": {
            "signature_id": {"type": "string", "format": "uuid"},
            "top_n": {"type": "integer"},
        },
    },
    "pathway_enrichment": {
        "type": "object",
        "required": ["signature_id"],
        "properties": {
            "signature_id": {"type": "string", "format": "uuid"},
            "top_pathways": {"type": "integer"},
        },
    },
    "free_text": {
        "type": "object",
        "properties": {
            "content": {"type": "string"},
        },
    },
    "image": {
        "type": "object",
        "properties": {
            "image_url": {"type": "string"},
            "image_base64": {"type": "string"},
            "caption": {"type": "string"},
        },
    },
    "table_of_contents": {
        "type": "object",
        "properties": {},
    },
    "appendix": {
        "type": "object",
        "properties": {
            "content": {"type": "string"},
        },
    },
}


# ============================================================================
# VALIDATION - P1 FIX
# ============================================================================

def validate_section_config(section_type: str, config: dict) -> List[str]:
    """
    Validate section config against schema.
    
    Returns list of validation errors (empty if valid).
    """
    errors = []
    
    if section_type not in SECTION_SCHEMAS:
        errors.append(f"Unknown section type: {section_type}")
        return errors
    
    schema = SECTION_SCHEMAS[section_type]
    
    # Check required fields
    required = schema.get("required", [])
    for field in required:
        if field not in config:
            errors.append(f"Section '{section_type}' missing required field: {field}")
    
    return errors


def validate_sections(sections: List[dict]) -> List[str]:
    """Validate all sections in a template."""
    errors = []
    for i, section in enumerate(sections):
        section_type = section.get("type")
        if not section_type:
            errors.append(f"Section {i}: missing 'type' field")
            continue
        
        config = section.get("config", {})
        section_errors = validate_section_config(section_type, config)
        for err in section_errors:
            errors.append(f"Section {i} ({section_type}): {err}")
    
    return errors


# ============================================================================
# TEMPLATE CRUD
# ============================================================================

def create_template(
    db: Session,
    name: str,
    sections: List[dict],
    description: Optional[str] = None,
    is_public: bool = False,
    program_id: Optional[UUID] = None,
    created_by_id: Optional[UUID] = None,
) -> ReportTemplate:
    """Create a new report template."""
    # Validate sections before save
    validation_errors = validate_sections(sections)
    if validation_errors:
        raise ValueError(f"Invalid sections: {'; '.join(validation_errors)}")
    
    template = ReportTemplate(
        name=name,
        description=description,
        sections=sections,
        is_public=is_public,
        program_id=program_id,
        created_by_id=created_by_id,
    )
    db.add(template)
    db.commit()
    db.refresh(template)
    logger.info(f"Created report template: {template.id} - {name}")
    return template


def get_template(db: Session, template_id: UUID) -> Optional[ReportTemplate]:
    """Get template by ID."""
    return db.query(ReportTemplate).filter(ReportTemplate.id == template_id).first()


def list_templates(
    db: Session,
    created_by_id: Optional[UUID] = None,
    program_id: Optional[UUID] = None,
    include_public: bool = True,
    limit: int = 100,
) -> List[ReportTemplate]:
    """List templates with optional filters."""
    query = db.query(ReportTemplate)
    
    if created_by_id and include_public:
        # User's templates + public templates
        query = query.filter(
            (ReportTemplate.created_by_id == created_by_id) | 
            (ReportTemplate.is_public.is_(True))
        )
    elif created_by_id:
        query = query.filter(ReportTemplate.created_by_id == created_by_id)
    
    if program_id:
        query = query.filter(
            (ReportTemplate.program_id == program_id) | 
            (ReportTemplate.program_id.is_(None))
        )
    
    return query.order_by(ReportTemplate.created_at.desc()).limit(limit).all()


def update_template(
    db: Session,
    template_id: UUID,
    name: Optional[str] = None,
    description: Optional[str] = None,
    sections: Optional[List[dict]] = None,
    is_public: Optional[bool] = None,
) -> ReportTemplate:
    """Update template fields."""
    template = get_template(db, template_id)
    if not template:
        raise ValueError(f"Template {template_id} not found")
    
    if sections is not None:
        validation_errors = validate_sections(sections)
        if validation_errors:
            raise ValueError(f"Invalid sections: {'; '.join(validation_errors)}")
        template.sections = sections
    
    if name is not None:
        template.name = name
    if description is not None:
        template.description = description
    if is_public is not None:
        template.is_public = is_public
    
    db.commit()
    db.refresh(template)
    return template


def delete_template(db: Session, template_id: UUID) -> bool:
    """Delete template."""
    template = get_template(db, template_id)
    if not template:
        return False
    
    db.delete(template)
    db.commit()
    logger.info(f"Deleted report template: {template_id}")
    return True


def clone_template(
    db: Session,
    template_id: UUID,
    new_name: str,
    created_by_id: Optional[UUID] = None,
) -> ReportTemplate:
    """
    P1 FIX: Clone an existing template.
    """
    original = get_template(db, template_id)
    if not original:
        raise ValueError(f"Template {template_id} not found")
    
    return create_template(
        db=db,
        name=new_name,
        sections=original.sections.copy() if original.sections else [],
        description=f"Cloned from: {original.name}",
        is_public=False,  # Clones are private by default
        program_id=original.program_id,
        created_by_id=created_by_id,
    )


# ============================================================================
# SECTION RENDERING - P1 FIX: Error handling wrapper
# ============================================================================

def render_section(section_type: str, config: dict, db: Session) -> str:
    """
    Render a single section to HTML.
    
    P1 FIX: Wrapped with error handling for graceful degradation.
    """
    if section_type not in RENDERERS:
        return f'<div class="section-error">Unknown section type: {section_type}</div>'
    
    try:
        return RENDERERS[section_type](config, db)
    except Exception as e:
        logger.error(f"Section render failed [{section_type}]: {e}")
        return f'<div class="section-error">Section unavailable: {section_type}</div>'


# Placeholder renderers - full implementation in Batch 2
def _render_title_page(config: dict, db: Session) -> str:
    title = config.get("title", "Report")
    subtitle = config.get("subtitle", "")
    return f"""
    <div class="section title-page">
        <h1>{title}</h1>
        {f'<p class="subtitle">{subtitle}</p>' if subtitle else ''}
        <p class="date">{datetime.now(timezone.utc).strftime('%Y-%m-%d')}</p>
    </div>
    """


def _render_free_text(config: dict, db: Session) -> str:
    import markdown
    content = config.get("content", "")
    html = markdown.markdown(content, extensions=["tables", "fenced_code"])
    return f'<div class="section free-text">{html}</div>'


def _render_placeholder(section_type: str):
    """Generate placeholder renderer for sections not yet implemented."""
    def _renderer(config: dict, db: Session) -> str:
        return f'<div class="section placeholder">[{section_type} - Coming in Batch 2]</div>'
    return _renderer


# Renderer registry
RENDERERS: Dict[str, callable] = {
    "title_page": _render_title_page,
    "free_text": _render_free_text,
    # Placeholders for Batch 2
    "executive_summary": _render_placeholder("executive_summary"),
    "compound_profile": _render_placeholder("compound_profile"),
    "compound_table": _render_placeholder("compound_table"),
    "experiment_summary": _render_placeholder("experiment_summary"),
    "dataset_stats": _render_placeholder("dataset_stats"),
    "activity_chart": _render_placeholder("activity_chart"),
    "admet_radar": _render_placeholder("admet_radar"),
    "signature_heatmap": _render_placeholder("signature_heatmap"),
    "pathway_enrichment": _render_placeholder("pathway_enrichment"),
    "image": _render_placeholder("image"),
    "table_of_contents": _render_placeholder("table_of_contents"),
    "appendix": _render_placeholder("appendix"),
}


# ============================================================================
# REPORT BUILDING
# ============================================================================

REPORT_CSS = """
<style>
body { font-family: Georgia, serif; max-width: 800px; margin: 0 auto; padding: 20px; }
h1, h2, h3 { font-family: Arial, sans-serif; color: #1f3b57; }
.title-page { text-align: center; padding: 60px 0; border-bottom: 2px solid #1f3b57; }
.title-page h1 { font-size: 32px; margin-bottom: 10px; }
.title-page .subtitle { font-size: 18px; color: #666; }
.title-page .date { font-size: 14px; color: #999; margin-top: 20px; }
.section { margin: 30px 0; }
.section-error { background: #fee; border: 1px solid #f00; padding: 10px; color: #900; }
.placeholder { background: #f5f5f5; padding: 20px; text-align: center; color: #666; font-style: italic; }
table { width: 100%; border-collapse: collapse; margin: 15px 0; }
th, td { border: 1px solid #ddd; padding: 8px; text-align: left; }
th { background: #f5f5f5; }
</style>
"""


def build_report(sections: List[dict], db: Session, title: str = "Report") -> str:
    """
    Build complete HTML report from sections.
    """
    html_parts = [
        "<!DOCTYPE html>",
        "<html>",
        "<head>",
        f"<title>{title}</title>",
        REPORT_CSS,
        "</head>",
        "<body>",
    ]
    
    # Sort sections by order
    sorted_sections = sorted(sections, key=lambda s: s.get("order", 0))
    
    for section in sorted_sections:
        section_type = section.get("type")
        config = section.get("config", {})
        section_html = render_section(section_type, config, db)
        html_parts.append(section_html)
    
    html_parts.extend(["</body>", "</html>"])
    
    return "\n".join(html_parts)


def export_to_pdf(html: str) -> bytes:
    """Convert HTML to PDF bytes."""
    from weasyprint import HTML
    return HTML(string=html).write_pdf()
