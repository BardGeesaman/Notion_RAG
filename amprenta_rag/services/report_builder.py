"""Report Builder service layer."""

from __future__ import annotations

import base64
import io
import logging
from datetime import datetime, timezone
from typing import Any, Callable, Dict, List, Optional
from uuid import UUID

import markdown
import pandas as pd
from sqlalchemy.orm import Session

from amprenta_rag.database.models import ReportTemplate, User

# Conditional imports for optional features
try:
    from rdkit import Chem
    from rdkit.Chem import Draw
    RDKIT_AVAILABLE = True
except ImportError:
    RDKIT_AVAILABLE = False

try:
    import plotly.graph_objects as go
    import plotly.express as px
    PLOTLY_AVAILABLE = True
except ImportError:
    PLOTLY_AVAILABLE = False

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


def get_section_registry() -> Dict[str, Dict[str, Any]]:
    """Get section registry as a dict keyed by section type."""
    return SECTION_TYPE_MAP


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


# ============================================================================
# SECTION RENDERERS - Full implementations
# ============================================================================

def _render_title_page(config: dict, db: Session) -> str:
    """Render title page section."""
    title = config.get("title", "Report")
    subtitle = config.get("subtitle", "")
    show_date = config.get("show_date", True)
    show_author = config.get("show_author", False)
    
    html = f"""
    <div class="section title-page">
        <h1>{title}</h1>
        {f'<p class="subtitle">{subtitle}</p>' if subtitle else ''}
    """
    
    if show_date:
        html += f'<p class="date">Generated: {datetime.now(timezone.utc).strftime("%Y-%m-%d %H:%M UTC")}</p>'
    
    html += "</div>"
    return html


def _render_executive_summary(config: dict, db: Session) -> str:
    """Render executive summary section."""
    content = config.get("content", "")
    
    if not content:
        content = "_No executive summary provided._"
    
    html_content = markdown.markdown(content, extensions=["tables"])
    
    return f"""
    <div class="section executive-summary">
        <h2>Executive Summary</h2>
        {html_content}
    </div>
    """


def _render_compound_profile(config: dict, db: Session) -> str:
    """Render compound profile with structure and properties."""
    from amprenta_rag.models.chemistry import Compound
    
    compound_id = config.get("compound_id")
    show_structure = config.get("show_structure", True)
    show_properties = config.get("show_properties", True)
    show_admet = config.get("show_admet", False)
    
    compound = db.query(Compound).filter(Compound.id == compound_id).first()
    if not compound:
        return f'<div class="section-error">Compound not found: {compound_id}</div>'
    
    html = f"""
    <div class="section compound-profile">
        <h2>Compound: {compound.compound_id}</h2>
    """
    
    # Structure image
    if show_structure and RDKIT_AVAILABLE and compound.smiles:
        try:
            mol = Chem.MolFromSmiles(compound.smiles)
            if mol:
                img = Draw.MolToImage(mol, size=(300, 200))
                buffer = io.BytesIO()
                img.save(buffer, format="PNG")
                img_b64 = base64.b64encode(buffer.getvalue()).decode()
                html += f'<div class="structure"><img src="data:image/png;base64,{img_b64}" alt="Structure" /></div>'
        except Exception:
            pass  # Skip structure on error
    
    # Properties table
    if show_properties:
        html += """
        <table class="properties">
            <tr><th>Property</th><th>Value</th></tr>
        """
        props = [
            ("SMILES", compound.smiles),
            ("Molecular Weight", f"{compound.molecular_weight:.2f}" if compound.molecular_weight else "-"),
            ("LogP", f"{compound.logp:.2f}" if compound.logp else "-"),
            ("HBD Count", compound.hbd_count),
            ("HBA Count", compound.hba_count),
            ("Rotatable Bonds", compound.rotatable_bonds),
        ]
        for name, value in props:
            html += f"<tr><td>{name}</td><td>{value if value is not None else '-'}</td></tr>"
        html += "</table>"
    
    html += "</div>"
    return html


def _render_compound_table(config: dict, db: Session) -> str:
    """Render table of multiple compounds."""
    from amprenta_rag.models.chemistry import Compound
    
    compound_ids = config.get("compound_ids", [])
    columns = config.get("columns", ["compound_id", "smiles", "molecular_weight", "logp"])
    
    if not compound_ids:
        return '<div class="section-error">No compounds specified</div>'
    
    compounds = db.query(Compound).filter(Compound.id.in_(compound_ids)).all()
    
    if not compounds:
        return '<div class="section-error">No compounds found</div>'
    
    # Build table
    html = """
    <div class="section compound-table">
        <h2>Compounds</h2>
        <table>
            <tr>
    """
    
    # Headers
    col_labels = {
        "compound_id": "ID",
        "smiles": "SMILES",
        "molecular_weight": "MW",
        "logp": "LogP",
        "hbd_count": "HBD",
        "hba_count": "HBA",
    }
    for col in columns:
        html += f"<th>{col_labels.get(col, col)}</th>"
    html += "</tr>"
    
    # Rows
    for compound in compounds:
        html += "<tr>"
        for col in columns:
            val = getattr(compound, col, "-")
            if isinstance(val, float):
                val = f"{val:.2f}"
            html += f"<td>{val if val is not None else '-'}</td>"
        html += "</tr>"
    
    html += "</table></div>"
    return html


def _render_experiment_summary(config: dict, db: Session) -> str:
    """Render experiment summary."""
    from amprenta_rag.database.models import Experiment, Dataset
    
    experiment_id = config.get("experiment_id")
    show_datasets = config.get("show_datasets", True)
    
    experiment = db.query(Experiment).filter(Experiment.id == experiment_id).first()
    if not experiment:
        return f'<div class="section-error">Experiment not found: {experiment_id}</div>'
    
    html = f"""
    <div class="section experiment-summary">
        <h2>Experiment: {experiment.name}</h2>
        <table>
            <tr><td><strong>Description</strong></td><td>{experiment.description or '-'}</td></tr>
            <tr><td><strong>Design Type</strong></td><td>{experiment.design_type or '-'}</td></tr>
            <tr><td><strong>Status</strong></td><td>{experiment.status or '-'}</td></tr>
            <tr><td><strong>Created</strong></td><td>{experiment.created_at.strftime('%Y-%m-%d') if experiment.created_at else '-'}</td></tr>
        </table>
    """
    
    if show_datasets:
        datasets = db.query(Dataset).filter(Dataset.experiment_id == experiment_id).all()
        if datasets:
            html += f"<h3>Datasets ({len(datasets)})</h3><ul>"
            for ds in datasets[:10]:  # Limit to 10
                html += f"<li>{ds.name} ({ds.omics_type or 'unknown'})</li>"
            if len(datasets) > 10:
                html += f"<li>... and {len(datasets) - 10} more</li>"
            html += "</ul>"
    
    html += "</div>"
    return html


def _render_dataset_stats(config: dict, db: Session) -> str:
    """Render dataset statistics."""
    from amprenta_rag.database.models import Dataset, Feature
    
    dataset_id = config.get("dataset_id")
    show_chart = config.get("show_chart", False)
    
    dataset = db.query(Dataset).filter(Dataset.id == dataset_id).first()
    if not dataset:
        return f'<div class="section-error">Dataset not found: {dataset_id}</div>'
    
    # Count features
    feature_count = db.query(Feature).filter(Feature.dataset_id == dataset_id).count()
    
    html = f"""
    <div class="section dataset-stats">
        <h2>Dataset: {dataset.name}</h2>
        <table>
            <tr><td><strong>Type</strong></td><td>{dataset.omics_type or '-'}</td></tr>
            <tr><td><strong>Features</strong></td><td>{feature_count}</td></tr>
            <tr><td><strong>Samples</strong></td><td>{dataset.sample_count or '-'}</td></tr>
            <tr><td><strong>Source</strong></td><td>{dataset.source or '-'}</td></tr>
        </table>
    </div>
    """
    return html


def _render_activity_chart(config: dict, db: Session) -> str:
    """Render activity data chart."""
    compound_id = config.get("compound_id")
    chart_type = config.get("chart_type", "bar")
    
    # Placeholder - would query HTSResult or BiochemicalResult
    html = f"""
    <div class="section activity-chart">
        <h2>Activity Data</h2>
        <p class="placeholder">[Activity chart for compound {str(compound_id)[:8]}... - requires HTS/biochemical data]</p>
    </div>
    """
    return html


def _render_admet_radar(config: dict, db: Session) -> str:
    """Render ADMET radar chart."""
    from amprenta_rag.models.chemistry import Compound
    
    compound_id = config.get("compound_id")
    endpoints = config.get("endpoints", ["hERG", "LogS", "LogP", "BBB", "Caco2"])
    
    compound = db.query(Compound).filter(Compound.id == compound_id).first()
    if not compound:
        return f'<div class="section-error">Compound not found: {compound_id}</div>'
    
    # ADMET radar requires predictions - show placeholder
    html = f"""
    <div class="section admet-radar">
        <h2>ADMET Profile: {compound.compound_id}</h2>
        <p>Endpoints: {', '.join(endpoints)}</p>
        <p class="placeholder">[ADMET radar chart - integrate with ADMET predictor service]</p>
    </div>
    """
    return html


def _render_signature_heatmap(config: dict, db: Session) -> str:
    """Render signature feature heatmap."""
    from amprenta_rag.database.models import Signature
    
    signature_id = config.get("signature_id")
    top_n = config.get("top_n", 20)
    
    signature = db.query(Signature).filter(Signature.id == signature_id).first()
    if not signature:
        return f'<div class="section-error">Signature not found: {signature_id}</div>'
    
    html = f"""
    <div class="section signature-heatmap">
        <h2>Signature: {signature.name}</h2>
        <p>Top {top_n} features</p>
        <p class="placeholder">[Signature heatmap - requires feature data integration]</p>
    </div>
    """
    return html


def _render_pathway_enrichment(config: dict, db: Session) -> str:
    """Render pathway enrichment results."""
    from amprenta_rag.database.models import Signature
    
    signature_id = config.get("signature_id")
    top_pathways = config.get("top_pathways", 10)
    
    signature = db.query(Signature).filter(Signature.id == signature_id).first()
    if not signature:
        return f'<div class="section-error">Signature not found: {signature_id}</div>'
    
    html = f"""
    <div class="section pathway-enrichment">
        <h2>Pathway Enrichment</h2>
        <p>Signature: {signature.name}</p>
        <p>Top {top_pathways} enriched pathways</p>
        <p class="placeholder">[Pathway enrichment - integrate with pathway analysis service]</p>
    </div>
    """
    return html


def _render_free_text(config: dict, db: Session) -> str:
    """Render free text (markdown) section."""
    content = config.get("content", "")
    html_content = markdown.markdown(content, extensions=["tables", "fenced_code"])
    return f'<div class="section free-text">{html_content}</div>'


def _render_image(config: dict, db: Session) -> str:
    """Render image section."""
    image_url = config.get("image_url", "")
    image_base64 = config.get("image_base64", "")
    caption = config.get("caption", "")
    
    if image_base64:
        src = f"data:image/png;base64,{image_base64}"
    elif image_url:
        src = image_url
    else:
        return '<div class="section-error">No image provided</div>'
    
    html = f"""
    <div class="section image">
        <img src="{src}" alt="{caption}" style="max-width: 100%;" />
        {f'<p class="caption">{caption}</p>' if caption else ''}
    </div>
    """
    return html


def _render_table_of_contents(config: dict, db: Session) -> str:
    """Render table of contents placeholder (populated after full render)."""
    return """
    <div class="section table-of-contents">
        <h2>Table of Contents</h2>
        <p class="placeholder">[Auto-generated from section headers]</p>
    </div>
    """


def _render_appendix(config: dict, db: Session) -> str:
    """Render appendix section."""
    content = config.get("content", "")
    html_content = markdown.markdown(content, extensions=["tables", "fenced_code"])
    return f"""
    <div class="section appendix">
        <h2>Appendix</h2>
        {html_content if content else '<p>No appendix content.</p>'}
    </div>
    """


# ============================================================================
# UPDATE RENDERERS DICT
# ============================================================================

RENDERERS: Dict[str, Callable[[dict, Session], str]] = {
    "title_page": _render_title_page,
    "executive_summary": _render_executive_summary,
    "compound_profile": _render_compound_profile,
    "compound_table": _render_compound_table,
    "experiment_summary": _render_experiment_summary,
    "dataset_stats": _render_dataset_stats,
    "activity_chart": _render_activity_chart,
    "admet_radar": _render_admet_radar,
    "signature_heatmap": _render_signature_heatmap,
    "pathway_enrichment": _render_pathway_enrichment,
    "free_text": _render_free_text,
    "image": _render_image,
    "table_of_contents": _render_table_of_contents,
    "appendix": _render_appendix,
}


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
        if not section_type:
            continue  # Skip sections without type
        config = section.get("config", {})
        section_html = render_section(str(section_type), config, db)
        html_parts.append(section_html)
    
    html_parts.extend(["</body>", "</html>"])
    
    return "\n".join(html_parts)


def export_to_pdf(html: str) -> bytes:
    """Convert HTML to PDF bytes."""
    from weasyprint import HTML
    return HTML(string=html).write_pdf()
