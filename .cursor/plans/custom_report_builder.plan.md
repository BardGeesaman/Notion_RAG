# Custom Report Builder

## Overview

Build a drag-and-drop report builder that allows users to create custom reports by selecting sections, entities, and visualizations. Users can save templates for reuse and export to PDF/HTML.

## Architecture

```
┌─────────────────────────────────────────────────────────────────┐
│                    Custom Report Builder                        │
├─────────────────────────────────────────────────────────────────┤
│  Dashboard UI (Streamlit)                                       │
│  ├── Section Library (draggable components)                     │
│  ├── Report Canvas (drop zone, reorder)                        │
│  ├── Section Config (entity picker, options)                   │
│  ├── Preview Pane (live preview)                               │
│  └── Export Controls (PDF/HTML/save template)                  │
├─────────────────────────────────────────────────────────────────┤
│  API Layer                                                      │
│  ├── POST /report-builder/templates (save)                     │
│  ├── GET /report-builder/templates (list)                      │
│  ├── GET /report-builder/templates/{id} (load)                 │
│  ├── DELETE /report-builder/templates/{id}                     │
│  ├── POST /report-builder/generate (render report)             │
│  └── GET /report-builder/sections (available sections)         │
├─────────────────────────────────────────────────────────────────┤
│  Service Layer                                                  │
│  ├── build_report() - Assemble sections into HTML              │
│  ├── render_section() - Render individual section              │
│  ├── export_pdf() - Convert HTML to PDF                        │
│  └── Template CRUD                                              │
├─────────────────────────────────────────────────────────────────┤
│  Database                                                       │
│  └── ReportTemplate (sections JSON, created_by, program_id)    │
└─────────────────────────────────────────────────────────────────┘
```

## Available Section Types

| Section Type | Description | Entity Required | Config Options |
|-------------|-------------|-----------------|----------------|
| **title_page** | Report title, date, author | No | title, subtitle, logo |
| **executive_summary** | Auto-generated summary | Optional | max_length |
| **compound_profile** | Compound properties, structure | compound_id | show_structure, show_admet |
| **compound_table** | Table of multiple compounds | compound_ids[] | columns[] |
| **experiment_summary** | Experiment overview | experiment_id | show_datasets |
| **dataset_stats** | Dataset statistics | dataset_id | show_chart |
| **activity_chart** | Activity data visualization | compound_id or experiment_id | chart_type |
| **admet_radar** | ADMET radar chart | compound_id | endpoints[] |
| **signature_heatmap** | Signature feature heatmap | signature_id | top_n |
| **pathway_enrichment** | Pathway analysis results | signature_id | top_pathways |
| **free_text** | Custom markdown text | No | content |
| **image** | Custom image/figure | No | image_url, caption |
| **table_of_contents** | Auto-generated TOC | No | - |
| **appendix** | References, methods | No | content |

## Implementation Batches

### Batch 1: Database Model + Service Foundation (Day 1)

**1.1 Create `amprenta_rag/models/reports.py`:**

```python
class ReportTemplate(Base):
    """User-created report template."""
    __tablename__ = "report_templates"
    
    id = Column(UUID, primary_key=True, default=uuid4)
    name = Column(String(255), nullable=False)
    description = Column(Text, nullable=True)
    sections = Column(JSON, nullable=False)  # [{type, config, order}]
    is_public = Column(Boolean, default=False)
    program_id = Column(UUID, ForeignKey("programs.id"), nullable=True)
    created_by_id = Column(UUID, ForeignKey("users.id"), nullable=True)
    created_at = Column(DateTime, default=datetime.utcnow)
    updated_at = Column(DateTime, onupdate=datetime.utcnow)
```

**1.2 Create `amprenta_rag/services/report_builder.py`:**

Core service with:
- `SECTION_REGISTRY` - Available section types with metadata
- `render_section(section_type, config, db)` - Render single section to HTML
- `build_report(sections, db)` - Assemble full report HTML
- `export_to_pdf(html)` - Convert to PDF bytes
- Template CRUD functions

**1.3 Alembic migration**

### Batch 2: Section Renderers (Day 2-3)

Implement renderers for each section type:

```python
# Section renderer pattern
def render_compound_profile(config: dict, db: Session) -> str:
    """Render compound profile section."""
    compound_id = config.get("compound_id")
    compound = db.query(Compound).filter(Compound.id == compound_id).first()
    if not compound:
        return "<p>Compound not found</p>"
    
    html = f"""
    <div class="section compound-profile">
        <h2>{compound.compound_id}</h2>
        <div class="structure">
            <img src="data:image/svg+xml;base64,{get_structure_svg_base64(compound.smiles)}" />
        </div>
        <table class="properties">
            <tr><td>SMILES</td><td>{compound.smiles}</td></tr>
            <tr><td>MW</td><td>{compound.molecular_weight}</td></tr>
            <tr><td>LogP</td><td>{compound.logp}</td></tr>
        </table>
    </div>
    """
    return html
```

Sections to implement:
1. `title_page` - Static HTML with placeholders
2. `executive_summary` - LLM-generated or template
3. `compound_profile` - RDKit structure + properties
4. `compound_table` - Pandas DataFrame to HTML
5. `experiment_summary` - Experiment metadata
6. `dataset_stats` - Feature counts, sample counts
7. `activity_chart` - Plotly to static image
8. `admet_radar` - Plotly radar chart
9. `signature_heatmap` - Seaborn/Plotly heatmap
10. `pathway_enrichment` - Bar chart + table
11. `free_text` - Markdown to HTML
12. `image` - Base64 embed
13. `table_of_contents` - Auto-extract h2/h3 tags
14. `appendix` - Markdown to HTML

### Batch 3: API Endpoints (Day 4)

**Create `amprenta_rag/api/routers/report_builder.py`:**

```python
@router.get("/sections")
def list_available_sections() -> List[SectionMetadata]:
    """List available section types with metadata."""

@router.post("/templates", status_code=201)
def create_template(data: ReportTemplateCreate) -> ReportTemplateResponse:
    """Save a report template."""

@router.get("/templates")
def list_templates(program_id: Optional[UUID] = None) -> List[ReportTemplateResponse]:
    """List user's templates."""

@router.get("/templates/{template_id}")
def get_template(template_id: UUID) -> ReportTemplateResponse:
    """Get template details."""

@router.put("/templates/{template_id}")
def update_template(template_id: UUID, data: ReportTemplateUpdate) -> ReportTemplateResponse:
    """Update template."""

@router.delete("/templates/{template_id}", status_code=204)
def delete_template(template_id: UUID):
    """Delete template."""

@router.post("/generate")
def generate_report(data: ReportGenerateRequest) -> ReportGenerateResponse:
    """Generate report from sections or template."""

@router.post("/preview")
def preview_section(data: SectionPreviewRequest) -> SectionPreviewResponse:
    """Preview a single section (for live preview)."""
```

**Schemas:**

```python
class SectionConfig(BaseSchema):
    type: str
    config: dict
    order: int

class ReportTemplateCreate(BaseSchema):
    name: str
    description: Optional[str]
    sections: List[SectionConfig]
    is_public: bool = False
    program_id: Optional[UUID]

class ReportGenerateRequest(BaseSchema):
    template_id: Optional[UUID]  # Use saved template
    sections: Optional[List[SectionConfig]]  # Or ad-hoc sections
    format: Literal["html", "pdf"] = "html"
    title: Optional[str]

class ReportGenerateResponse(BaseSchema):
    content: Optional[str]  # HTML content or base64 PDF
    download_url: Optional[str]
```

### Batch 4: Dashboard UI (Day 5-6)

**Create `scripts/dashboard/pages/report_builder.py`:**

5-tab interface:
1. **Build Report** - Main builder with drag-drop
2. **My Templates** - Saved templates list
3. **Preview** - Full report preview
4. **Export** - Download PDF/HTML
5. **Section Library** - Browse available sections

**UI Components:**

```python
def render_report_builder_page():
    st.header("📄 Custom Report Builder")
    
    tabs = st.tabs(["Build Report", "My Templates", "Preview", "Export", "Section Library"])
    
    with tabs[0]:  # Build Report
        col1, col2 = st.columns([1, 2])
        
        with col1:
            # Section Library (draggable)
            st.subheader("Available Sections")
            for section in SECTION_REGISTRY:
                if st.button(f"➕ {section['name']}", key=f"add_{section['type']}"):
                    add_section_to_report(section['type'])
        
        with col2:
            # Report Canvas
            st.subheader("Report Sections")
            for i, section in enumerate(st.session_state.report_sections):
                with st.expander(f"{i+1}. {section['type']}", expanded=True):
                    render_section_config(section, i)
                    col1, col2 = st.columns(2)
                    with col1:
                        if st.button("⬆️", key=f"up_{i}") and i > 0:
                            move_section(i, i-1)
                    with col2:
                        if st.button("🗑️", key=f"del_{i}"):
                            remove_section(i)
    
    with tabs[1]:  # My Templates
        render_templates_tab()
    
    with tabs[2]:  # Preview
        render_preview_tab()
    
    with tabs[3]:  # Export
        render_export_tab()
    
    with tabs[4]:  # Section Library
        render_section_library_tab()
```

**Entity Pickers:**

For sections requiring entities:
- Compound picker (searchable dropdown)
- Experiment picker
- Dataset picker
- Signature picker

**Live Preview:**

```python
def render_preview_tab():
    if st.button("🔄 Refresh Preview"):
        sections = st.session_state.report_sections
        html = build_report_preview(sections)
        st.session_state.preview_html = html
    
    if st.session_state.get("preview_html"):
        st.components.v1.html(st.session_state.preview_html, height=800, scrolling=True)
```

### Batch 5: Tests (Day 7)

**Unit Tests (`amprenta_rag/tests/services/test_report_builder_service.py`):**

1. `test_render_section_title_page` - Title page renders correctly
2. `test_render_section_compound_profile` - Compound profile with structure
3. `test_render_section_compound_not_found` - Graceful handling
4. `test_render_section_free_text` - Markdown conversion
5. `test_build_report_multiple_sections` - Full report assembly
6. `test_build_report_empty_sections` - Edge case
7. `test_export_to_pdf` - PDF generation
8. `test_template_crud` - Create, read, update, delete

**API Tests (`amprenta_rag/tests/api/test_report_builder_api.py`):**

1. `test_list_available_sections` - Returns section registry
2. `test_create_template` - Template creation
3. `test_create_template_validation` - Invalid sections rejected
4. `test_list_templates` - User's templates
5. `test_get_template` - Template details
6. `test_update_template` - Template update
7. `test_delete_template` - Template deletion
8. `test_generate_report_from_template` - Generation from saved template
9. `test_generate_report_adhoc` - Ad-hoc generation
10. `test_preview_section` - Single section preview

**E2E Tests (`amprenta_rag/tests/e2e/test_report_builder_e2e.py`):**

1. `test_page_loads` - Builder page accessible
2. `test_add_section` - Add section to canvas
3. `test_reorder_sections` - Move sections up/down
4. `test_remove_section` - Delete section
5. `test_save_template` - Save as template
6. `test_load_template` - Load saved template
7. `test_export_pdf` - PDF download

## Verification Commands

```bash
cd /Users/bard/Documents/RAG && conda activate myenv

# Model imports
python -c "from amprenta_rag.models.reports import ReportTemplate; print('Model OK')"

# Service imports
python -c "from amprenta_rag.services.report_builder import build_report, SECTION_REGISTRY; print('Service OK')"

# API imports
python -c "from amprenta_rag.api.routers.report_builder import router; print('API OK')"

# Run tests
pytest amprenta_rag/tests/services/test_report_builder_service.py -v
pytest amprenta_rag/tests/api/test_report_builder_api.py -v
```

## Dependencies

- `weasyprint` - Already installed (PDF export)
- `markdown` - Already installed (Markdown to HTML)
- `plotly` - Already installed (charts)
- `rdkit` - Already installed (molecule images)
- `pandas` - Already installed (tables)

## UI/UX Considerations

1. **Drag-and-Drop**: Use Streamlit session state for section ordering (no actual drag-drop, but reorder buttons)
2. **Entity Search**: Autocomplete dropdowns for compound/experiment/dataset pickers
3. **Live Preview**: Debounced preview updates
4. **Template Sharing**: Public templates visible to team
5. **Program Scoping**: Templates can be program-specific

## Success Criteria

- [ ] Users can build reports by adding/removing/reordering sections
- [ ] At least 10 section types available
- [ ] Reports export to PDF and HTML
- [ ] Templates can be saved and reused
- [ ] Live preview shows report structure
- [ ] 25+ tests passing

## Estimated Effort

| Batch | Effort | Description |
|-------|--------|-------------|
| 1 | 0.5 days | Database model + service foundation |
| 2 | 1.5 days | Section renderers (10+ types) |
| 3 | 0.5 days | API endpoints |
| 4 | 1.5 days | Dashboard UI |
| 5 | 1 day | Tests |
| **Total** | **5 days** | |

---

*Plan created: 2025-01-03*

