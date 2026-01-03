"""Custom Report Builder Dashboard."""

from __future__ import annotations

import streamlit as st
from uuid import UUID

from amprenta_rag.database.session import db_session
from amprenta_rag.database.models import Experiment, Dataset, Signature, Program
from amprenta_rag.models.chemistry import Compound
from amprenta_rag.services import report_builder as service
from amprenta_rag.auth.session import get_current_user


# ============================================================================
# HELPER FUNCTIONS
# ============================================================================

def get_compounds(limit: int = 100):
    """Fetch compounds for dropdown."""
    with db_session() as db:
        compounds = db.query(Compound).order_by(Compound.compound_id).limit(limit).all()
        return [(str(c.id), c.compound_id) for c in compounds]


def get_experiments(limit: int = 100):
    """Fetch experiments for dropdown."""
    with db_session() as db:
        experiments = db.query(Experiment).order_by(Experiment.name).limit(limit).all()
        return [(str(e.id), e.name) for e in experiments]


def get_datasets(limit: int = 100):
    """Fetch datasets for dropdown."""
    with db_session() as db:
        datasets = db.query(Dataset).order_by(Dataset.name).limit(limit).all()
        return [(str(d.id), d.name) for d in datasets]


def get_signatures(limit: int = 100):
    """Fetch signatures for dropdown."""
    with db_session() as db:
        signatures = db.query(Signature).order_by(Signature.name).limit(limit).all()
        return [(str(s.id), s.name) for s in signatures]


def get_programs(limit: int = 50):
    """Fetch programs for dropdown."""
    with db_session() as db:
        programs = db.query(Program).order_by(Program.name).limit(limit).all()
        return [(str(p.id), p.name) for p in programs]


def add_section(section_type: str, config: dict = None):
    """Add a section to the report."""
    section_info = service.SECTION_TYPE_MAP.get(section_type, {})
    st.session_state.report_sections.append({
        "type": section_type,
        "config": config or {},
        "order": len(st.session_state.report_sections),
        "name": section_info.get("name", section_type),
    })


def remove_section(index: int):
    """Remove a section from the report."""
    st.session_state.report_sections.pop(index)
    # Reorder
    for i, s in enumerate(st.session_state.report_sections):
        s["order"] = i


def move_section(from_idx: int, to_idx: int):
    """Move a section up or down."""
    sections = st.session_state.report_sections
    sections[from_idx], sections[to_idx] = sections[to_idx], sections[from_idx]
    for i, s in enumerate(sections):
        s["order"] = i


def save_template(name: str, description: str, is_public: bool):
    """Save current report as template."""
    try:
        user = get_current_user()
        template_data = {
            "name": name,
            "description": description,
            "sections": st.session_state.report_sections,
            "is_public": is_public,
            "program_id": user.get("program_id") if user else None,
            "created_by_id": user.get("id") if user else None
        }
        
        with db_session() as db:
            template = service.create_template(db, template_data)
            st.success(f"Template saved with ID: {template.id}")
    except Exception as e:
        st.error(f"Error saving template: {e}")


def load_template(template_id: str):
    """Load template into current report."""
    try:
        with db_session() as db:
            template = service.get_template(db, UUID(template_id))
            if template:
                st.session_state.report_sections = template.sections
                st.success(f"Loaded template: {template.name}")
            else:
                st.error("Template not found")
    except Exception as e:
        st.error(f"Error loading template: {e}")


def delete_template(template_id: str):
    """Delete a template."""
    try:
        with db_session() as db:
            service.delete_template(db, UUID(template_id))
            st.success("Template deleted")
    except Exception as e:
        st.error(f"Error deleting template: {e}")


def get_templates():
    """Get user's templates."""
    try:
        user = get_current_user()
        with db_session() as db:
            templates = service.list_templates(db, created_by_id=UUID(user["id"]) if user else None)
            return [(t.id, t.name, t.description, t.is_public, t.created_at) for t in templates]
    except Exception as e:
        st.error(f"Error loading templates: {e}")
        return []


def generate_preview():
    """Generate HTML preview of current report."""
    try:
        with db_session() as db:
            html = service.generate_report(db, st.session_state.report_sections)
            st.session_state.preview_html = html
    except Exception as e:
        st.error(f"Error generating preview: {e}")


def export_html():
    """Export report as HTML file."""
    try:
        with db_session() as db:
            html = service.generate_report(db, st.session_state.report_sections)
            st.download_button(
                label="Download HTML",
                data=html,
                file_name="report.html",
                mime="text/html"
            )
    except Exception as e:
        st.error(f"Error exporting HTML: {e}")


def export_pdf():
    """Export report as PDF file."""
    try:
        with db_session() as db:
            pdf_bytes = service.generate_pdf_report(db, st.session_state.report_sections)
            st.download_button(
                label="Download PDF",
                data=pdf_bytes,
                file_name="report.pdf",
                mime="application/pdf"
            )
    except Exception as e:
        st.error(f"Error exporting PDF: {e}")


def preview_section(section_type: str):
    """Preview a single section."""
    try:
        with db_session() as db:
            html = service.render_section(db, section_type, {})
            st.components.v1.html(html, height=300, scrolling=True)
    except Exception as e:
        st.error(f"Error previewing section: {e}")


# ============================================================================
# REGISTER IN PAGE_REGISTRY
# ============================================================================

def main():
    """Entry point for page registry."""
    st.set_page_config(page_title="Report Builder", page_icon="📄", layout="wide")
    st.title("📄 Custom Report Builder")

    # Initialize session state
    if "report_sections" not in st.session_state:
        st.session_state.report_sections = []
    if "preview_html" not in st.session_state:
        st.session_state.preview_html = None

    # ============================================================================
    # TABS
    # ============================================================================

    tab1, tab2, tab3, tab4, tab5 = st.tabs([
        "🔧 Build Report",
        "📁 My Templates", 
        "👀 Preview",
        "📤 Export",
        "📚 Section Library"
    ])

    # ============================================================================
    # TAB 1: BUILD REPORT
    # ============================================================================

    with tab1:
        st.header("Build Report")
        
        # Section configuration
        col1, col2 = st.columns([1, 1])
        
        with col1:
            section_type = st.selectbox(
                "Section Type",
                options=list(service.get_section_registry().keys()),
                help="Choose the type of section to add"
            )
            
        with col2:
            if section_type:
                registry = service.get_section_registry()
                section_info = registry[section_type]
                st.info(f"**{section_info['name']}**: {section_info['description']}")
        
        # Section configuration based on type
        config = {}
        
        if section_type == "compound_table":
            compound_options = get_compounds()
            selected_compounds = st.multiselect(
                "Select Compounds",
                options=[opt[0] for opt in compound_options],
                format_func=lambda x: next(opt[1] for opt in compound_options if opt[0] == x)
            )
            config = {"compound_ids": selected_compounds}
            
        elif section_type == "experiment_summary":
            experiment_options = get_experiments()
            selected_experiment = st.selectbox(
                "Select Experiment",
                options=[opt[0] for opt in experiment_options] if experiment_options else [],
                format_func=lambda x: next(opt[1] for opt in experiment_options if opt[0] == x) if experiment_options else ""
            )
            if selected_experiment:
                config = {"experiment_id": selected_experiment}
                
        elif section_type == "dataset_stats":
            dataset_options = get_datasets()
            selected_dataset = st.selectbox(
                "Select Dataset", 
                options=[opt[0] for opt in dataset_options] if dataset_options else [],
                format_func=lambda x: next(opt[1] for opt in dataset_options if opt[0] == x) if dataset_options else ""
            )
            if selected_dataset:
                config = {"dataset_id": selected_dataset}
                
        elif section_type == "signature_heatmap":
            signature_options = get_signatures()
            selected_signature = st.selectbox(
                "Select Signature",
                options=[opt[0] for opt in signature_options] if signature_options else [],
                format_func=lambda x: next(opt[1] for opt in signature_options if opt[0] == x) if signature_options else ""
            )
            if selected_signature:
                config = {"signature_id": selected_signature}
        
        # Add section button
        if st.button("Add Section", type="primary"):
            add_section(section_type, config)
            st.success(f"Added {section_type} section")
            st.rerun()
        
        # Current sections
        st.subheader("Current Sections")
        if st.session_state.report_sections:
            for i, section in enumerate(st.session_state.report_sections):
                col1, col2, col3 = st.columns([3, 1, 1])
                with col1:
                    st.write(f"{i+1}. {section['name']}")
                with col2:
                    if st.button("↑", key=f"up_{i}", disabled=i==0):
                        move_section(i, i-1)
                        st.rerun()
                    if st.button("↓", key=f"down_{i}", disabled=i==len(st.session_state.report_sections)-1):
                        move_section(i, i+1)
                        st.rerun()
                with col3:
                    if st.button("🗑️", key=f"delete_{i}"):
                        remove_section(i)
                        st.rerun()
        else:
            st.info("No sections added yet. Add sections above to build your report.")

    # ============================================================================
    # TAB 2: MY TEMPLATES
    # ============================================================================

    with tab2:
        st.header("My Templates")
        
        col1, col2 = st.columns([2, 1])
        
        with col1:
            # Save current report as template
            st.subheader("Save Current Report")
            template_name = st.text_input("Template Name")
            template_description = st.text_area("Description")
            is_public = st.checkbox("Make Public", help="Allow other users to use this template")
            
            has_sections = bool(st.session_state.report_sections)
            if st.button("💾 Save Current as Template", disabled=not has_sections):
                if template_name:
                    save_template(template_name, template_description, is_public)
                    st.success(f"Template '{template_name}' saved!")
                    st.rerun()
                else:
                    st.error("Please enter a template name")
            
            if not has_sections:
                st.info("Add sections to your report first, then save as template")
        
        with col2:
            # Load template
            st.subheader("📁 Saved Templates")
            templates = get_templates()
            
            if templates:
                template_options = [(str(t[0]), f"{t[1]} ({'Public' if t[3] else 'Private'})") for t in templates]
                selected_template = st.selectbox(
                    "Choose Template",
                    options=[opt[0] for opt in template_options],
                    format_func=lambda x: next(opt[1] for opt in template_options if opt[0] == x)
                )
                
                if st.button("Load Template"):
                    load_template(selected_template)
                    st.success("Template loaded!")
                    st.rerun()
            else:
                st.info("No templates available")
        
        # Template management
        st.subheader("Manage Templates")
        templates = get_templates()
        if templates:
            for template_id, name, description, is_public, created_at in templates:
                with st.expander(f"{name} ({'Public' if is_public else 'Private'})"):
                    st.write(f"**Description:** {description or 'No description'}")
                    st.write(f"**Created:** {created_at}")
                    
                    col1, col2 = st.columns(2)
                    with col1:
                        if st.button(f"Load", key=f"load_{template_id}"):
                            load_template(str(template_id))
                            st.success("Template loaded!")
                            st.rerun()
                    with col2:
                        if st.button(f"Delete", key=f"delete_{template_id}"):
                            delete_template(str(template_id))
                            st.success("Template deleted!")
                            st.rerun()

    # ============================================================================
    # TAB 3: PREVIEW
    # ============================================================================

    with tab3:
        st.header("👁️ Report Preview")
        
        has_sections = bool(st.session_state.report_sections)
        if st.button("🔄 Generate Preview", disabled=not has_sections):
            generate_preview()
            
        if st.session_state.preview_html:
            st.components.v1.html(st.session_state.preview_html, height=600, scrolling=True)
        elif not has_sections:
            st.info("Add sections to preview the report")

    # ============================================================================
    # TAB 4: EXPORT
    # ============================================================================

    with tab4:
        st.header("📥 Export Report")
        
        # Report Title
        report_title = st.text_input("Report Title", value="My Report")
        
        # Format selection
        format_option = st.selectbox("Format", ["HTML", "PDF", "Word"])
        
        # Generate & Download button
        has_sections = bool(st.session_state.report_sections)
        if st.button("📥 Generate & Download", disabled=not has_sections):
            if format_option == "HTML":
                export_html()
            elif format_option == "PDF":
                export_pdf()
            else:
                st.info("Word export coming soon!")
        
        if not has_sections:
            st.info("Add sections to export the report")

    # ============================================================================
    # TAB 5: SECTION LIBRARY
    # ============================================================================

    with tab5:
        st.header("📚 Available Section Types")
        st.write("Browse all available section types")
        
        registry = service.get_section_registry()
        
        for section_type, section_info in registry.items():
            with st.expander(f"{section_info['name']} ({section_type})"):
                st.write(f"**Description:** {section_info['description']}")
                st.write(f"**Category:** {section_info.get('category', 'General')}")
                
                # Show preview
                if st.button(f"Preview {section_info['name']}", key=f"preview_{section_type}"):
                    preview_section(section_type)
                    
                # Quick add button
                if st.button(f"Add {section_info['name']}", key=f"add_{section_type}"):
                    add_section(section_type)
                    st.success(f"Added {section_info['name']} to report")
                    st.rerun()


# Add to PAGE_REGISTRY in scripts/dashboard/core/__init__.py:
# "Report Builder": {"module": "scripts.dashboard.pages.report_builder", "icon": "📄", "group": "Collaboration"}
