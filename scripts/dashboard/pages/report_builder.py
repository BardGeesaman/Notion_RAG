"""Custom Report Builder Dashboard."""

from __future__ import annotations

import streamlit as st
from uuid import UUID

from amprenta_rag.database.session import db_session
from amprenta_rag.database.models import Experiment, Dataset, Signature, Program
from amprenta_rag.models.chemistry import Compound
from amprenta_rag.services import report_builder as service
from amprenta_rag.auth.session import get_current_user


st.set_page_config(page_title="Report Builder", page_icon="📄", layout="wide")
st.title("📄 Custom Report Builder")

# Initialize session state
if "report_sections" not in st.session_state:
    st.session_state.report_sections = []
if "preview_html" not in st.session_state:
    st.session_state.preview_html = None


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


def add_section(section_type: str):
    """Add a section to the report."""
    section_info = service.SECTION_TYPE_MAP.get(section_type, {})
    st.session_state.report_sections.append({
        "type": section_type,
        "config": {},
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


def render_section_config(section: dict, index: int):
    """Render configuration form for a section."""
    section_type = section["type"]
    config = section.get("config", {})
    
    if section_type == "title_page":
        config["title"] = st.text_input("Report Title", value=config.get("title", ""), key=f"title_{index}")
        config["subtitle"] = st.text_input("Subtitle", value=config.get("subtitle", ""), key=f"subtitle_{index}")
        config["show_date"] = st.checkbox("Show Date", value=config.get("show_date", True), key=f"date_{index}")
    
    elif section_type == "executive_summary":
        config["content"] = st.text_area("Summary Content (Markdown)", value=config.get("content", ""), key=f"summary_{index}", height=150)
    
    elif section_type == "compound_profile":
        compounds = get_compounds()
        options = [("", "-- Select Compound --")] + compounds
        selected = st.selectbox(
            "Compound",
            options=[c[0] for c in options],
            format_func=lambda x: dict(options).get(x, x),
            key=f"compound_{index}"
        )
        if selected:
            config["compound_id"] = selected
        config["show_structure"] = st.checkbox("Show Structure", value=config.get("show_structure", True), key=f"struct_{index}")
        config["show_properties"] = st.checkbox("Show Properties", value=config.get("show_properties", True), key=f"props_{index}")
    
    elif section_type == "compound_table":
        compounds = get_compounds()
        options = [c[0] for c in compounds]
        labels = {c[0]: c[1] for c in compounds}
        selected = st.multiselect(
            "Select Compounds",
            options=options,
            format_func=lambda x: labels.get(x, x),
            key=f"compounds_{index}"
        )
        config["compound_ids"] = selected
    
    elif section_type == "experiment_summary":
        experiments = get_experiments()
        options = [("", "-- Select Experiment --")] + experiments
        selected = st.selectbox(
            "Experiment",
            options=[e[0] for e in options],
            format_func=lambda x: dict(options).get(x, x),
            key=f"experiment_{index}"
        )
        if selected:
            config["experiment_id"] = selected
        config["show_datasets"] = st.checkbox("Show Datasets", value=config.get("show_datasets", True), key=f"datasets_{index}")
    
    elif section_type == "dataset_stats":
        datasets = get_datasets()
        options = [("", "-- Select Dataset --")] + datasets
        selected = st.selectbox(
            "Dataset",
            options=[d[0] for d in options],
            format_func=lambda x: dict(options).get(x, x),
            key=f"dataset_{index}"
        )
        if selected:
            config["dataset_id"] = selected
    
    elif section_type in ("signature_heatmap", "pathway_enrichment"):
        signatures = get_signatures()
        options = [("", "-- Select Signature --")] + signatures
        selected = st.selectbox(
            "Signature",
            options=[s[0] for s in options],
            format_func=lambda x: dict(options).get(x, x),
            key=f"signature_{index}"
        )
        if selected:
            config["signature_id"] = selected
        if section_type == "signature_heatmap":
            config["top_n"] = st.number_input("Top N Features", value=config.get("top_n", 20), min_value=5, max_value=100, key=f"topn_{index}")
    
    elif section_type in ("activity_chart", "admet_radar"):
        compounds = get_compounds()
        options = [("", "-- Select Compound --")] + compounds
        selected = st.selectbox(
            "Compound",
            options=[c[0] for c in options],
            format_func=lambda x: dict(options).get(x, x),
            key=f"viz_compound_{index}"
        )
        if selected:
            config["compound_id"] = selected
    
    elif section_type == "free_text":
        config["content"] = st.text_area("Content (Markdown)", value=config.get("content", ""), key=f"text_{index}", height=200)
    
    elif section_type == "image":
        config["image_url"] = st.text_input("Image URL", value=config.get("image_url", ""), key=f"url_{index}")
        config["caption"] = st.text_input("Caption", value=config.get("caption", ""), key=f"caption_{index}")
    
    elif section_type == "appendix":
        config["content"] = st.text_area("Appendix Content (Markdown)", value=config.get("content", ""), key=f"appendix_{index}", height=150)
    
    # Update config
    section["config"] = config


# ============================================================================
# TABS
# ============================================================================

tab1, tab2, tab3, tab4, tab5 = st.tabs([
    "🔨 Build Report",
    "📁 My Templates", 
    "👁️ Preview",
    "📥 Export",
    "📚 Section Library"
])


# ============================================================================
# TAB 1: BUILD REPORT
# ============================================================================

with tab1:
    col_library, col_canvas = st.columns([1, 2])
    
    with col_library:
        st.subheader("➕ Add Sections")
        
        # Group sections by type
        content_sections = ["title_page", "executive_summary", "free_text", "image", "appendix", "table_of_contents"]
        data_sections = ["compound_profile", "compound_table", "experiment_summary", "dataset_stats"]
        viz_sections = ["activity_chart", "admet_radar", "signature_heatmap", "pathway_enrichment"]
        
        st.markdown("**Content**")
        for section_type in content_sections:
            info = service.SECTION_TYPE_MAP.get(section_type, {})
            if st.button(f"{info.get('icon', '📄')} {info.get('name', section_type)}", key=f"add_{section_type}", use_container_width=True):
                add_section(section_type)
                st.rerun()
        
        st.markdown("**Data**")
        for section_type in data_sections:
            info = service.SECTION_TYPE_MAP.get(section_type, {})
            if st.button(f"{info.get('icon', '📊')} {info.get('name', section_type)}", key=f"add_{section_type}", use_container_width=True):
                add_section(section_type)
                st.rerun()
        
        st.markdown("**Visualizations**")
        for section_type in viz_sections:
            info = service.SECTION_TYPE_MAP.get(section_type, {})
            if st.button(f"{info.get('icon', '📈')} {info.get('name', section_type)}", key=f"add_{section_type}", use_container_width=True):
                add_section(section_type)
                st.rerun()
    
    with col_canvas:
        st.subheader("📋 Report Sections")
        
        if not st.session_state.report_sections:
            st.info("Add sections from the left panel to build your report.")
        else:
            for i, section in enumerate(st.session_state.report_sections):
                section_info = service.SECTION_TYPE_MAP.get(section["type"], {})
                
                with st.expander(f"{i+1}. {section_info.get('icon', '📄')} {section_info.get('name', section['type'])}", expanded=True):
                    # Section config based on type
                    render_section_config(section, i)
                    
                    # Actions
                    col1, col2, col3, col4 = st.columns(4)
                    with col1:
                        if st.button("⬆️", key=f"up_{i}", disabled=i==0):
                            move_section(i, i-1)
                            st.rerun()
                    with col2:
                        if st.button("⬇️", key=f"down_{i}", disabled=i==len(st.session_state.report_sections)-1):
                            move_section(i, i+1)
                            st.rerun()
                    with col3:
                        if st.button("🗑️", key=f"del_{i}"):
                            remove_section(i)
                            st.rerun()
            
            st.markdown("---")
            st.caption(f"Total sections: {len(st.session_state.report_sections)}")


# ============================================================================
# TAB 2: MY TEMPLATES
# ============================================================================

with tab2:
    st.subheader("📁 Saved Templates")
    
    user = get_current_user()
    user_id = user.get("id") if user else None
    
    col1, col2 = st.columns([3, 1])
    with col2:
        if st.button("💾 Save Current as Template", type="primary", disabled=not st.session_state.report_sections):
            st.session_state.show_save_dialog = True
    
    # Save dialog
    if st.session_state.get("show_save_dialog"):
        with st.form("save_template_form"):
            st.subheader("Save Template")
            template_name = st.text_input("Template Name *")
            template_desc = st.text_area("Description")
            is_public = st.checkbox("Make public (visible to team)")
            
            programs = get_programs()
            prog_options = [("", "-- No Program --")] + programs
            program_id = st.selectbox(
                "Associate with Program",
                options=[p[0] for p in prog_options],
                format_func=lambda x: dict(prog_options).get(x, x)
            )
            
            col1, col2 = st.columns(2)
            with col1:
                if st.form_submit_button("Save", type="primary"):
                    if template_name:
                        try:
                            with db_session() as db:
                                service.create_template(
                                    db=db,
                                    name=template_name,
                                    description=template_desc if template_desc else None,
                                    sections=st.session_state.report_sections,
                                    is_public=is_public,
                                    program_id=UUID(program_id) if program_id else None,
                                    created_by_id=UUID(user_id) if user_id else None,
                                )
                            st.success(f"✅ Template '{template_name}' saved!")
                            st.session_state.show_save_dialog = False
                            st.rerun()
                        except ValueError as e:
                            st.error(str(e))
                    else:
                        st.error("Template name is required")
            with col2:
                if st.form_submit_button("Cancel"):
                    st.session_state.show_save_dialog = False
                    st.rerun()
    
    # List templates
    with db_session() as db:
        templates = service.list_templates(
            db=db,
            created_by_id=UUID(user_id) if user_id else None,
            include_public=True,
            limit=50,
        )
        
        if templates:
            for template in templates:
                with st.expander(f"{'🌐' if template.is_public else '🔒'} {template.name}"):
                    st.write(template.description or "_No description_")
                    st.caption(f"Sections: {len(template.sections)} | Created: {template.created_at.strftime('%Y-%m-%d') if template.created_at else '-'}")
                    
                    col1, col2, col3 = st.columns(3)
                    with col1:
                        if st.button("📂 Load", key=f"load_{template.id}"):
                            st.session_state.report_sections = template.sections.copy()
                            st.success(f"Loaded '{template.name}'")
                            st.rerun()
                    with col2:
                        if st.button("📋 Clone", key=f"clone_{template.id}"):
                            try:
                                with db_session() as db2:
                                    service.clone_template(
                                        db=db2,
                                        template_id=template.id,
                                        new_name=f"{template.name} (Copy)",
                                        created_by_id=UUID(user_id) if user_id else None,
                                    )
                                st.success("Template cloned!")
                                st.rerun()
                            except ValueError as e:
                                st.error(str(e))
                    with col3:
                        if st.button("🗑️ Delete", key=f"delete_{template.id}"):
                            with db_session() as db2:
                                service.delete_template(db2, template.id)
                            st.success("Template deleted")
                            st.rerun()
        else:
            st.info("No templates yet. Build a report and save it as a template!")


# ============================================================================
# TAB 3: PREVIEW
# ============================================================================

with tab3:
    st.subheader("👁️ Report Preview")
    
    if st.button("🔄 Generate Preview", type="primary", disabled=not st.session_state.report_sections):
        with db_session() as db:
            html = service.build_report(
                st.session_state.report_sections,
                db,
                title="Report Preview"
            )
            st.session_state.preview_html = html
    
    if st.session_state.preview_html:
        st.components.v1.html(st.session_state.preview_html, height=800, scrolling=True)
    else:
        st.info("Click 'Generate Preview' to see your report.")


# ============================================================================
# TAB 4: EXPORT
# ============================================================================

with tab4:
    st.subheader("📥 Export Report")
    
    if not st.session_state.report_sections:
        st.warning("Add sections to your report before exporting.")
    else:
        export_title = st.text_input("Report Title", value="Custom Report")
        export_format = st.radio("Format", ["HTML", "PDF"], horizontal=True)
        
        if st.button("📥 Generate & Download", type="primary"):
            with db_session() as db:
                html = service.build_report(
                    st.session_state.report_sections,
                    db,
                    title=export_title
                )
                
                if export_format == "HTML":
                    st.download_button(
                        "⬇️ Download HTML",
                        data=html,
                        file_name=f"{export_title.replace(' ', '_')}.html",
                        mime="text/html",
                    )
                else:
                    try:
                        pdf_bytes = service.export_to_pdf(html)
                        st.download_button(
                            "⬇️ Download PDF",
                            data=pdf_bytes,
                            file_name=f"{export_title.replace(' ', '_')}.pdf",
                            mime="application/pdf",
                        )
                    except Exception as e:
                        st.error(f"PDF export failed: {e}")
                        st.info("Try downloading as HTML instead.")


# ============================================================================
# TAB 5: SECTION LIBRARY
# ============================================================================

with tab5:
    st.subheader("📚 Available Section Types")
    st.markdown("Browse all available section types you can add to your reports.")
    
    for section in service.SECTION_REGISTRY:
        with st.expander(f"{section['icon']} {section['name']}"):
            st.write(section["description"])
            if section["requires_entity"]:
                st.caption(f"📌 Requires: {section.get('entity_type', 'entity')}")
            else:
                st.caption("📌 No entity required")
            
            if st.button("➕ Add to Report", key=f"lib_add_{section['type']}"):
                add_section(section["type"])
                st.success(f"Added {section['name']} to report")
                st.rerun()


# ============================================================================
# REGISTER IN PAGE_REGISTRY
# ============================================================================

def main():
    """Entry point for page registry."""
    return  # Page renders on import


# Add to PAGE_REGISTRY in scripts/dashboard/core/__init__.py:
# "Report Builder": {"module": "scripts.dashboard.pages.report_builder", "icon": "📄", "group": "Collaboration"}
