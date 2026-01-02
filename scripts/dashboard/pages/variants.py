"""Variant Tracking page for genomic variants."""

from __future__ import annotations

import streamlit as st
import pandas as pd
from amprenta_rag.database.models import GeneticVariant
from amprenta_rag.database.session import db_session


def render_variants_page() -> None:
    """Render the Variant Tracking page."""
    st.title("🧬 Variant Tracking")
    st.markdown("Track genetic variants in experiments and cell lines")

    tab_browse, tab_add, tab_vcf = st.tabs(["Browse", "Add Variant", "📤 Upload VCF"])

    with tab_browse:
        _render_browse_tab()

    with tab_add:
        _render_add_tab()

    with tab_vcf:
        _render_vcf_upload_tab()


def _render_browse_tab() -> None:
    """Render the Browse tab."""
    st.subheader("Browse Variants")

    # Filters
    col1, col2 = st.columns(2)
    with col1:
        gene_filter = st.text_input("Gene", placeholder="Filter by gene name", key="browse_gene_filter")
    with col2:
        organism_filter = st.text_input("Cell Line/Organism", placeholder="Filter by organism", key="browse_organism_filter")

    # Query variants
    with db_session() as db:
        query = db.query(GeneticVariant)

        if gene_filter:
            query = query.filter(GeneticVariant.gene.ilike(f"%{gene_filter}%"))
        if organism_filter:
            query = query.filter(GeneticVariant.organism.ilike(f"%{organism_filter}%"))

        variants = query.order_by(GeneticVariant.created_at.desc()).limit(500).all()

        if variants:
            st.metric("Variants Found", len(variants))
            if len(variants) == 500:
                st.info("Showing first 500 variants. Use filters to narrow down.")

            # Prepare data for display
            variant_data = []
            for v in variants:
                variant_data.append({
                    "Gene": v.gene,
                    "Variant": v.variant,
                    "Zygosity": v.zygosity or "N/A",
                    "Organism": v.organism,
                    "Experiment": v.experiment.name if v.experiment else "N/A",
                    "Notes": v.notes[:100] + "..." if v.notes and len(v.notes) > 100 else (v.notes or ""),
                    "Created": v.created_at.strftime("%Y-%m-%d %H:%M") if v.created_at else "N/A",
                })

            df = pd.DataFrame(variant_data)
            st.dataframe(df, use_container_width=True, hide_index=True)
        else:
            st.info("No variants found. Add variants using the 'Add Variant' tab.")


def _render_add_tab() -> None:
    """Render the Add Variant tab."""
    st.subheader("Add Variant")

    with st.form("add_variant_form"):
        gene = st.text_input("Gene*", placeholder="e.g., TP53, BRCA1")
        variant = st.text_input("Variant*", placeholder="e.g., p.R273H, c.742C>T")
        zygosity = st.selectbox(
            "Zygosity",
            ["", "Homozygous", "Heterozygous", "Hemizygous"],
            index=0,
        )
        organism = st.text_input("Cell Line/Organism*", placeholder="e.g., HeLa, HEK293, mouse")
        notes = st.text_area("Notes", placeholder="Additional information about this variant")

        submitted = st.form_submit_button("💾 Save Variant", type="primary")

        if submitted:
            # Validation
            if not gene or not gene.strip():
                st.error("Gene is required")
                return
            if not variant or not variant.strip():
                st.error("Variant is required")
                return
            if not organism or not organism.strip():
                st.error("Cell Line/Organism is required")
                return

            # Get experiment if needed (optional)
            experiment_id = None
            # Could add experiment selector here if needed

            # Save variant
            try:
                with db_session() as db:
                    new_variant = GeneticVariant(
                        gene=gene.strip(),
                        variant=variant.strip(),
                        zygosity=zygosity.strip() if zygosity else None,
                        organism=organism.strip(),
                        experiment_id=experiment_id,
                        notes=notes.strip() if notes else None,
                    )
                    db.add(new_variant)
                    db.commit()
                    st.success(f"Variant {gene} {variant} saved successfully!")
                    st.rerun()
            except Exception as e:
                st.error(f"Failed to save variant: {e}")
                st.exception(e)


def _render_vcf_upload_tab() -> None:
    """VCF file upload and import tab."""
    st.subheader("📤 Upload VCF File")
    st.caption("Import genetic variants from VCF files")
    
    # File uploader
    uploaded_file = st.file_uploader(
        "Choose VCF file",
        type=["vcf"],  # Note: .vcf.gz may need special handling
        help="Max 50MB. Supports standard VCF format.",
        key="vcf_file_uploader"
    )
    
    # Optional experiment selector
    col1, col2 = st.columns([2, 1])
    
    with col1:
        # Get experiments for selector
        from amprenta_rag.database.session import db_session
        from amprenta_rag.database.models import Experiment
        
        experiments = []
        with db_session() as db:
            exps = db.query(Experiment).order_by(Experiment.name).limit(100).all()
            experiments = [(str(e.id), e.name) for e in exps]
        
        experiment_options = [("", "No experiment (standalone)")] + experiments
        experiment_id = st.selectbox(
            "Link to Experiment (optional)",
            options=[x[0] for x in experiment_options],
            format_func=lambda x: dict(experiment_options).get(x, x),
            key="vcf_experiment"
        )
    
    with col2:
        preview_only = st.checkbox("Preview only", value=True, key="vcf_preview_mode")
    
    if uploaded_file:
        st.info(f"📁 File: {uploaded_file.name} ({uploaded_file.size / 1024:.1f} KB)")
        
        action_label = "🔍 Preview Variants" if preview_only else "📥 Import Variants"
        
        if st.button(action_label, key="vcf_action_btn", type="primary"):
            with st.spinner("Processing VCF file..."):
                import httpx
                import os
                
                API_BASE = os.environ.get("API_URL", "http://localhost:8000")
                
                try:
                    params = f"?preview_only={str(preview_only).lower()}"
                    if experiment_id:
                        params += f"&experiment_id={experiment_id}"
                    
                    with httpx.Client(timeout=120) as client:
                        files = {"file": (uploaded_file.name, uploaded_file.getvalue(), "text/plain")}
                        response = client.post(
                            f"{API_BASE}/api/v1/genomics/variants/upload{params}",
                            files=files,
                        )
                        result = response.json()
                    
                    if result.get("success"):
                        if preview_only:
                            st.success(f"✅ Parsed {result.get('variants_count', 0)} variants")
                            
                            # Show preview table
                            variants = result.get("variants", [])
                            if variants:
                                import pandas as pd
                                df = pd.DataFrame(variants)
                                st.dataframe(df, use_container_width=True, hide_index=True)
                                
                                if result.get("variants_count", 0) > len(variants):
                                    st.info(f"Showing first {len(variants)} of {result.get('variants_count')} variants")
                        else:
                            st.success(f"✅ Imported {result.get('imported_count', 0)} variants to database")
                            st.balloons()
                    else:
                        st.error(f"❌ Error: {result.get('error', 'Unknown error')}")
                
                except Exception as e:
                    st.error(f"❌ Upload failed: {str(e)}")

