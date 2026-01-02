"""Variant Tracking page for genomic variants."""

from __future__ import annotations

import streamlit as st
import pandas as pd
from amprenta_rag.models.misc import GeneticVariant
from amprenta_rag.database.session import db_session


def _get_auth_headers() -> dict:
    """Get authorization headers for API calls."""
    token = st.session_state.get("auth_token")
    return {"Authorization": f"Bearer {token}"} if token else {}


def _get_annotation_stats() -> dict:
    """Fetch annotation statistics from API."""
    import httpx
    import os
    API_BASE = os.environ.get("API_URL", "http://localhost:8000")
    try:
        with httpx.Client(timeout=30) as client:
            response = client.get(
                f"{API_BASE}/api/v1/genomics/variants/annotations/stats",
                headers=_get_auth_headers()
            )
            if response.status_code == 200:
                return response.json()
    except Exception:
        pass
    return {}


def _get_variant_annotations(variant_id: str) -> list:
    """Fetch annotations for a specific variant."""
    import httpx
    import os
    API_BASE = os.environ.get("API_URL", "http://localhost:8000")
    try:
        with httpx.Client(timeout=30) as client:
            response = client.get(
                f"{API_BASE}/api/v1/genomics/variants/{variant_id}/annotations",
                headers=_get_auth_headers()
            )
            if response.status_code == 200:
                return response.json()
    except Exception:
        pass
    return []


def _trigger_batch_annotation(variant_ids: list[str]) -> dict:
    """Trigger batch annotation job."""
    import httpx
    import os
    API_BASE = os.environ.get("API_URL", "http://localhost:8000")
    try:
        with httpx.Client(timeout=30) as client:
            response = client.post(
                f"{API_BASE}/api/v1/genomics/variants/annotate/batch",
                headers=_get_auth_headers(),
                json={"variant_ids": variant_ids}
            )
            if response.status_code == 200:
                return response.json()
    except Exception:
        pass
    return {}


def _render_impact_badge(impact: str) -> None:
    """Render color-coded impact badge."""
    colors = {
        "HIGH": ("🔴", "error"),
        "MODERATE": ("🟠", "warning"),
        "LOW": ("🟡", "info"),
        "MODIFIER": ("🟢", "success"),
    }
    emoji, style = colors.get(impact, ("⚪", "info"))
    st.markdown(f"{emoji} **{impact}**")


def _render_prediction_badge(pred_type: str, prediction: str, score: float = None) -> None:
    """Render SIFT/PolyPhen prediction badge."""
    if not prediction:
        return
    
    # Color based on prediction
    if prediction.lower() in ["deleterious", "probably_damaging", "possibly_damaging"]:
        color = "🔴"
    elif prediction.lower() in ["tolerated", "benign"]:
        color = "🟢"
    else:
        color = "🟡"
    
    score_text = f" ({score:.3f})" if score is not None else ""
    st.markdown(f"{color} **{pred_type}**: {prediction}{score_text}")


def render_variants_page() -> None:
    """Render the Variant Tracking page."""
    st.title("🧬 Variant Tracking")
    st.markdown("Track genetic variants in experiments and cell lines")

    tab_browse, tab_add, tab_vcf, tab_annotations = st.tabs(["Browse", "Add Variant", "📤 Upload VCF", "🔬 Annotations"])

    with tab_browse:
        _render_browse_tab()

    with tab_add:
        _render_add_tab()

    with tab_vcf:
        _render_vcf_upload_tab()
    
    with tab_annotations:
        _render_annotations_tab()


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
            
            # Batch annotation section
            st.divider()
            st.subheader("🔬 Batch Annotation")
            st.markdown("Select variants to annotate with VEP (Variant Effect Predictor)")
            
            # Check which variants don't have annotations
            with db_session() as db:
                # Get variant IDs that need annotation (this is simplified - in reality we'd need to map GeneticVariant to Variant)
                # For now, just show a sample selection interface
                variant_ids_for_annotation = [str(v.id) for v in variants[:10]]  # Sample first 10
                
                if variant_ids_for_annotation:
                    # Multi-select for variants to annotate
                    col1, col2 = st.columns([3, 1])
                    
                    with col1:
                        selected_variants = st.multiselect(
                            "Select variants to annotate:",
                            options=variant_ids_for_annotation,
                            format_func=lambda x: f"Variant {x[:8]}...",  # Show first 8 chars of UUID
                            key="batch_annotation_selection"
                        )
                    
                    with col2:
                        st.markdown("**Annotation Info**")
                        st.info("💡 Max 1000 variants per batch")
                        if selected_variants:
                            st.success(f"✅ {len(selected_variants)} selected")
                    
                    # Annotation trigger button
                    if selected_variants:
                        if st.button("🚀 Start Batch Annotation", type="primary", key="trigger_batch_annotation"):
                            with st.spinner("Submitting annotation job..."):
                                result = _trigger_batch_annotation(selected_variants)
                                
                                if result.get("success"):
                                    st.success("✅ Annotation job queued successfully!")
                                    st.info(f"**Job ID:** {result.get('job_id')}")
                                    st.info(f"**Queued variants:** {result.get('queued_count')}")
                                    st.markdown(result.get("message", ""))
                                    st.balloons()
                                else:
                                    st.error("❌ Failed to submit annotation job. Check API connectivity.")
                    else:
                        st.info("Select variants above to enable batch annotation.")
                else:
                    st.info("No variants available for annotation.")
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


def _render_annotations_tab() -> None:
    """Render the Annotations tab with VEP results."""
    st.subheader("🔬 Variant Annotations")
    
    # Stats display
    stats = _get_annotation_stats()
    if stats:
        st.markdown("### 📊 Annotation Statistics")
        col1, col2, col3, col4 = st.columns(4)
        
        with col1:
            st.metric("Total Variants", stats.get("total_variants", 0))
        with col2:
            st.metric("Annotated", stats.get("annotated_variants", 0))
        with col3:
            st.metric("Missing", stats.get("unannotated_variants", 0))
        with col4:
            st.metric("High Impact", stats.get("high_impact", 0))
        
        # Impact breakdown
        st.markdown("### Impact Distribution")
        col1, col2, col3 = st.columns(3)
        with col1:
            st.metric("🔴 HIGH", stats.get("high_impact", 0))
        with col2:
            st.metric("🟠 MODERATE", stats.get("moderate_impact", 0))
        with col3:
            st.metric("🟡 LOW", stats.get("low_impact", 0))
    else:
        st.warning("Unable to fetch annotation statistics. Check API connectivity.")
    
    st.divider()
    
    # Annotation browser
    st.markdown("### 🔍 Browse Annotations")
    
    # Get variants with annotations
    with db_session() as db:
        from amprenta_rag.database.models import Variant, VariantAnnotation
        
        # Query variants that have annotations
        variants_with_annotations = db.query(Variant).join(VariantAnnotation).distinct().limit(100).all()
        
        if variants_with_annotations:
            # Create variant options for dropdown
            variant_options = {}
            for v in variants_with_annotations:
                label = f"{v.gene_symbol or 'Unknown'} | {v.chromosome}:{v.position} {v.ref_allele}>{v.alt_allele}"
                variant_options[str(v.id)] = label
            
            # Variant selector
            selected_variant_id = st.selectbox(
                "Select variant to view annotations:",
                options=list(variant_options.keys()),
                format_func=lambda x: variant_options[x],
                key="annotation_variant_selector"
            )
            
            if selected_variant_id:
                # Fetch annotations for selected variant
                annotations = _get_variant_annotations(selected_variant_id)
                
                if annotations:
                    st.markdown("### 📋 Annotation Details")
                    
                    for i, ann in enumerate(annotations):
                        with st.expander(f"Annotation {i+1} - {ann.get('source', 'Unknown')} | {ann.get('annotated_at', 'Unknown date')}", expanded=i==0):
                            col1, col2 = st.columns(2)
                            
                            with col1:
                                st.markdown("**Basic Information**")
                                if ann.get("impact"):
                                    st.markdown("**Impact:**")
                                    _render_impact_badge(ann["impact"])
                                
                                if ann.get("consequence"):
                                    st.markdown(f"**Consequence:** {ann['consequence']}")
                                
                                if ann.get("symbol"):
                                    st.markdown(f"**Gene:** {ann['symbol']}")
                                
                                if ann.get("gene_id"):
                                    st.markdown(f"**Gene ID:** {ann['gene_id']}")
                            
                            with col2:
                                st.markdown("**Predictions**")
                                if ann.get("sift_prediction"):
                                    _render_prediction_badge("SIFT", ann["sift_prediction"], ann.get("sift_score"))
                                
                                if ann.get("polyphen_prediction"):
                                    _render_prediction_badge("PolyPhen", ann["polyphen_prediction"], ann.get("polyphen_score"))
                                
                                if ann.get("clin_sig"):
                                    st.markdown(f"**Clinical Significance:** {ann['clin_sig']}")
                            
                            # Additional details
                            if ann.get("transcript_id") or ann.get("amino_acids") or ann.get("codons"):
                                st.markdown("**Transcript Details**")
                                if ann.get("transcript_id"):
                                    st.markdown(f"**Transcript:** {ann['transcript_id']}")
                                if ann.get("amino_acids"):
                                    st.markdown(f"**Amino Acid Change:** {ann['amino_acids']}")
                                if ann.get("codons"):
                                    st.markdown(f"**Codon Change:** {ann['codons']}")
                                if ann.get("protein_position"):
                                    st.markdown(f"**Protein Position:** {ann['protein_position']}")
                else:
                    st.info("No annotations found for this variant.")
        else:
            st.info("No annotated variants found. Upload variants and run annotation jobs to see results here.")

