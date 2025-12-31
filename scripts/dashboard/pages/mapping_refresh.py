"""ID Mapping Refresh Dashboard - Monitor and manage ID mappings."""

from __future__ import annotations

import streamlit as st


def render_mapping_refresh_page() -> None:
    """Render the ID Mapping Refresh management page."""
    try:
        from scripts.dashboard.core.auth import require_auth
        user = require_auth()
    except Exception:
        # Fallback for testing/development
        user = {"id": "test", "username": "test", "role": "admin"}
    
    st.header("🔄 ID Mapping Refresh")
    st.markdown("Monitor UniProt/KEGG ID mappings and trigger manual refreshes.")
    
    # 4 tabs as specified in plan
    tab1, tab2, tab3, tab4 = st.tabs(["Status", "Statistics", "Lookup", "Jobs"])
    
    with tab1:
        render_status_tab(user)
    with tab2:
        render_statistics_tab()
    with tab3:
        render_lookup_tab()
    with tab4:
        render_jobs_tab(user)


def render_status_tab(user: dict) -> None:
    """Display mapping status and refresh controls."""
    st.subheader("📊 Mapping Status")
    st.info("Demo mode - API endpoints not yet implemented")
    
    # Show demo data
    col1, col2, col3 = st.columns(3)
    with col1:
        st.metric("Total Mappings", "12,345")
    with col2:
        st.metric("Expired", "234")
    with col3:
        st.metric("Types", "3")
    
    st.subheader("Breakdown by Source Type")
    st.write("**Gene**: 5,432")
    st.write("**Protein**: 4,321") 
    st.write("**Metabolite**: 2,592")
    
    # Admin-only refresh button
    if user.get("role") == "admin":
        st.divider()
        st.subheader("🔧 Manual Refresh (Admin Only)")
        st.info("API endpoints not yet implemented. This is a demo interface.")
        
        if st.button("🔄 Refresh UniProt Mappings", key="refresh_uniprot", type="primary"):
            st.info("🚧 Demo mode: Refresh functionality will be available when API endpoints are implemented.")
    else:
        st.info("🔒 Admin access required to trigger manual refreshes.")


def render_statistics_tab() -> None:
    """Display detailed mapping statistics."""
    st.subheader("📈 Detailed Statistics")
    st.info("Demo mode - showing sample data")
    
    # Show demo data
    col1, col2, col3 = st.columns(3)
    with col1:
        st.metric("Total Mappings", "12,345")
    with col2:
        st.metric("Permanent Mappings", "8,901")
    with col3:
        st.metric("TTL Mappings", "3,444")
    
    # Demo charts
    st.subheader("📊 Coverage by Source Type")
    demo_source_data = {"gene": 5432, "protein": 4321, "metabolite": 2592}
    st.bar_chart(demo_source_data)
    
    st.subheader("🎯 Coverage by Target Type")
    demo_target_data = {"uniprot": 6543, "kegg_gene": 3456, "kegg_compound": 2346}
    st.bar_chart(demo_target_data)
    
    st.success("✅ No expired mappings found.")


def render_lookup_tab() -> None:
    """Interactive ID lookup functionality."""
    st.subheader("🔍 ID Lookup")
    st.info("Demo mode - lookup functionality will be available when API is implemented")
    
    # Single ID lookup section
    st.markdown("### Single ID Lookup")
    
    col1, col2 = st.columns(2)
    with col1:
        source_type = st.selectbox("Source Type", ["gene", "protein", "metabolite"], key="single_source_type")
    with col2:
        source_id = st.text_input("Source ID", placeholder="e.g., TP53, BRCA1", key="single_source_id")
    
    _fallback = st.checkbox("Enable API fallback", value=True, 
                           help="If disabled, only checks database. If enabled, falls back to external APIs if not found in database.")
    
    if st.button("🔍 Look Up", key="lookup_single") and source_id:
        # Demo lookup results
        if source_id.upper() == "TP53":
            st.success(f"✅ Found 2 mapping(s) for {source_type} '{source_id}' (demo):")
            st.write("**UNIPROT**: `P04637`")
            st.write("**KEGG_GENE**: `hsa:7157`")
        elif source_id.upper() in ["BRCA1", "EGFR", "MYC"]:
            st.success(f"✅ Found 1 mapping(s) for {source_type} '{source_id}' (demo):")
            st.write("**UNIPROT**: `P12345`")
        else:
            st.warning(f"⚠️ No mappings found for {source_type} '{source_id}' (demo)")
    
    # Batch lookup section
    st.divider()
    st.markdown("### Batch Lookup")
    st.info("Look up multiple IDs at once (maximum 1000 IDs per request)")
    
    col1, col2 = st.columns(2)
    with col1:
        _batch_source_type = st.selectbox("Source Type", ["gene", "protein", "metabolite"], key="batch_source_type")
    with col2:
        _batch_target_type = st.selectbox("Target Type", ["uniprot", "kegg_gene", "kegg_compound"], key="batch_target_type")
    
    ids_input = st.text_area("IDs (one per line)", placeholder="TP53\nBRCA1\nEGFR", height=150, key="batch_ids")
    
    if st.button("🔍 Batch Look Up", key="lookup_batch") and ids_input:
        ids = [id_val.strip() for id_val in ids_input.split("\n") if id_val.strip()]
        
        if len(ids) > 1000:
            st.error("❌ Maximum 1000 IDs per request. Please reduce the number of IDs.")
        elif len(ids) == 0:
            st.error("❌ Please enter at least one ID.")
        else:
            # Demo batch results
            known_genes = {"TP53", "BRCA1", "EGFR", "MYC", "AKT1", "PIK3CA"}
            demo_found = [id_val for id_val in ids if id_val.upper() in known_genes]
            demo_not_found = [id_val for id_val in ids if id_val.upper() not in known_genes]
            
            st.success(f"✅ Batch lookup complete (demo): {len(demo_found)} found, {len(demo_not_found)} not found")
            
            col1, col2 = st.columns(2)
            with col1:
                st.metric("Found", len(demo_found))
            with col2:
                st.metric("Not Found", len(demo_not_found))
            
            if demo_found or demo_not_found:
                st.subheader("Detailed Results (Demo)")
                
                if demo_found:
                    st.markdown("**Found Mappings:**")
                    for gene in demo_found:
                        st.write(f"• {gene} → `P{hash(gene.upper()) % 10000:05d}`")
                
                if demo_not_found:
                    st.markdown("**Not Found:**")
                    for gene in demo_not_found:
                        st.write(f"• {gene}")


def render_jobs_tab(user: dict) -> None:
    """Display recent refresh job history."""
    st.subheader("🔄 Refresh Jobs")
    st.info("Demo mode - job management will be available when API is implemented")
    
    # Show current job status info
    st.markdown("### Current Status")
    st.markdown("""
    - **Scheduled Jobs**: UniProt mappings refresh weekly (Sunday 2:00 AM)
    - **Cleanup Jobs**: Expired mappings cleanup daily (3:00 AM)
    - **Manual Jobs**: Triggered via Status tab (admin only)
    """)
    
    # Admin quick actions
    if user.get("role") == "admin":
        st.divider()
        st.markdown("### 🚀 Quick Actions (Admin)")
        
        col1, col2 = st.columns(2)
        
        with col1:
            if st.button("🔄 Trigger UniProt Refresh", key="jobs_refresh", type="primary"):
                st.info("🚧 Demo mode: API endpoints not yet implemented. Job functionality will be available when backend is complete.")
        
        with col2:
            if st.button("📊 Refresh Statistics", key="refresh_stats"):
                st.info("🔄 Statistics are refreshed automatically when you view them in other tabs.")
                st.rerun()
    else:
        st.info("🔒 Admin access required for manual job triggers.")
