"""Target Management Dashboard."""

import os
import uuid
from typing import Dict, Any, List, Optional

import pandas as pd
import streamlit as st
import requests
from datetime import datetime

# API Configuration
API_BASE = os.environ.get("API_URL", "http://localhost:8000")


def get_auth_headers() -> Dict[str, str]:
    """Get authentication headers for API requests."""
    # In production, this would use actual session tokens
    return {"Authorization": "Bearer dummy-token"}


def api_get(path: str) -> Dict[str, Any]:
    """Make authenticated GET request to API."""
    try:
        response = requests.get(f"{API_BASE}{path}", headers=get_auth_headers(), timeout=30.0)
        response.raise_for_status()
        return response.json()
    except Exception as e:
        st.error(f"API request failed: {str(e)}")
        return {}


def api_post(path: str, payload: Dict[str, Any]) -> Dict[str, Any]:
    """Make authenticated POST request to API."""
    try:
        response = requests.post(f"{API_BASE}{path}", json=payload, headers=get_auth_headers(), timeout=30.0)
        response.raise_for_status()
        return response.json()
    except Exception as e:
        st.error(f"API request failed: {str(e)}")
        return {}


def api_put(path: str, payload: Dict[str, Any]) -> Dict[str, Any]:
    """Make authenticated PUT request to API."""
    try:
        response = requests.put(f"{API_BASE}{path}", json=payload, headers=get_auth_headers(), timeout=30.0)
        response.raise_for_status()
        return response.json()
    except Exception as e:
        st.error(f"API request failed: {str(e)}")
        return {}


def get_druggability_color(score: Optional[float]) -> str:
    """Get color for druggability score."""
    if score is None:
        return "gray"
    elif score > 0.7:
        return "green"
    elif score >= 0.4:
        return "orange"
    else:
        return "red"


def render_target_list_tab():
    """Render Tab 1: Target List."""
    st.subheader("Target List")
    
    # Filters
    col1, col2, col3, col4 = st.columns(4)
    
    with col1:
        target_class = st.selectbox(
            "Target Class",
            ["All", "kinase", "gpcr", "ion_channel", "enzyme", "nuclear_receptor", "other"],
            index=0
        )
    
    with col2:
        target_family = st.text_input("Target Family", placeholder="e.g., Tyrosine kinase")
    
    with col3:
        lifecycle_status = st.selectbox(
            "Status",
            ["active", "archived", "deprecated"],
            index=0
        )
    
    with col4:
        search_query = st.text_input("Search", placeholder="Target name, gene symbol...")
    
    # Search vs List
    if search_query:
        with st.spinner("Searching targets..."):
            targets_data = api_get(f"/api/v1/targets/search/{search_query}")
            if isinstance(targets_data, list):
                targets = targets_data
            else:
                targets = []
    else:
        # Build query parameters
        params = {"lifecycle_status": lifecycle_status}
        if target_class != "All":
            params["target_class"] = target_class
        if target_family:
            params["target_family"] = target_family
        
        query_string = "&".join([f"{k}={v}" for k, v in params.items()])
        
        with st.spinner("Loading targets..."):
            targets = api_get(f"/api/v1/targets?{query_string}")
    
    if not targets:
        st.info("No targets found. Create your first target below.")
    else:
        # Convert to DataFrame for display
        df_data = []
        for target in targets:
            df_data.append({
                "Name": target["name"],
                "Gene Symbol": target.get("gene_symbol", ""),
                "Class": target.get("target_class", ""),
                "Family": target.get("target_family", ""),
                "Druggability": target.get("druggability_score", 0.0),
                "Status": target["lifecycle_status"],
                "ID": target["id"]
            })
        
        df = pd.DataFrame(df_data)
        
        # Display metrics
        col1, col2, col3, col4 = st.columns(4)
        with col1:
            st.metric("Total Targets", len(targets))
        with col2:
            validated_count = sum(1 for t in targets if t.get("validation_status") in ["validated", "clinical"])
            st.metric("Validated", validated_count)
        with col3:
            druggable_count = sum(1 for t in targets if (t.get("druggability_score") or 0) > 0.7)
            st.metric("Highly Druggable", druggable_count)
        with col4:
            kinase_count = sum(1 for t in targets if t.get("target_class") == "kinase")
            st.metric("Kinases", kinase_count)
        
        # Display table
        st.dataframe(
            df.drop("ID", axis=1),  # Hide ID column
            use_container_width=True,
            hide_index=True
        )
        
        # Store targets in session state for Tab 2
        st.session_state["targets_list"] = targets
    
    # Create Target Form
    st.divider()
    st.subheader("Create New Target")
    
    with st.form("create_target_form"):
        col1, col2 = st.columns(2)
        
        with col1:
            name = st.text_input("Target Name*", placeholder="e.g., EGFR")
            gene_symbol = st.text_input("Gene Symbol", placeholder="e.g., EGFR")
            target_class = st.selectbox(
                "Target Class",
                ["", "kinase", "gpcr", "ion_channel", "enzyme", "nuclear_receptor", "transporter", "other"]
            )
        
        with col2:
            description = st.text_area("Description", placeholder="Target description...")
            uniprot_id = st.text_input("UniProt ID", placeholder="e.g., P00533")
            target_family = st.text_input("Target Family", placeholder="e.g., Tyrosine kinase")
        
        submitted = st.form_submit_button("Create Target", type="primary")
        
        if submitted:
            if not name:
                st.error("Target name is required")
            else:
                payload = {
                    "name": name,
                    "description": description if description else None,
                    "gene_symbol": gene_symbol if gene_symbol else None,
                    "uniprot_id": uniprot_id if uniprot_id else None,
                    "target_class": target_class if target_class else None,
                    "target_family": target_family if target_family else None
                }
                
                with st.spinner("Creating target..."):
                    result = api_post("/api/v1/targets", payload)
                    if result:
                        st.success(f"Created target: {result['name']}")
                        st.rerun()
                    else:
                        st.error("Failed to create target")


def render_target_detail_tab():
    """Render Tab 2: Target Detail."""
    st.subheader("Target Details")
    
    # Target selection
    targets = st.session_state.get("targets_list", [])
    if not targets:
        st.info("No targets available. Please load targets in the Target List tab first.")
        return
    
    target_options = {f"{t['name']} ({t.get('gene_symbol', 'N/A')})": t["id"] for t in targets}
    selected_target_key = st.selectbox("Select Target", list(target_options.keys()))
    
    if not selected_target_key:
        return
    
    target_id = target_options[selected_target_key]
    
    # Load target details
    with st.spinner("Loading target details..."):
        target = api_get(f"/api/v1/targets/{target_id}")
    
    if not target:
        st.error("Failed to load target details")
        return
    
    # Display target information
    col1, col2 = st.columns(2)
    
    with col1:
        st.markdown("#### Basic Information")
        st.text(f"Name: {target['name']}")
        st.text(f"Gene Symbol: {target.get('gene_symbol', 'N/A')}")
        st.text(f"UniProt ID: {target.get('uniprot_id', 'N/A')}")
        st.text(f"Description: {target.get('description', 'N/A')}")
        
        st.markdown("#### Classification")
        st.text(f"Class: {target.get('target_class', 'N/A')}")
        st.text(f"Family: {target.get('target_family', 'N/A')}")
        st.text(f"Synthetic: {'Yes' if target.get('is_synthetic') else 'No'}")
    
    with col2:
        st.markdown("#### Druggability")
        score = target.get('druggability_score')
        if score is not None:
            color = get_druggability_color(score)
            st.metric("Druggability Score", f"{score:.2f}", delta_color=color)
        else:
            st.metric("Druggability Score", "Not calculated")
        
        st.text(f"Source: {target.get('druggability_source', 'N/A')}")
        st.text(f"Pocket Count: {target.get('pocket_count', 'N/A')}")
        
        st.markdown("#### Validation")
        st.text(f"Status: {target.get('validation_status', 'N/A')}")
        st.text(f"Lifecycle: {target['lifecycle_status']}")
    
    # External Links
    st.markdown("#### External Resources")
    col1, col2, col3 = st.columns(3)
    
    with col1:
        if target.get('uniprot_id'):
            st.link_button(
                "UniProt",
                f"https://www.uniprot.org/uniprot/{target['uniprot_id']}"
            )
    
    with col2:
        if target.get('chembl_id'):
            st.link_button(
                "ChEMBL",
                f"https://www.ebi.ac.uk/chembl/target_report_card/{target['chembl_id']}"
            )
    
    with col3:
        if target.get('ensembl_id'):
            st.link_button(
                "Ensembl",
                f"https://www.ensembl.org/Homo_sapiens/Gene/Summary?g={target['ensembl_id']}"
            )
    
    # Action buttons
    st.divider()
    col1, col2 = st.columns(2)
    
    with col1:
        if st.button("Calculate Druggability", type="primary"):
            with st.spinner("Calculating druggability..."):
                result = api_post(f"/api/v1/targets/{target_id}/druggability", {})
                if result:
                    st.success(f"Druggability calculated: {result['score']:.2f}")
                    if result.get('factors'):
                        st.json(result['factors'])
                    st.rerun()
    
    with col2:
        if st.button("Edit Target"):
            st.session_state["edit_target_id"] = target_id
            st.session_state["edit_target_data"] = target
            st.rerun()
    
    # Edit form (if in edit mode)
    if st.session_state.get("edit_target_id") == target_id:
        st.divider()
        st.subheader("Edit Target")
        
        with st.form("edit_target_form"):
            edit_data = st.session_state.get("edit_target_data", {})
            
            col1, col2 = st.columns(2)
            
            with col1:
                new_name = st.text_input("Name", value=edit_data.get("name", ""))
                new_gene_symbol = st.text_input("Gene Symbol", value=edit_data.get("gene_symbol", "") or "")
                new_target_class = st.selectbox(
                    "Target Class",
                    ["", "kinase", "gpcr", "ion_channel", "enzyme", "nuclear_receptor", "transporter", "other"],
                    index=0 if not edit_data.get("target_class") else 
                          ["", "kinase", "gpcr", "ion_channel", "enzyme", "nuclear_receptor", "transporter", "other"].index(edit_data.get("target_class"))
                )
                new_validation_status = st.selectbox(
                    "Validation Status",
                    ["", "hypothesis", "validated", "clinical"],
                    index=0 if not edit_data.get("validation_status") else
                          ["", "hypothesis", "validated", "clinical"].index(edit_data.get("validation_status"))
                )
            
            with col2:
                new_description = st.text_area("Description", value=edit_data.get("description", "") or "")
                new_target_family = st.text_input("Target Family", value=edit_data.get("target_family", "") or "")
                new_lifecycle_status = st.selectbox(
                    "Lifecycle Status",
                    ["active", "archived", "deprecated"],
                    index=["active", "archived", "deprecated"].index(edit_data.get("lifecycle_status", "active"))
                )
                new_is_synthetic = st.checkbox("Synthetic Target", value=edit_data.get("is_synthetic", False))
            
            col1, col2 = st.columns(2)
            with col1:
                update_submitted = st.form_submit_button("Update Target", type="primary")
            with col2:
                cancel_edit = st.form_submit_button("Cancel")
            
            if cancel_edit:
                del st.session_state["edit_target_id"]
                del st.session_state["edit_target_data"]
                st.rerun()
            
            if update_submitted:
                update_payload = {}
                if new_name != edit_data.get("name"):
                    update_payload["name"] = new_name
                if new_description != (edit_data.get("description") or ""):
                    update_payload["description"] = new_description or None
                if new_gene_symbol != (edit_data.get("gene_symbol") or ""):
                    update_payload["gene_symbol"] = new_gene_symbol or None
                if new_target_class != (edit_data.get("target_class") or ""):
                    update_payload["target_class"] = new_target_class or None
                if new_target_family != (edit_data.get("target_family") or ""):
                    update_payload["target_family"] = new_target_family or None
                if new_validation_status != (edit_data.get("validation_status") or ""):
                    update_payload["validation_status"] = new_validation_status or None
                if new_lifecycle_status != edit_data.get("lifecycle_status"):
                    update_payload["lifecycle_status"] = new_lifecycle_status
                if new_is_synthetic != edit_data.get("is_synthetic", False):
                    update_payload["is_synthetic"] = new_is_synthetic
                
                if update_payload:
                    with st.spinner("Updating target..."):
                        result = api_put(f"/api/v1/targets/{target_id}", update_payload)
                        if result:
                            st.success("Target updated successfully")
                            del st.session_state["edit_target_id"]
                            del st.session_state["edit_target_data"]
                            st.rerun()
                else:
                    st.info("No changes detected")


def render_competitive_landscape_tab():
    """Render Tab 3: Competitive Landscape."""
    st.subheader("Competitive Landscape")
    
    # Target selection
    targets = st.session_state.get("targets_list", [])
    if not targets:
        st.info("No targets available. Please load targets in the Target List tab first.")
        return
    
    target_options = {f"{t['name']} ({t.get('gene_symbol', 'N/A')})": t["id"] for t in targets}
    selected_target_key = st.selectbox("Select Target for Landscape Analysis", list(target_options.keys()))
    
    if not selected_target_key:
        return
    
    target_id = target_options[selected_target_key]
    
    # Load landscape data
    with st.spinner("Loading competitive landscape..."):
        landscape = api_get(f"/api/v1/targets/{target_id}/landscape")
    
    if landscape:
        st.info("🚧 Competitive landscape analysis coming soon!")
        
        # Show placeholder metrics
        col1, col2, col3, col4 = st.columns(4)
        
        with col1:
            st.metric("Published Compounds", landscape.get("published_compounds", 0))
        with col2:
            st.metric("Clinical Trials", landscape.get("clinical_trials", 0))
        with col3:
            st.metric("Approved Drugs", landscape.get("approved_drugs", 0))
        with col4:
            st.metric("Patent Families", 0)
        
        st.markdown("### Planned Features")
        st.markdown("""
        - **Literature Mining**: PubMed articles mentioning target
        - **Clinical Trials**: ClinicalTrials.gov integration
        - **Patent Landscape**: Patent family analysis
        - **Market Intelligence**: Competitive drug analysis
        - **Key Players**: Companies working on target
        """)
    else:
        st.error("Failed to load competitive landscape data")


def render_assay_linkage_tab():
    """Render Tab 4: Assay Linkage."""
    st.subheader("Assay Linkage")
    
    # Target selection
    targets = st.session_state.get("targets_list", [])
    if not targets:
        st.info("No targets available. Please load targets in the Target List tab first.")
        return
    
    target_options = {f"{t['name']} ({t.get('gene_symbol', 'N/A')})": t["id"] for t in targets}
    selected_target_key = st.selectbox("Select Target for Assay Data", list(target_options.keys()))
    
    if not selected_target_key:
        return
    
    target_id = target_options[selected_target_key]
    
    # Load assays and compounds
    col1, col2 = st.columns(2)
    
    with col1:
        st.markdown("#### Linked Assays")
        with st.spinner("Loading assays..."):
            assays = api_get(f"/api/v1/targets/{target_id}/assays")
        
        if assays:
            assay_df = pd.DataFrame(assays)
            st.dataframe(assay_df[["name", "experiment_type", "created_at"]], hide_index=True, use_container_width=True)
        else:
            st.info("No assays linked to this target yet.")
    
    with col2:
        st.markdown("#### Compound Activity")
        with st.spinner("Loading compounds..."):
            compounds = api_get(f"/api/v1/targets/{target_id}/compounds")
        
        if compounds:
            # Create DataFrame with activity data
            compound_data = []
            for comp in compounds:
                compound_data.append({
                    "Compound": comp.get("compound_name", comp["compound_id"]),
                    "SMILES": comp["smiles"][:50] + "..." if len(comp["smiles"]) > 50 else comp["smiles"],
                    "Activity": f"{comp.get('activity_value', 'N/A')} {comp.get('activity_units', '')}".strip(),
                    "Type": comp.get("activity_type", "N/A")
                })
            
            comp_df = pd.DataFrame(compound_data)
            st.dataframe(comp_df, hide_index=True, use_container_width=True)
            
            # Summary metrics
            st.markdown("#### Activity Summary")
            col1, col2, col3 = st.columns(3)
            
            with col1:
                st.metric("Total Compounds", len(compounds))
            
            with col2:
                activity_types = [c.get("activity_type") for c in compounds if c.get("activity_type")]
                unique_types = len(set(activity_types)) if activity_types else 0
                st.metric("Activity Types", unique_types)
            
            with col3:
                # Calculate potent compounds (assuming IC50 < 100 nM is potent)
                potent_count = 0
                for comp in compounds:
                    if (comp.get("activity_type") in ["IC50", "Ki", "Kd"] and 
                        comp.get("activity_value") and 
                        comp.get("activity_units") == "nM" and
                        comp["activity_value"] < 100):
                        potent_count += 1
                st.metric("Potent Compounds", potent_count)
        else:
            st.info("No compounds with activity data for this target.")
    
    # Future: Add "Link New Assay" button and form


def main():
    """Main Target Management page."""
    st.title("🎯 Target Management")
    st.markdown("Manage biological targets for drug discovery programs.")
    
    # Create tabs
    tab1, tab2, tab3, tab4 = st.tabs([
        "Target List", 
        "Target Details", 
        "Competitive Landscape", 
        "Assay Linkage"
    ])
    
    with tab1:
        render_target_list_tab()
    
    with tab2:
        render_target_detail_tab()
    
    with tab3:
        render_competitive_landscape_tab()
    
    with tab4:
        render_assay_linkage_tab()


if __name__ == "__main__":
    main()
