"""Chemical Structure Sketcher page."""

from __future__ import annotations

import streamlit as st

from scripts.dashboard.components.ketcher_sketcher import render_ketcher_editor


def render_chemical_sketcher_page() -> None:
    """Render the Chemical Structure Sketcher page."""
    from scripts.dashboard.auth import require_auth
    require_auth()
    
    st.header("✏️ Chemical Structure Sketcher")
    st.caption("Draw chemical structures with Ketcher and register compounds.")
    
    tab1, tab2 = st.tabs(["Draw Structure", "Structure Search"])
    
    with tab1:
        render_draw_structure_tab()
    
    with tab2:
        render_structure_search_tab()


def render_draw_structure_tab() -> None:
    """Render the structure drawing tab."""
    st.subheader("Draw Chemical Structure")
    
    # SMILES input for quick loading
    smiles_input = st.text_input(
        "Load SMILES (optional)",
        placeholder="e.g., CCO for ethanol",
        help="Enter SMILES to load into editor",
    )
    
    # Ketcher editor
    render_ketcher_editor(
        key="main_sketcher",
        initial_smiles=smiles_input if smiles_input else None,
        width=800,
        height=600,
    )
    
    st.markdown("---")
    st.markdown("### Instructions")
    st.markdown("""
    1. **Draw** your structure in the editor above
    2. Click **Get SMILES** button in the editor to export
    3. **Copy** the SMILES from the display
    4. **Paste** below to register or search
    """)
    
    # SMILES registration
    st.markdown("### Register Compound")
    
    manual_smiles = st.text_input(
        "SMILES to Register",
        placeholder="Paste SMILES from editor above",
        help="SMILES string exported from Ketcher",
        key="register_smiles",
    )
    
    compound_name = st.text_input(
        "Compound Name (optional)",
        placeholder="e.g., My Novel Compound",
    )
    
    if st.button("📝 Register Compound", type="primary", disabled=not manual_smiles):
        try:
            from amprenta_rag.chemistry.registration import register_compound
            
            with st.spinner("Registering compound..."):
                result = register_compound(
                    smiles=manual_smiles,
                    name=compound_name if compound_name else None,
                )
                
                if result:
                    st.success(f"Compound registered! ID: {result.get('compound_id', 'unknown')}")
                    
                    # Show basic info
                    if result.get('molecular_weight'):
                        st.metric("Molecular Weight", f"{result['molecular_weight']:.2f}")
                    if result.get('canonical_smiles'):
                        st.code(result['canonical_smiles'], language=None)
                else:
                    st.warning("Registration returned no result")
        except Exception as e:
            st.error(f"Registration failed: {e}")


def render_structure_search_tab() -> None:
    """Render the structure search tab."""
    st.subheader("Structure Search")
    
    search_smiles = st.text_input(
        "Query SMILES",
        placeholder="Paste SMILES from sketcher",
        help="SMILES to use as search query",
        key="search_smiles",
    )
    
    search_type = st.radio(
        "Search Type",
        options=["Substructure", "Similarity"],
        horizontal=True,
    )
    
    if search_type == "Similarity":
        threshold = st.slider(
            "Similarity Threshold",
            min_value=0.0,
            max_value=1.0,
            value=0.7,
            step=0.05,
            help="Tanimoto similarity threshold",
        )
    
    if st.button("🔍 Search", type="primary", disabled=not search_smiles):
        try:
            with st.spinner(f"Performing {search_type.lower()} search..."):
                if search_type == "Substructure":
                    from amprenta_rag.chemistry.search import substructure_search
                    
                    results = substructure_search(search_smiles, limit=50)
                else:
                    from amprenta_rag.chemistry.search import similarity_search
                    
                    results = similarity_search(search_smiles, threshold=threshold, limit=50)
                
                st.session_state["search_results"] = results
        except Exception as e:
            st.error(f"Search failed: {e}")
    
    # Display results
    results = st.session_state.get("search_results", [])
    
    if results:
        st.success(f"Found {len(results)} matching compounds")
        
        for i, compound in enumerate(results[:10], 1):  # Show first 10
            with st.expander(f"{i}. {compound.get('compound_id', 'Unknown')}"):
                st.code(compound.get('smiles', ''), language=None)
                
                if compound.get('name'):
                    st.markdown(f"**Name:** {compound['name']}")
                if compound.get('similarity'):
                    st.metric("Similarity", f"{compound['similarity']:.3f}")
        
        if len(results) > 10:
            st.caption(f"...and {len(results) - 10} more results")
    elif "search_results" in st.session_state:
        st.info("No matching compounds found.")


if __name__ == "__main__":
    render_chemical_sketcher_page()

