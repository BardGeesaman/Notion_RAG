"""Compare page for comparing entities side-by-side."""
from __future__ import annotations

import streamlit as st

from amprenta_rag.utils.comparison import compare_experiments, compare_compounds
from scripts.dashboard.db_session import db_session


def render_compare_page() -> None:
    """Render the Compare page."""
    st.header("🔍 Compare")
    st.markdown("Compare two experiments or compounds side-by-side.")
    
    # Entity type selection
    entity_type = st.radio(
        "Compare",
        ["Experiments", "Compounds"],
        horizontal=True,
        key="compare_entity_type"
    )
    
    with db_session() as db:
        if entity_type == "Experiments":
            from amprenta_rag.database.models import Experiment
            
            experiments = db.query(Experiment).order_by(Experiment.name).all()
            exp_options = {exp.name: str(exp.id) for exp in experiments}
            
            if not exp_options:
                st.info("No experiments available to compare.")
                return
            
            col1, col2 = st.columns(2)
            with col1:
                item_a_name = st.selectbox("Item A", list(exp_options.keys()), key="compare_exp_a")
            with col2:
                item_b_name = st.selectbox("Item B", list(exp_options.keys()), key="compare_exp_b")
            
            if st.button("Compare", type="primary"):
                if item_a_name == item_b_name:
                    st.error("Please select two different experiments.")
                else:
                    try:
                        result = compare_experiments(exp_options[item_a_name], exp_options[item_b_name], db)
                        
                        st.markdown("---")
                        st.subheader("Comparison Results")
                        
                        col1, col2 = st.columns(2)
                        
                        with col1:
                            st.markdown("### Item A")
                            for key, value in result["item1"].items():
                                is_different = key in result["differences"]
                                prefix = "🔴 " if is_different else ""
                                st.markdown(f"**{key}:** {prefix}{value or 'N/A'}")
                        
                        with col2:
                            st.markdown("### Item B")
                            for key, value in result["item2"].items():
                                is_different = key in result["differences"]
                                prefix = "🔴 " if is_different else ""
                                st.markdown(f"**{key}:** {prefix}{value or 'N/A'}")
                        
                        if result["differences"]:
                            st.info(f"Found {len(result['differences'])} difference(s): {', '.join(result['differences'])}")
                        else:
                            st.success("No differences found - items are identical.")
                    
                    except Exception as e:
                        st.error(f"Comparison failed: {e}")
        
        else:  # Compounds
            from amprenta_rag.database.models import Compound
            
            compounds = db.query(Compound).order_by(Compound.compound_id).all()
            comp_options = {comp.compound_id: str(comp.id) for comp in compounds}
            
            if not comp_options:
                st.info("No compounds available to compare.")
                return
            
            col1, col2 = st.columns(2)
            with col1:
                item_a_id = st.selectbox("Item A", list(comp_options.keys()), key="compare_comp_a")
            with col2:
                item_b_id = st.selectbox("Item B", list(comp_options.keys()), key="compare_comp_b")
            
            if st.button("Compare", type="primary"):
                if item_a_id == item_b_id:
                    st.error("Please select two different compounds.")
                else:
                    try:
                        result = compare_compounds(comp_options[item_a_id], comp_options[item_b_id], db)
                        
                        st.markdown("---")
                        st.subheader("Comparison Results")
                        
                        col1, col2 = st.columns(2)
                        
                        with col1:
                            st.markdown("### Item A")
                            for key, value in result["item1"].items():
                                is_different = key in result["differences"]
                                prefix = "🔴 " if is_different else ""
                                if key == "smiles" and value:
                                    # Truncate long SMILES
                                    display_value = value[:50] + "..." if len(value) > 50 else value
                                    st.markdown(f"**{key}:** {prefix}`{display_value}`")
                                else:
                                    st.markdown(f"**{key}:** {prefix}{value or 'N/A'}")
                        
                        with col2:
                            st.markdown("### Item B")
                            for key, value in result["item2"].items():
                                is_different = key in result["differences"]
                                prefix = "🔴 " if is_different else ""
                                if key == "smiles" and value:
                                    # Truncate long SMILES
                                    display_value = value[:50] + "..." if len(value) > 50 else value
                                    st.markdown(f"**{key}:** {prefix}`{display_value}`")
                                else:
                                    st.markdown(f"**{key}:** {prefix}{value or 'N/A'}")
                        
                        if result["differences"]:
                            st.info(f"Found {len(result['differences'])} difference(s): {', '.join(result['differences'])}")
                        else:
                            st.success("No differences found - items are identical.")
                    
                    except Exception as e:
                        st.error(f"Comparison failed: {e}")


if __name__ == "__main__":
    render_compare_page()
