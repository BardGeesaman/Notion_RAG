"""Generative Chemistry Dashboard - De Novo Molecular Design."""

from __future__ import annotations

import io
import json
import pandas as pd
import requests
import streamlit as st
from typing import Dict, List, Optional

from scripts.dashboard.auth import require_auth


def render_generative_chemistry_page() -> None:
    """Render the Generative Chemistry page."""
    require_auth()
    
    st.header("🧪 Generative Chemistry")
    st.caption("AI-powered de novo molecular design and optimization.")
    
    # Initialize session state
    if "generative_results" not in st.session_state:
        st.session_state.generative_results = []
    
    # Create tabs
    tab1, tab2, tab3, tab4 = st.tabs(["🎲 Generate", "🔄 Interpolate", "⚡ Optimize", "📊 Results"])
    
    with tab1:
        render_generate_tab()
    
    with tab2:
        render_interpolate_tab()
    
    with tab3:
        render_optimize_tab()
    
    with tab4:
        render_results_tab()


def render_generate_tab() -> None:
    """Render the molecular generation tab."""
    st.subheader("Random Molecular Generation")
    st.caption("Generate novel molecules by sampling from the VAE latent space.")
    
    col1, col2 = st.columns([1, 2])
    
    with col1:
        st.markdown("### Parameters")
        
        n_samples = st.slider(
            "Number of molecules",
            min_value=1,
            max_value=100,
            value=10,
            help="How many molecules to generate"
        )
        
        temperature = st.slider(
            "Temperature",
            min_value=0.1,
            max_value=2.0,
            value=1.0,
            step=0.1,
            help="Higher temperature = more diverse structures"
        )
        
        max_length = st.slider(
            "Max SMILES length",
            min_value=20,
            max_value=200,
            value=100,
            help="Maximum length of generated SMILES"
        )
        
        generate_button = st.button("🎲 Generate Molecules", type="primary")
    
    with col2:
        st.markdown("### Results")
        
        if generate_button:
            with st.spinner("Generating molecules..."):
                try:
                    # Call API
                    response = requests.post(
                        f"{get_api_base_url()}/api/v1/generative/sample",
                        json={
                            "n_samples": n_samples,
                            "temperature": temperature,
                            "max_length": max_length,
                        },
                        headers=get_auth_headers(),
                        timeout=30
                    )
                    
                    if response.status_code == 200:
                        data = response.json()
                        molecules = data["molecules"]
                        
                        # Add to session state
                        for mol in molecules:
                            mol["source"] = "Generate"
                            mol["timestamp"] = pd.Timestamp.now().strftime("%Y-%m-%d %H:%M:%S")
                        
                        st.session_state.generative_results.extend(molecules)
                        
                        # Display results
                        st.success(f"Generated {len(molecules)} molecules!")
                        display_molecules_table(molecules)
                        
                    elif response.status_code == 503:
                        st.error("Generative model not available. Please contact administrator.")
                    else:
                        st.error(f"Generation failed: {response.text}")
                        
                except Exception as e:
                    st.error(f"Error generating molecules: {str(e)}")


def render_interpolate_tab() -> None:
    """Render the molecular interpolation tab."""
    st.subheader("Molecular Interpolation")
    st.caption("Explore chemical space between two molecules.")
    
    col1, col2 = st.columns([1, 2])
    
    with col1:
        st.markdown("### Input Molecules")
        
        smiles_start = st.text_input(
            "Start SMILES",
            placeholder="e.g., CCO",
            help="Starting molecule SMILES"
        )
        
        smiles_end = st.text_input(
            "End SMILES", 
            placeholder="e.g., CCC",
            help="Ending molecule SMILES"
        )
        
        steps = st.slider(
            "Interpolation steps",
            min_value=2,
            max_value=50,
            value=10,
            help="Number of intermediate molecules"
        )
        
        interpolation_type = st.selectbox(
            "Interpolation method",
            ["linear", "spherical"],
            help="Linear or spherical interpolation in latent space"
        )
        
        interpolate_button = st.button("🔄 Interpolate", type="primary")
    
    with col2:
        st.markdown("### Trajectory")
        
        if interpolate_button:
            if not smiles_start or not smiles_end:
                st.error("Please provide both start and end SMILES")
            else:
                with st.spinner("Interpolating molecules..."):
                    try:
                        response = requests.post(
                            f"{get_api_base_url()}/api/v1/generative/interpolate",
                            json={
                                "smiles_start": smiles_start,
                                "smiles_end": smiles_end,
                                "steps": steps,
                                "interpolation_type": interpolation_type,
                            },
                            headers=get_auth_headers(),
                            timeout=30
                        )
                        
                        if response.status_code == 200:
                            data = response.json()
                            molecules = data["molecules"]
                            
                            # Add metadata
                            for mol in molecules:
                                mol["source"] = f"Interpolate ({smiles_start} → {smiles_end})"
                                mol["timestamp"] = pd.Timestamp.now().strftime("%Y-%m-%d %H:%M:%S")
                            
                            st.session_state.generative_results.extend(molecules)
                            
                            st.success(f"Generated {len(molecules)} interpolation steps!")
                            display_interpolation_trajectory(molecules)
                            
                        elif response.status_code == 503:
                            st.error("Generative model not available.")
                        else:
                            st.error(f"Interpolation failed: {response.text}")
                            
                    except Exception as e:
                        st.error(f"Error during interpolation: {str(e)}")


def render_optimize_tab() -> None:
    """Render the property optimization tab."""
    st.subheader("Property-Guided Optimization")
    st.caption("Optimize molecules for desired properties.")
    
    col1, col2 = st.columns([1, 2])
    
    with col1:
        st.markdown("### Optimization Setup")
        
        seed_smiles = st.text_input(
            "Seed molecule SMILES",
            placeholder="e.g., CCO",
            help="Starting molecule for optimization"
        )
        
        st.markdown("#### Property Constraints")
        
        # Initialize constraints in session state
        if "opt_constraints" not in st.session_state:
            st.session_state.opt_constraints = []
        
        # Add constraint form
        with st.expander("Add New Constraint", expanded=len(st.session_state.opt_constraints) == 0):
            prop_name = st.selectbox(
                "Property",
                ["logp", "mw", "tpsa", "hbd", "hba", "herg", "logs"],
                help="Property to constrain"
            )
            
            constraint_type = st.radio(
                "Constraint type",
                ["Range", "Target"],
                horizontal=True
            )
            
            if constraint_type == "Range":
                col_min, col_max = st.columns(2)
                with col_min:
                    min_val = st.number_input("Min value", value=0.0, step=0.1)
                with col_max:
                    max_val = st.number_input("Max value", value=5.0, step=0.1)
                target_val = None
            else:
                target_val = st.number_input("Target value", value=2.5, step=0.1)
                min_val = max_val = None
            
            weight = st.slider("Weight", min_value=0.1, max_value=5.0, value=1.0, step=0.1)
            
            if st.button("Add Constraint"):
                constraint = {
                    "name": prop_name,
                    "min_value": min_val,
                    "max_value": max_val,
                    "target_value": target_val,
                    "weight": weight
                }
                st.session_state.opt_constraints.append(constraint)
                st.rerun()
        
        # Display current constraints
        if st.session_state.opt_constraints:
            st.markdown("#### Current Constraints")
            for i, constraint in enumerate(st.session_state.opt_constraints):
                col_info, col_remove = st.columns([3, 1])
                with col_info:
                    if constraint["target_value"] is not None:
                        st.text(f"{constraint['name']}: target={constraint['target_value']:.2f}")
                    else:
                        st.text(f"{constraint['name']}: {constraint['min_value']:.2f}-{constraint['max_value']:.2f}")
                with col_remove:
                    if st.button("🗑️", key=f"remove_{i}"):
                        st.session_state.opt_constraints.pop(i)
                        st.rerun()
        
        # Optimization parameters
        st.markdown("#### Optimization Parameters")
        n_iterations = st.slider("Iterations", min_value=1, max_value=1000, value=100)
        n_samples_per_iter = st.slider("Samples per iteration", min_value=1, max_value=50, value=10)
        learning_rate = st.slider("Learning rate", min_value=0.01, max_value=1.0, value=0.1, step=0.01)
        
        optimize_button = st.button("⚡ Optimize", type="primary")
    
    with col2:
        st.markdown("### Optimization Results")
        
        if optimize_button:
            if not seed_smiles:
                st.error("Please provide a seed SMILES")
            elif not st.session_state.opt_constraints:
                st.error("Please add at least one constraint")
            else:
                with st.spinner("Optimizing molecules..."):
                    try:
                        response = requests.post(
                            f"{get_api_base_url()}/api/v1/generative/optimize",
                            json={
                                "seed_smiles": seed_smiles,
                                "constraints": st.session_state.opt_constraints,
                                "n_iterations": n_iterations,
                                "n_samples_per_iter": n_samples_per_iter,
                                "learning_rate": learning_rate,
                            },
                            headers=get_auth_headers(),
                            timeout=120  # Longer timeout for optimization
                        )
                        
                        if response.status_code == 200:
                            data = response.json()
                            molecules = data["optimized"]
                            seed_properties = data["seed_properties"]
                            best_score = data["best_score"]
                            
                            # Add metadata
                            for mol in molecules:
                                mol["source"] = f"Optimize ({seed_smiles})"
                                mol["timestamp"] = pd.Timestamp.now().strftime("%Y-%m-%d %H:%M:%S")
                            
                            st.session_state.generative_results.extend(molecules)
                            
                            st.success(f"Optimization complete! Best score: {best_score:.3f}")
                            
                            # Show before/after comparison
                            display_optimization_results(seed_smiles, seed_properties, molecules)
                            
                        elif response.status_code == 503:
                            st.error("Generative model not available.")
                        else:
                            st.error(f"Optimization failed: {response.text}")
                            
                    except Exception as e:
                        st.error(f"Error during optimization: {str(e)}")


def render_results_tab() -> None:
    """Render the results management tab."""
    st.subheader("Session Results")
    st.caption("View and manage generated molecules from this session.")
    
    if not st.session_state.generative_results:
        st.info("No molecules generated yet. Use the other tabs to generate molecules.")
        return
    
    # Results summary
    col1, col2, col3 = st.columns(3)
    with col1:
        st.metric("Total Molecules", len(st.session_state.generative_results))
    with col2:
        sources = [mol.get("source", "Unknown") for mol in st.session_state.generative_results]
        unique_sources = len(set(sources))
        st.metric("Generation Methods", unique_sources)
    with col3:
        valid_smiles = sum(1 for mol in st.session_state.generative_results if mol.get("smiles"))
        st.metric("Valid SMILES", valid_smiles)
    
    # Results table
    st.markdown("### All Generated Molecules")
    
    # Convert to DataFrame
    df = pd.DataFrame(st.session_state.generative_results)
    
    # Add selection column
    df.insert(0, "Select", False)
    
    # Display editable dataframe
    edited_df = st.data_editor(
        df,
        use_container_width=True,
        hide_index=True,
        column_config={
            "Select": st.column_config.CheckboxColumn("Select", default=False),
            "smiles": st.column_config.TextColumn("SMILES", width="medium"),
            "score": st.column_config.NumberColumn("Score", format="%.3f"),
            "properties": st.column_config.TextColumn("Properties", width="large"),
            "source": st.column_config.TextColumn("Source", width="medium"),
            "timestamp": st.column_config.TextColumn("Generated", width="medium"),
        }
    )
    
    # Action buttons
    col1, col2, col3 = st.columns(3)
    
    with col1:
        if st.button("🧪 Register Selected", type="primary"):
            selected_rows = edited_df[edited_df["Select"] == True]
            if len(selected_rows) > 0:
                st.success(f"Would register {len(selected_rows)} molecules (feature not implemented)")
            else:
                st.warning("No molecules selected")
    
    with col2:
        if st.button("📁 Export CSV"):
            csv_buffer = io.StringIO()
            df.drop("Select", axis=1).to_csv(csv_buffer, index=False)
            st.download_button(
                label="💾 Download CSV",
                data=csv_buffer.getvalue(),
                file_name=f"generative_chemistry_results_{pd.Timestamp.now().strftime('%Y%m%d_%H%M%S')}.csv",
                mime="text/csv"
            )
    
    with col3:
        if st.button("🗑️ Clear All"):
            st.session_state.generative_results = []
            st.rerun()


def display_molecules_table(molecules: List[Dict]) -> None:
    """Display molecules in a formatted table."""
    if not molecules:
        return
    
    df = pd.DataFrame(molecules)
    
    # Format properties column
    if "properties" in df.columns:
        df["properties_str"] = df["properties"].apply(lambda x: json.dumps(x, indent=1) if x else "{}")
    
    st.dataframe(
        df[["smiles", "score", "properties_str"]].rename(columns={"properties_str": "properties"}),
        use_container_width=True,
        hide_index=True
    )


def display_interpolation_trajectory(molecules: List[Dict]) -> None:
    """Display interpolation trajectory."""
    if not molecules:
        return
    
    st.markdown("#### Interpolation Steps")
    
    # Create trajectory visualization
    df = pd.DataFrame(molecules)
    df["step_label"] = df.apply(lambda x: f"Step {x.get('step', 0)}: {x['smiles']}", axis=1)
    
    st.dataframe(
        df[["step", "smiles"]],
        use_container_width=True,
        hide_index=True
    )


def display_optimization_results(seed_smiles: str, seed_properties: Dict, optimized: List[Dict]) -> None:
    """Display optimization before/after comparison."""
    if not optimized:
        return
    
    st.markdown("#### Before vs After")
    
    col1, col2 = st.columns(2)
    
    with col1:
        st.markdown("**Original Molecule**")
        st.code(seed_smiles)
        if seed_properties:
            st.json(seed_properties)
    
    with col2:
        st.markdown("**Best Optimized Molecule**")
        best_mol = max(optimized, key=lambda x: x.get("score", 0))
        st.code(best_mol["smiles"])
        if best_mol.get("properties"):
            st.json(best_mol["properties"])
        st.metric("Optimization Score", f"{best_mol.get('score', 0):.3f}")


def get_api_base_url() -> str:
    """Get API base URL from environment or default."""
    import os
    return os.getenv("API_BASE_URL", "http://localhost:8000")


def get_auth_headers() -> Dict[str, str]:
    """Get authentication headers."""
    # In a real implementation, this would get the auth token
    return {"Content-Type": "application/json"}
