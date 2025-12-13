"""Compounds, registration, and structure search tabs."""

from __future__ import annotations

import pandas as pd
import streamlit as st


def render_compounds_tab(
    tab_compounds,
    tab_register,
    tab_structure,
    *,
    db_session,
    compound_model,
    register_compound,
    export_compounds,
    substructure_search,
    similarity_search,
    pharmacophore_search,
    get_pharmacophore_features,
) -> None:
    """Render Compounds, Register, and Structure Search tabs."""
    _render_compounds(tab_compounds, db_session, compound_model, export_compounds)
    _render_register(tab_register, register_compound)
    _render_structure_search(
        tab_structure,
        db_session=db_session,
        substructure_search=substructure_search,
        similarity_search=similarity_search,
        pharmacophore_search=pharmacophore_search,
        get_pharmacophore_features=get_pharmacophore_features,
        compound_model=compound_model,
    )


def _render_compounds(tab, db_session, compound_model, export_compounds):
    with tab:
        st.subheader("Compounds")
        search_term = st.text_input("Search compounds by ID or SMILES", "", key="compound_search")

        with db_session() as db:
            query = db.query(compound_model)
            if search_term:
                query = query.filter(
                    (compound_model.compound_id.ilike(f"%{search_term}%"))
                    | (compound_model.smiles.ilike(f"%{search_term}%"))
                )

            compounds = query.limit(100).all()

            st.metric("Compounds Found", len(compounds))
            if len(compounds) == 100:
                st.info("Showing first 100 compounds. Use search to narrow down.")

            if compounds:
                compound_data = []
                for compound in compounds:
                    hts_count = len(compound.hts_results) if compound.hts_results else 0
                    bio_count = len(compound.biochemical_results) if compound.biochemical_results else 0
                    compound_data.append(
                        {
                            "Compound ID": compound.compound_id,
                            "SMILES": compound.smiles[:50] + "..." if len(compound.smiles) > 50 else compound.smiles,
                            "InChI Key": compound.inchi_key or "-",
                            "Molecular Weight": compound.molecular_weight or "-",
                            "HTS Results": hts_count,
                            "Biochemical Results": bio_count,
                            "Created": compound.created_at.strftime("%Y-%m-%d"),
                        }
                    )
                df_compounds = pd.DataFrame(compound_data)
                st.dataframe(df_compounds, width="stretch", hide_index=True)

                st.markdown("### Export")
                col1, col2 = st.columns([1, 3])
                with col1:
                    export_format = st.selectbox("Format", ["csv", "json", "excel"], key="compound_export_format")
                with col2:
                    compound_ids = [c.compound_id for c in compounds]
                    mime_types = {
                        "csv": "text/csv",
                        "json": "application/json",
                        "excel": "application/vnd.openxmlformats-officedocument.spreadsheetml.sheet",
                    }
                    file_extensions = {"csv": "csv", "json": "json", "excel": "xlsx"}

                    try:
                        export_data = export_compounds(compound_ids, export_format, db)
                        st.download_button(
                            label=f"📥 Download Compounds ({export_format.upper()})",
                            data=export_data,
                            file_name=f"compounds.{file_extensions[export_format]}",
                            mime=mime_types[export_format],
                            key="export_compounds",
                        )
                    except Exception as e:
                        st.error(f"Export failed: {e}")

                st.markdown("---")
                st.subheader("Compound Details")

                for compound in compounds[:10]:
                    with st.expander(f"**{compound.compound_id}**"):
                        col1, col2 = st.columns(2)
                        with col1:
                            st.write(f"**ID:** `{compound.id}`")
                            st.write(f"**Compound ID:** {compound.compound_id}")
                            st.write(f"**SMILES:** `{compound.smiles}`")
                            if compound.inchi_key:
                                st.write(f"**InChI Key:** `{compound.inchi_key}`")
                            if compound.canonical_smiles:
                                st.write(f"**Canonical SMILES:** `{compound.canonical_smiles}`")
                        with col2:
                            st.write(f"**Created:** {compound.created_at.strftime('%Y-%m-%d %H:%M')}")
                            if compound.molecular_weight:
                                st.write(f"**Molecular Weight:** {compound.molecular_weight}")
                            if compound.molecular_formula:
                                st.write(f"**Molecular Formula:** {compound.molecular_formula}")
                            if compound.logp is not None:
                                st.write(f"**LogP:** {compound.logp}")
                            if compound.hbd_count is not None:
                                st.write(f"**HBD Count:** {compound.hbd_count}")
                            if compound.hba_count is not None:
                                st.write(f"**HBA Count:** {compound.hba_count}")
            else:
                st.info("No compounds found.")


def _render_register(tab, register_compound):
    with tab:
        st.subheader("Register Compound")
        name = st.text_input("Compound Name")
        smiles = st.text_input("SMILES")
        salt_form = st.selectbox(
            "Salt Form",
            ["None", "HCl", "sodium", "potassium", "free base", "other"],
            index=0,
        )
        if st.button("Register", type="primary"):
            if not name or not smiles:
                st.error("Name and SMILES are required.")
            else:
                try:
                    corp_id = register_compound(
                        name=name,
                        smiles=smiles,
                        salt_form=None if salt_form == "None" else salt_form,
                    )
                    st.success(f"Registered compound ID: {corp_id}")
                except Exception as e:
                    st.error(f"Registration failed: {e}")


def _render_structure_search(
    tab,
    *,
    db_session,
    substructure_search,
    similarity_search,
    pharmacophore_search,
    get_pharmacophore_features,
    compound_model,
):
    with tab:
        st.subheader("Structure Search")
        mode = st.radio("Search Mode", ["Substructure", "Similarity"], horizontal=True)
        query = st.text_input("SMARTS Pattern" if mode == "Substructure" else "Query SMILES")
        threshold = st.slider("Similarity Threshold", min_value=0.5, max_value=1.0, value=0.7, step=0.01)

        if st.button("Search", type="primary"):
            try:
                if mode == "Substructure":
                    results = substructure_search(query)
                else:
                    results = similarity_search(query, threshold=threshold)

                if not results:
                    st.info("No matches found.")
                else:
                    st.dataframe(results, hide_index=True, use_container_width=True)
            except Exception as e:
                st.error(f"Search failed: {e}")

        st.divider()
        st.markdown("### Pharmacophore Search")
        st.markdown("Search compounds by pharmacophore features (HBD, HBA, aromatic rings).")

        col1, col2, col3 = st.columns(3)
        with col1:
            st.markdown("**Hydrogen Bond Donors (HBD)**")
            min_hbd = st.number_input("Min HBD", min_value=0, value=0, key="min_hbd")
            max_hbd = st.number_input("Max HBD", min_value=0, value=10, key="max_hbd")

        with col2:
            st.markdown("**Hydrogen Bond Acceptors (HBA)**")
            min_hba = st.number_input("Min HBA", min_value=0, value=0, key="min_hba")
            max_hba = st.number_input("Max HBA", min_value=0, value=10, key="max_hba")

        with col3:
            st.markdown("**Aromatic Rings**")
            min_aromatic = st.number_input("Min Aromatic", min_value=0, value=0, key="min_aromatic")
            max_aromatic = st.number_input("Max Aromatic", min_value=0, value=5, key="max_aromatic")

        if st.button("Search by Pharmacophore", type="primary", key="pharmacophore_search"):
            try:
                with db_session() as db:
                    results = pharmacophore_search(
                        min_hbd=min_hbd if min_hbd > 0 else None,
                        max_hbd=max_hbd if max_hbd < 10 else None,
                        min_hba=min_hba if min_hba > 0 else None,
                        max_hba=max_hba if max_hba < 10 else None,
                        min_aromatic=min_aromatic if min_aromatic > 0 else None,
                        max_aromatic=max_aromatic if max_aromatic < 5 else None,
                        db=db,
                        limit=100,
                    )

                    if not results:
                        st.info("No compounds found matching the criteria.")
                    else:
                        results_data = []
                        for compound in results:
                            features = get_pharmacophore_features(compound.smiles) if compound.smiles else {}
                            results_data.append(
                                {
                                    "Compound ID": compound.compound_id,
                                    "SMILES": compound.smiles[:50] if compound.smiles else "",
                                    "HBD": features.get("hbd", compound.hbd_count or 0),
                                    "HBA": features.get("hba", compound.hba_count or 0),
                                    "Aromatic Rings": features.get("aromatic_rings", 0),
                                    "LogP": compound.logp,
                                    "MW": compound.molecular_weight,
                                }
                            )

                        st.dataframe(pd.DataFrame(results_data), hide_index=True, use_container_width=True)
                        st.info(f"Found {len(results)} compounds")
            except Exception as e:
                st.error(f"Pharmacophore search failed: {e}")
                st.exception(e)

