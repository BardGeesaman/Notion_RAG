"""Signature linking, lead optimization, and procurement tabs."""

from __future__ import annotations

import json

import pandas as pd
import streamlit as st


def render_lead_optimization_tab(
    tab_signature,
    tab_lead_opt,
    tab_procurement,
    *,
    db_session,
    compound_model,
    adme_model,
    pk_model,
    tox_model,
    link_compound_to_signature,
    search_vendors,
    get_vendor_info,
) -> None:
    """Render signature linking, lead optimization, and procurement tabs."""
    _render_signature_links(tab_signature, db_session, compound_model, link_compound_to_signature)
    _render_lead_optimization(
        tab_lead_opt,
        db_session=db_session,
        compound_model=compound_model,
        adme_model=adme_model,
        pk_model=pk_model,
        tox_model=tox_model,
    )
    _render_procurement_tab(tab_procurement, search_vendors=search_vendors, get_vendor_info=get_vendor_info)


def _render_signature_links(tab, db_session, compound_model, link_compound_to_signature):
    with tab:
        st.subheader("Signature Links")

        def _load_compound_ids() -> list[str]:
            try:
                with db_session() as db:
                    compounds = db.query(compound_model).order_by(compound_model.created_at.desc()).limit(500).all()
                    return [c.compound_id for c in compounds]
            except Exception:
                return []

        compound_options = _load_compound_ids()

        st.markdown("### Link Compound to Signature")
        col_l, col_r = st.columns(2)
        with col_l:
            compound_sel = st.selectbox("Compound", compound_options, key="siglink_compound_sel")
            signature_id = st.text_input("Signature ID (UUID)", key="siglink_signature_id")
            effect_type = st.selectbox(
                "Effect Type",
                ["reverses", "mimics", "partial", "unknown"],
                key="siglink_effect_type",
            )
        with col_r:
            correlation = st.number_input("Correlation", value=0.0, step=0.01, format="%.3f", key="siglink_corr")
            p_value = st.number_input("p-value", value=1.0, step=0.001, format="%.4f", key="siglink_pval")
            evidence_source = st.text_input("Evidence Source", key="siglink_evidence")

        if st.button("Create Link", key="siglink_create_btn"):
            if compound_sel and signature_id:
                try:
                    link_compound_to_signature(
                        compound_id=compound_sel,
                        signature_id=signature_id,
                        effect_type=effect_type,
                        correlation=correlation,
                        p_value=p_value,
                        evidence_source=evidence_source or None,
                    )
                    st.success(f"Linked compound {compound_sel} to signature {signature_id} ({effect_type}).")
                except Exception as e:
                    st.error(f"Error creating link: {e}")
            else:
                st.warning("Compound and Signature ID are required.")

        st.markdown("---")
        st.info("Use linking above; advanced signature searches removed.")


def _render_lead_optimization(
    tab,
    *,
    db_session,
    compound_model,
    adme_model,
    pk_model,
    tox_model,
):
    with tab:
        st.subheader("Lead Optimization")
        try:
            with db_session() as db:
                compound_choices = [
                    c.compound_id for c in db.query(compound_model).order_by(compound_model.created_at.desc()).limit(500).all()
                ]
        except Exception:
            compound_choices = []

        col_lo1, col_lo2 = st.columns(2)

        with col_lo1:
            st.markdown("### ADME Data")
            adme_comp = st.selectbox("Compound (ADME)", compound_choices, key="lo_adme_comp")
            adme_assay = st.selectbox("Assay Type", ["permeability", "stability", "cyp_inhibition", "other"], key="lo_adme_assay")
            adme_value = st.number_input("Value", value=0.0, step=0.1, key="lo_adme_value")
            adme_unit = st.text_input("Unit", value="uM", key="lo_adme_unit")
            adme_conditions = st.text_area("Conditions (JSON)", key="lo_adme_conditions")
            if st.button("Save ADME", key="lo_adme_save"):
                try:
                    with db_session() as db:
                        cond = None
                        if adme_conditions.strip():
                            cond = json.loads(adme_conditions)
                        comp_obj = (
                            db.query(compound_model).filter(compound_model.compound_id == adme_comp).first()
                            if adme_comp
                            else None
                        )
                        rec = adme_model(
                            compound_id=comp_obj.id if comp_obj else None,
                            assay_type=adme_assay,
                            value=adme_value,
                            unit=adme_unit or None,
                            conditions=cond,
                        )
                        db.add(rec)
                        db.commit()
                        st.success("ADME entry saved.")
                except Exception as e:
                    st.error(f"Failed to save ADME: {e}")

        with col_lo2:
            st.markdown("### PK Studies")
            pk_comp = st.selectbox("Compound (PK)", compound_choices, key="lo_pk_comp")
            pk_species = st.text_input("Species", value="mouse", key="lo_pk_species")
            pk_route = st.text_input("Route", value="IV", key="lo_pk_route")
            pk_dose = st.number_input("Dose", value=0.0, step=0.1, key="lo_pk_dose")
            pk_dose_unit = st.text_input("Dose Unit", value="mg/kg", key="lo_pk_dose_unit")
            pk_auc = st.number_input("AUC", value=0.0, step=0.1, key="lo_pk_auc")
            pk_cmax = st.number_input("Cmax", value=0.0, step=0.1, key="lo_pk_cmax")
            pk_tmax = st.number_input("Tmax", value=0.0, step=0.1, key="lo_pk_tmax")
            pk_half = st.number_input("Half-life", value=0.0, step=0.1, key="lo_pk_half")
            pk_f = st.number_input("Bioavailability (%)", value=0.0, step=0.1, key="lo_pk_f")
            pk_cl = st.number_input("Clearance", value=0.0, step=0.1, key="lo_pk_cl")
            pk_vd = st.number_input("Volume of Distribution", value=0.0, step=0.1, key="lo_pk_vd")
            if st.button("Save PK", key="lo_pk_save"):
                try:
                    with db_session() as db:
                        comp_obj = (
                            db.query(compound_model).filter(compound_model.compound_id == pk_comp).first()
                            if pk_comp
                            else None
                        )
                        rec = pk_model(
                            compound_id=comp_obj.id if comp_obj else None,
                            species=pk_species or None,
                            route=pk_route or None,
                            dose=pk_dose or None,
                            dose_unit=pk_dose_unit or None,
                            auc=pk_auc or None,
                            c_max=pk_cmax or None,
                            t_max=pk_tmax or None,
                            half_life=pk_half or None,
                            bioavailability=pk_f or None,
                            clearance=pk_cl or None,
                            vd=pk_vd or None,
                        )
                        db.add(rec)
                        db.commit()
                        st.success("PK entry saved.")
                except Exception as e:
                    st.error(f"Failed to save PK: {e}")

        st.markdown("---")
        col_tox, col_profile = st.columns(2)

        with col_tox:
            st.markdown("### Toxicology")
            tox_comp = st.selectbox("Compound (Tox)", compound_choices, key="lo_tox_comp")
            tox_assay = st.text_input("Assay Type", value="herg", key="lo_tox_assay")
            tox_value = st.number_input("Value", value=0.0, step=0.1, key="lo_tox_value")
            tox_unit = st.text_input("Unit", value="uM", key="lo_tox_unit")
            tox_result = st.selectbox("Result", ["negative", "positive", "inconclusive"], key="lo_tox_result")
            tox_conditions = st.text_area("Conditions (JSON)", key="lo_tox_conditions")
            if st.button("Save Toxicology", key="lo_tox_save"):
                try:
                    with db_session() as db:
                        cond = None
                        if tox_conditions.strip():
                            cond = json.loads(tox_conditions)
                        comp_obj = (
                            db.query(compound_model).filter(compound_model.compound_id == tox_comp).first()
                            if tox_comp
                            else None
                        )
                        rec = tox_model(
                            compound_id=comp_obj.id if comp_obj else None,
                            assay_type=tox_assay or None,
                            value=tox_value or None,
                            unit=tox_unit or None,
                            result=tox_result or None,
                            conditions=cond,
                        )
                        db.add(rec)
                        db.commit()
                        st.success("Toxicology entry saved.")
                except Exception as e:
                    st.error(f"Failed to save Toxicology: {e}")

        with col_profile:
            st.markdown("### Compound Profile")
            profile_comp = st.selectbox("Compound (Profile)", compound_choices, key="lo_profile_comp")
            if profile_comp:
                with db_session() as db:
                    comp = db.query(compound_model).filter(compound_model.compound_id == profile_comp).first()
                    if comp:
                        adme = db.query(adme_model).filter(adme_model.compound_id == comp.id).all()
                        pk = db.query(pk_model).filter(pk_model.compound_id == comp.id).all()
                        tox = db.query(tox_model).filter(tox_model.compound_id == comp.id).all()

                        if adme:
                            st.markdown("**ADME Results**")
                            st.dataframe(
                                [
                                    {
                                        "assay_type": a.assay_type,
                                        "value": a.value,
                                        "unit": a.unit,
                                        "conditions": a.conditions,
                                    }
                                    for a in adme
                                ],
                                hide_index=True,
                                use_container_width=True,
                            )
                        if pk:
                            st.markdown("**PK Studies**")
                            st.dataframe(
                                [
                                    {
                                        "species": p.species,
                                        "route": p.route,
                                        "auc": p.auc,
                                        "c_max": p.c_max,
                                        "t_max": p.t_max,
                                        "half_life": p.half_life,
                                        "bioavailability": p.bioavailability,
                                    }
                                    for p in pk
                                ],
                                hide_index=True,
                                use_container_width=True,
                            )
                        if tox:
                            st.markdown("**Toxicology Results**")
                            st.dataframe(
                                [
                                    {
                                        "assay_type": t.assay_type,
                                        "value": t.value,
                                        "unit": t.unit,
                                        "result": t.result,
                                    }
                                    for t in tox
                                ],
                                hide_index=True,
                                use_container_width=True,
                            )

                        warnings_list = []
                        for t in tox:
                            if t.assay_type and t.assay_type.lower() == "herg" and t.value is not None and t.value < 10:
                                warnings_list.append("Potential hERG liability (IC50 < 10 uM)")
                            if t.assay_type and t.assay_type.lower() == "ames" and t.result and t.result.lower() == "positive":
                                warnings_list.append("Ames positive")
                        for p in pk:
                            if p.bioavailability is not None and p.bioavailability < 10:
                                warnings_list.append("Low bioavailability (<10%)")

                        if warnings_list:
                            for w in warnings_list:
                                st.error(w)
                        else:
                            st.success("No major liabilities detected.")


def _render_procurement_tab(tab, *, search_vendors, get_vendor_info):
    with tab:
        st.subheader("🔍 Vendor Search")
        st.markdown("Search for compounds across multiple vendor catalogs.")

        st.markdown("### Supported Vendors")
        vendors = get_vendor_info()
        vendor_cols = st.columns(min(len(vendors), 4)) if vendors else []
        for idx, vendor in enumerate(vendors):
            col_idx = idx % len(vendor_cols)
            with vendor_cols[col_idx]:
                st.markdown(f"**{vendor['name']}**")
                st.caption(f"[Visit Site]({vendor['url']})")

        st.markdown("---")

        search_query = st.text_input(
            "Search Compound",
            placeholder="Enter SMILES, CAS number, or compound name",
            key="procurement_search",
        )

        if st.button("🔍 Search Vendors", type="primary"):
            if not search_query or not search_query.strip():
                st.warning("Please enter a search query.")
            else:
                with st.spinner("Searching vendor catalogs..."):
                    try:
                        results = search_vendors(search_query.strip())

                        if results:
                            st.success(f"Found {len(results)} vendor results")

                            results_data = []
                            for result in results:
                                results_data.append(
                                    {
                                        "Vendor": result.get("vendor", "Unknown"),
                                        "Catalog ID": result.get("catalog_id", "N/A"),
                                        "Price": result.get("price", "N/A"),
                                        "Availability": result.get("availability", "Unknown"),
                                        "Link": result.get("url", ""),
                                    }
                                )

                            df_results = pd.DataFrame(results_data)
                            st.dataframe(df_results, use_container_width=True, hide_index=True)

                            st.markdown("---")
                            st.subheader("Vendor Links")
                            for result in results:
                                vendor_name = result.get("vendor", "Unknown")
                                url = result.get("url", "")
                                catalog_id = result.get("catalog_id", "N/A")
                                if url:
                                    st.markdown(f"- **{vendor_name}** ({catalog_id}): [View Product]({url})")
                        else:
                            st.info("No results found. Try a different search query.")

                    except Exception as e:
                        st.error(f"Error searching vendors: {e}")
                        st.exception(e)

