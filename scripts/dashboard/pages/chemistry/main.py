"""Chemistry page orchestration."""

from __future__ import annotations

import streamlit as st

import amprenta_rag.database.models as db_models
from amprenta_rag.database.session import db_session
from amprenta_rag.chemistry.compound_linking import link_compound_to_signature
from amprenta_rag.chemistry.registration import register_compound
from amprenta_rag.chemistry.structure_search import substructure_search, similarity_search
from amprenta_rag.chemistry.sar_analysis import (
    get_compound_properties,
    get_activity_data,
    calculate_lipinski,
    detect_activity_cliffs,
)
from amprenta_rag.chemistry.rgroup import find_common_core, decompose_rgroups, get_rgroup_statistics
from amprenta_rag.chemistry.pharmacophore import get_pharmacophore_features, pharmacophore_search
from amprenta_rag.chemistry.procurement import search_vendors, get_vendor_info
from amprenta_rag.utils.data_export import export_compounds

from .compounds import render_compounds_tab
from .sar_analysis import render_sar_tab
from .assay_results import render_assays_tab
from .lead_optimization import render_lead_optimization_tab


Compound = db_models.Compound
HTSCampaign = db_models.HTSCampaign
BiochemicalResult = db_models.BiochemicalResult
ADMEResult = db_models.ADMEResult
PKStudy = db_models.PKStudy
ToxicologyResult = db_models.ToxicologyResult


def render_chemistry_page() -> None:
    """Render the Chemistry page with modular tab handlers."""
    st.header("⚗️ Chemistry")

    tab_compounds, tab_reg, tab_struct, tab_sar, tab_hts, tab_bio, tab_sig, tab_lo, tab_proc = st.tabs(
        [
            "Compounds",
            "Register Compound",
            "Structure Search",
            "SAR Analysis",
            "HTS Campaigns",
            "Biochemical Results",
            "Signature Links",
            "Lead Optimization",
            "Procurement",
        ]
    )

    render_compounds_tab(
        tab_compounds,
        tab_reg,
        tab_struct,
        db_session=db_session,
        compound_model=Compound,
        register_compound=register_compound,
        export_compounds=export_compounds,
        substructure_search=substructure_search,
        similarity_search=similarity_search,
        pharmacophore_search=pharmacophore_search,
        get_pharmacophore_features=get_pharmacophore_features,
    )

    render_sar_tab(
        tab_sar,
        db_session=db_session,
        compound_model=Compound,
        get_compound_properties=get_compound_properties,
        calculate_lipinski=calculate_lipinski,
        get_activity_data=get_activity_data,
        detect_activity_cliffs=detect_activity_cliffs,
        find_common_core=find_common_core,
        decompose_rgroups=decompose_rgroups,
        get_rgroup_statistics=get_rgroup_statistics,
    )

    render_assays_tab(
        tab_hts,
        tab_bio,
        db_session=db_session,
        hts_campaign_model=HTSCampaign,
        biochemical_result_model=BiochemicalResult,
    )

    render_lead_optimization_tab(
        tab_sig,
        tab_lo,
        tab_proc,
        db_session=db_session,
        compound_model=Compound,
        adme_model=ADMEResult,
        pk_model=PKStudy,
        tox_model=ToxicologyResult,
        link_compound_to_signature=link_compound_to_signature,
        search_vendors=search_vendors,
        get_vendor_info=get_vendor_info,
    )

