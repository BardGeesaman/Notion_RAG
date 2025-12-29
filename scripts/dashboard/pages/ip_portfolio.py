"""IP Portfolio Dashboard page."""

from __future__ import annotations

import os
from typing import Any

import httpx
import pandas as pd
import streamlit as st


API_BASE = os.environ.get("API_URL", "http://localhost:8000")


def _api_get(path: str, *, timeout: int = 30) -> Any:
    """Make GET request to API."""
    with httpx.Client(timeout=timeout) as client:
        r = client.get(f"{API_BASE}{path}")
    r.raise_for_status()
    return r.json()


def _api_post(path: str, payload: dict, *, timeout: int = 30) -> Any:
    """Make POST request to API."""
    with httpx.Client(timeout=timeout) as client:
        r = client.post(f"{API_BASE}{path}", json=payload)
    r.raise_for_status()
    return r.json()


def render_ip_portfolio_page() -> None:
    """Render the IP Portfolio page."""
    from scripts.dashboard.auth import require_auth
    require_auth()
    
    st.header("📋 IP Portfolio")
    st.caption("Manage invention disclosures, patents, and evidence links.")
    
    tab1, tab2, tab3 = st.tabs(["Disclosures", "Patents", "Evidence Links"])
    
    with tab1:
        render_disclosures_tab()
    
    with tab2:
        render_patents_tab()
    
    with tab3:
        render_evidence_links_tab()


def render_disclosures_tab() -> None:
    """Render disclosures management tab."""
    st.subheader("Invention Disclosures")
    
    # New disclosure form
    with st.expander("➕ New Disclosure", expanded=False):
        title = st.text_input("Title", placeholder="Brief title of invention")
        description = st.text_area("Description", placeholder="Technical description")
        technical_field = st.text_input("Technical Field", placeholder="e.g., Chemistry, Biology")
        
        st.markdown("**Inventors**")
        user_id = st.text_input("Inventor User ID", placeholder="UUID")
        is_primary = st.checkbox("Primary Inventor", value=True)
        contribution = st.slider("Contribution %", 0, 100, 100)
        
        creator_id = st.text_input("Your User ID", placeholder="Creating user UUID")
        company_id = st.text_input("Company ID", placeholder="Company UUID")
        
        if st.button("Submit Disclosure", disabled=not all([title, description, user_id, creator_id, company_id])):
            try:
                result = _api_post(
                    "/api/v1/ip/disclosures",
                    {
                        "title": title,
                        "description": description,
                        "technical_field": technical_field,
                        "inventors": [{"user_id": user_id, "contribution_percentage": contribution, "is_primary": is_primary}],
                        "user_id": creator_id,
                        "company_id": company_id,
                    },
                )
                st.success(f"Disclosure created! ID: {result.get('id')}")
            except Exception as e:
                st.error(f"Failed to create disclosure: {e}")
    
    # List disclosures
    st.markdown("### Disclosure Portfolio")
    
    status_filter = st.selectbox(
        "Filter by Status",
        options=["All", "Draft", "Submitted", "Under Review", "Filed", "Granted", "Rejected"],
    )
    
    company_id_filter = st.text_input("Company ID", placeholder="UUID for filtering")
    
    if st.button("Load Disclosures") and company_id_filter:
        try:
            status_param = None if status_filter == "All" else status_filter.lower().replace(" ", "_")
            url = f"/api/v1/ip/disclosures?company_id={company_id_filter}"
            if status_param:
                url += f"&status={status_param}"
            
            disclosures = _api_get(url)
            st.session_state["disclosures"] = disclosures
        except Exception as e:
            st.error(f"Failed to load disclosures: {e}")
    
    disclosures = st.session_state.get("disclosures")
    
    if disclosures:
        df = pd.DataFrame(disclosures)
        st.dataframe(df, use_container_width=True, hide_index=True)


def render_patents_tab() -> None:
    """Render patents tab."""
    st.subheader("Patent Applications")
    
    # Filter controls
    jurisdiction_filter = st.selectbox("Jurisdiction", options=["All", "US", "EP", "PCT", "JP", "CN"])
    
    if st.button("Load Patents"):
        try:
            patents = _api_get("/api/v1/ip/patents")
            st.session_state["patents"] = patents
        except Exception as e:
            st.error(f"Failed to load patents: {e}")
    
    patents = st.session_state.get("patents")
    
    if patents:
        df = pd.DataFrame(patents)
        st.dataframe(df, use_container_width=True, hide_index=True)
        
        # View claims
        patent_ids = [p.get("id") for p in patents if p.get("id")]
        if patent_ids:
            selected_patent = st.selectbox("View Claims", patent_ids)
            
            if st.button("Load Claims"):
                try:
                    claims = _api_get(f"/api/v1/ip/patents/{selected_patent}/claims")
                    st.session_state["patent_claims"] = claims
                except Exception as e:
                    st.error(f"Failed to load claims: {e}")
    
    claims = st.session_state.get("patent_claims")
    
    if claims:
        st.markdown("#### Patent Claims")
        for claim in claims:
            claim_type_badge = "🔹 Independent" if claim.get("claim_type") == "independent" else "🔸 Dependent"
            st.markdown(f"**Claim {claim.get('claim_number')}** {claim_type_badge}")
            st.caption(claim.get("claim_text", "")[:200] + "...")


def render_evidence_links_tab() -> None:
    """Render evidence linking tab."""
    st.subheader("Evidence Links")
    
    # Link creation form
    with st.expander("🔗 Link Entity to Disclosure"):
        entity_type = st.selectbox("Entity Type", options=["Experiment", "Compound", "Dataset"], key="link_entity_type")
        entity_id = st.text_input("Entity ID", placeholder="UUID", key="link_entity_id")
        disclosure_id = st.text_input("Disclosure ID", placeholder="UUID", key="link_disclosure_id")
        link_type = st.selectbox("Link Type", options=["Supporting Evidence", "Prior Art", "Related Work"], key="link_type")
        notes = st.text_area("Notes", placeholder="Optional notes", key="link_notes")
        
        if st.button("Create Link", disabled=not all([entity_id, disclosure_id])):
            try:
                result = _api_post(
                    "/api/v1/ip/links",
                    {
                        "entity_type": entity_type.lower(),
                        "entity_id": entity_id,
                        "disclosure_id": disclosure_id,
                        "link_type": link_type.lower().replace(" ", "_"),
                        "notes": notes if notes else None,
                    },
                )
                st.success(f"Link created! ID: {result.get('id')}")
            except Exception as e:
                st.error(f"Failed to create link: {e}")
    
    # View entity links
    st.markdown("### Entity IP Links")
    
    col1, col2 = st.columns(2)
    
    with col1:
        view_entity_type = st.selectbox("Entity Type", options=["Experiment", "Compound", "Dataset"], key="view_entity_type")
    
    with col2:
        view_entity_id = st.text_input("Entity ID", key="view_entity_id")
    
    if st.button("Load Links") and view_entity_id:
        try:
            links = _api_get(f"/api/v1/ip/links/{view_entity_type.lower()}/{view_entity_id}")
            st.session_state["entity_links"] = links
        except Exception as e:
            st.error(f"Failed to load links: {e}")
    
    links = st.session_state.get("entity_links")
    
    if links:
        df = pd.DataFrame(links)
            st.dataframe(df, use_container_width=True, hide_index=True)
        else:
            st.info("No IP links found for this entity")


if __name__ == "__main__":
    render_ip_portfolio_page()

