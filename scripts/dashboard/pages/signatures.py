"""Signatures page for the Streamlit dashboard."""

from __future__ import annotations

import pandas as pd
import streamlit as st

from amprenta_rag.database.models import Signature
from amprenta_rag.signatures.signature_validation import validate_signature_against_all_datasets
from scripts.dashboard.db_session import db_session


def render_signatures_page() -> None:
    """
    Render the Signatures page with search, summary, and detailed views.

    Features:
    - Search signatures by name
    - Summary table with export
    - Modalities distribution visualization
    - Detailed signature views
    """
    st.header("📋 Signatures")

    with db_session() as db:
        search_term = st.text_input("Search signatures by name", "")

        query = db.query(Signature)
        if search_term:
            query = query.filter(Signature.name.ilike(f"%{search_term}%"))

        signatures = query.order_by(Signature.created_at.desc()).all()

        st.metric("Total Signatures", len(signatures))

        if signatures:
            # Summary table
            signature_data = []
            for sig in signatures:
                validated = validate_signature_against_all_datasets(sig.id)
                cov = validated.metrics.coverage
                mean_score = validated.metrics.mean_score or 0.0
                # quality_label: simple threshold
                if cov >= 0.7:
                    quality = "High"
                elif cov >= 0.4:
                    quality = "Medium"
                elif cov > 0:
                    quality = "Low"
                else:
                    quality = "None"
                signature_data.append(
                    {
                        "Name": sig.name,
                        "Description": sig.description or "",
                        "Modalities": ", ".join(sig.modalities) if sig.modalities else "",
                        "Components": len(sig.components),
                        "Coverage": cov,
                        "Mean Score": mean_score,
                        "Quality": quality,
                        "Created": sig.created_at.strftime("%Y-%m-%d"),
                    }
                )
            df_signatures = pd.DataFrame(signature_data)
            st.dataframe(df_signatures, use_container_width=True, hide_index=True)

            # Modalities distribution
            if len(df_signatures) > 0 and any(df_signatures["Modalities"]):
                st.subheader("Signature Modalities Distribution")
                # Count modalities
                all_modalities = []
                for mods in df_signatures["Modalities"]:
                    if mods:
                        all_modalities.extend([m.strip() for m in mods.split(",")])
                if all_modalities:
                    modality_counts = pd.Series(all_modalities).value_counts()
                    st.bar_chart(modality_counts)

            # Export button
            csv_signatures = df_signatures.to_csv(index=False)
            st.download_button(
                label="📥 Download Signatures (CSV)",
                data=csv_signatures,
                file_name="signatures.csv",
                mime="text/csv",
            )

            st.markdown("---")
            st.subheader("Signature Details")

            for signature in signatures:
                validated = validate_signature_against_all_datasets(signature.id)
                with st.expander(f"**{signature.name}**"):
                    col1, col2 = st.columns(2)
                    with col1:
                        st.write(f"**ID:** `{signature.id}`")
                        if signature.description:
                            st.write(f"**Description:** {signature.description}")
                        if signature.modalities:
                            st.write(f"**Modalities:** {', '.join(signature.modalities)}")
                    with col2:
                        st.write(f"**Created:** {signature.created_at.strftime('%Y-%m-%d %H:%M')}")
                        component_count = len(signature.components)
                        st.write(f"**Components:** {component_count}")
                    # Inline validation/quality summary
                    m = validated.metrics
                    st.markdown(
                        f"**Validation:** Coverage {m.num_matched_datasets}/{m.num_total_datasets} ({m.coverage:.1%}); Mean Score = {m.mean_score:.2f if m.mean_score else 0.0}\n\n{validated.summary}"
                    )
        else:
            st.info("No signatures found.")
