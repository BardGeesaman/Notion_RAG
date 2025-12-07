"""Features page for the Streamlit dashboard."""

from __future__ import annotations

import pandas as pd
import streamlit as st

from amprenta_rag.analysis.cross_feature_mapping import get_cross_omics_feature_neighbors
from amprenta_rag.database.models import Feature
from scripts.dashboard.db_session import db_session


def render_features_page() -> None:
    """
    Render the Features page with filtering, search, and distribution views.

    Features:
    - Filter by feature type (gene, protein, metabolite, lipid)
    - Search features by name
    - Feature type distribution visualization
    - Export functionality (CSV)
    """
    st.header("🧪 Features")

    with db_session() as db:
        col1, col2 = st.columns(2)
        with col1:
            feature_type_filter = st.selectbox(
                "Filter by Feature Type",
                ["All", "gene", "protein", "metabolite", "lipid"],
            )
        with col2:
            search_term = st.text_input("Search features by name", "")

        query = db.query(Feature)

        if feature_type_filter != "All":
            query = query.filter(Feature.feature_type == feature_type_filter)

        if search_term:
            query = query.filter(Feature.name.ilike(f"%{search_term}%"))

        features = query.order_by(Feature.name).limit(100).all()

        st.metric("Features Found", len(features))
        if len(features) == 100:
            st.info("Showing first 100 features. Use filters to narrow down.")

        if features:
            # Group by feature type
            feature_data = []
            for feature in features:
                # Count datasets linked to this feature
                dataset_count = len(feature.datasets) if hasattr(feature, "datasets") and feature.datasets else 0

                feature_data.append(
                    {
                        "Name": feature.name,
                        "Type": feature.feature_type,
                        "Normalized Name": feature.normalized_name or "-",
                        "Linked Datasets": dataset_count,
                    }
                )
            df_features = pd.DataFrame(feature_data)
            st.dataframe(df_features, use_container_width=True)

            # Feature type distribution
            if len(df_features) > 0:
                type_counts = df_features["Type"].value_counts()
                if len(type_counts) > 1:
                    st.subheader("Feature Type Distribution")
                    col1, col2 = st.columns(2)
                    with col1:
                        st.bar_chart(type_counts)
                    with col2:
                        st.dataframe(type_counts.reset_index().rename(columns={"index": "Type", "Type": "Count"}))

                # Export button
                csv_features = df_features.to_csv(index=False)
                st.download_button(
                    label="📥 Download Features (CSV)",
                    data=csv_features,
                    file_name="features.csv",
                    mime="text/csv",
                )

            # Display feature details
            st.markdown("---")
            st.subheader("Feature Details")
            with db_session() as db:
                for feature in features:
                    with st.expander(f"**{feature.name}** ({feature.feature_type})"):
                        st.write(f"**ID:** `{feature.id}`")
                        st.write(f"**Type:** {feature.feature_type}")
                        if hasattr(feature, "normalized_name") and feature.normalized_name:
                            st.write(f"**Normalized:** {feature.normalized_name}")
                        st.write(f"**Linked Datasets:** {dataset_count}")
                        # Cross-omics neighbors panel
                        neighbors = get_cross_omics_feature_neighbors(feature.id, db)
                        st.write("**Cross-Omics Neighbors:**")
                        if neighbors.get("pathways"):
                            st.caption(f"Pathways: {', '.join(str(x) for x in neighbors['pathways'])}")
                        # Add genes ↔ proteins expansion as needed
                        if neighbors.get("genes"):
                            st.caption(f"Mapped Genes: {', '.join(neighbors['genes'])}")
                        if neighbors.get("proteins"):
                            st.caption(f"Mapped Proteins: {', '.join(neighbors['proteins'])}")
        else:
            st.info("No features found.")
