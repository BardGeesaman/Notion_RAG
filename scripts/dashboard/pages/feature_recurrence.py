from collections import defaultdict

import pandas as pd
import streamlit as st

from amprenta_rag.database.models import Feature
from scripts.dashboard.auth import require_auth
from scripts.dashboard.db_session import db_session


def render_feature_recurrence_page():
    require_auth()
    st.header("Feature Recurrence")
    st.info(
        """The recurrence counts how many datasets each feature is found in. Highly-recurrent features may be core biology or technical background. Use filters to differentiate."""
    )
    feature_type = st.selectbox("Feature type", ["gene", "protein", "metabolite", "lipid"])
    min_rec = st.slider("Min recurrence (datasets)", min_value=1, max_value=20, value=2)
    disease_filter = st.text_input("Disease filter (optional)", value="")
    with db_session() as db:
        feats = db.query(Feature).filter(Feature.feature_type == feature_type).all()
        # Build feature -> set of dataset ids
        feat_to_ds = defaultdict(set)
        for f in feats:
            for ds in f.datasets if hasattr(f, "datasets") else []:
                if not disease_filter or (hasattr(ds, "disease") and disease_filter in (ds.disease or [])):
                    feat_to_ds[f.name].add(ds.id)
        rows = []
        for fname, dsset in feat_to_ds.items():
            rows.append({"feature_name": fname, "n_datasets": len(dsset), "feature_type": feature_type})
        rows = [r for r in rows if r["n_datasets"] >= min_rec]
        if rows:
            df = pd.DataFrame(rows).sort_values("n_datasets", ascending=False)
            st.bar_chart(df.set_index("feature_name")["n_datasets"].head(50))
            st.dataframe(df, width='stretch')
            st.caption("Bar shows top 50 features by dataset recurrence.")
        else:
            st.info("No features at min threshold.")
