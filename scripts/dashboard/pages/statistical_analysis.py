from __future__ import annotations

import pandas as pd
import streamlit as st
from sqlalchemy.orm import Session

from amprenta_rag.analysis import statistical_tests as stests
from amprenta_rag.database.models import Dataset, Feature, dataset_feature_assoc
from scripts.dashboard.db_session import db_session


def _list_datasets(db: Session):
    return [
        {"id": str(d.id), "name": d.name, "omics_type": d.omics_type or "unknown"}
        for d in db.query(Dataset).order_by(Dataset.updated_at.desc()).limit(200).all()
    ]


def _feature_value(feat: Feature):
    ext = feat.external_ids or {}
    for key in ("log2FC", "log2fc", "value", "fold_change"):
        if key in ext and ext[key] is not None:
            try:
                return float(ext[key])
            except Exception:
                continue
    return None


def _collect_feature_values(db: Session, dataset_ids, feature_name: str):
    """
    Return mapping dataset_id -> value for the given feature name.
    """
    if not dataset_ids or not feature_name:
        return {}
    rows = (
        db.query(dataset_feature_assoc.c.dataset_id, Feature)
        .join(Feature, Feature.id == dataset_feature_assoc.c.feature_id)
        .filter(dataset_feature_assoc.c.dataset_id.in_(dataset_ids))
        .filter(Feature.name == feature_name)
        .all()
    )
    values = {}
    for ds_id, feat in rows:
        val = _feature_value(feat)
        if val is not None:
            values[str(ds_id)] = val
    return values


def render_statistical_analysis_page() -> None:
    st.title("📊 Statistical Analysis")
    st.caption("Built-in tests: T-test, ANOVA, Mann-Whitney, Correlation.")

    test_type = st.selectbox(
        "Test type",
        ["T-test", "ANOVA", "Mann-Whitney", "Correlation"],
        key="stats_test_type",
    )

    correction = st.selectbox(
        "Multiple testing correction",
        ["fdr_bh", "bonferroni", "holm", "sidak", "none"],
        index=0,
        key="stats_correction",
    )

    with db_session() as db:
        datasets = _list_datasets(db)

    if not datasets:
        st.info("No datasets found.")
        return

    ds_options = {f"{d['name']} ({d['omics_type']})": d["id"] for d in datasets}
    selected = st.multiselect(
        "Datasets",
        list(ds_options.keys()),
        default=list(ds_options.keys())[:2],
        key="stats_datasets",
    )
    selected_ids = [ds_options[s] for s in selected]

    feature_x = st.text_input("Feature X (e.g., TP53)", key="stats_feature_x")
    feature_y = st.text_input("Feature Y (for correlation, optional)", key="stats_feature_y")

    default_split = max(1, len(selected_ids) // 2) if selected_ids else 1
    group_a_labels = st.multiselect(
        "Group A datasets",
        list(ds_options.keys()),
        default=list(ds_options.keys())[:default_split],
        key="stats_group_a",
    )
    group_b_labels = st.multiselect(
        "Group B datasets",
        list(ds_options.keys()),
        default=list(ds_options.keys())[default_split:default_split * 2],
        key="stats_group_b",
    )
    group_a_ids = [ds_options[lbl] for lbl in group_a_labels]
    group_b_ids = [ds_options[lbl] for lbl in group_b_labels]

    if st.button("Run test", key="stats_run"):
        if test_type == "Correlation":
            if not feature_x or not feature_y:
                st.warning("Enter both Feature X and Feature Y for correlation.")
                return
        else:
            if not feature_x:
                st.warning("Enter Feature X.")
                return

        if test_type != "Correlation" and (not group_a_ids or not group_b_ids):
            st.warning("Select at least one dataset for both Group A and Group B.")
            return

        with st.spinner("Running statistical test..."):
            results = []
            with db_session() as db:
                if test_type == "Correlation":
                    vals_x = _collect_feature_values(db, selected_ids, feature_x)
                    vals_y = _collect_feature_values(db, selected_ids, feature_y)
                    common_ids = list(set(vals_x.keys()) & set(vals_y.keys()))
                    if len(common_ids) < 2:
                        st.warning("Not enough overlapping datasets with both features.")
                        return
                    x = [vals_x[i] for i in common_ids]
                    y = [vals_y[i] for i in common_ids]
                    res = stests.pearson_corr(x, y)
                    results.append({"feature_x": feature_x, "feature_y": feature_y, **res})
                else:
                    g1_vals = _collect_feature_values(db, group_a_ids, feature_x)
                    g2_vals = _collect_feature_values(db, group_b_ids, feature_x)
                    g1 = list(g1_vals.values())
                    g2 = list(g2_vals.values())
                    if not g1 or not g2:
                        st.warning("Groups missing values for the selected feature.")
                        return
                    if test_type == "T-test":
                        res = stests.ttest_independent(g1, g2)
                    elif test_type == "ANOVA":
                        res = stests.anova_oneway(g1, g2)
                    else:
                        res = stests.mann_whitney(g1, g2)
                    results.append({"feature": feature_x, **res})

            df = pd.DataFrame(results)
            if correction != "none" and not df.empty:
                adj, reject = stests.adjust_pvalues(df["pvalue"].tolist(), method=correction)
                df["pvalue_adj"] = adj
                df["reject"] = reject

            st.dataframe(df, width='stretch', hide_index=True)
            st.download_button(
                "Export CSV",
                data=df.to_csv(index=False),
                file_name="stats_results.csv",
                mime="text/csv",
                key="stats_export",
            )

