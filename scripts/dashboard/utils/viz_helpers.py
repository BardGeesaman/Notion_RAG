"""
Visualization helpers for Amprenta dashboard.

Provides consistent theming, color palettes, and chart factory functions.
"""

from typing import Any, Dict, List, Optional

import pandas as pd
import plotly.express as px
import plotly.graph_objects as go
from sqlalchemy.orm import Session

# =============================================================================
# COLOR PALETTES
# =============================================================================

AMPRENTA_COLORS = {
    "primary": "#2E86AB",
    "secondary": "#A23B72",
    "success": "#28A745",
    "warning": "#F18F01",
    "danger": "#C73E1D",
    "neutral": "#6C757D",
    "light": "#E9ECEF",
    "dark": "#212529",
}

OMICS_COLORS = {
    "transcriptomics": "#4E79A7",
    "proteomics": "#F28E2B",
    "metabolomics": "#59A14F",
    "lipidomics": "#E15759",
    "unknown": "#76B7B2",
}

SIGNIFICANCE_COLORS = {
    True: "#FF6B6B",
    False: "#4E79A7",
}

CATEGORICAL_PALETTE = [
    "#4E79A7",
    "#F28E2B",
    "#E15759",
    "#76B7B2",
    "#59A14F",
    "#EDC948",
    "#B07AA1",
    "#FF9DA7",
    "#9C755F",
    "#BAB0AC",
]

# =============================================================================
# PLOTLY THEME
# =============================================================================


def get_amprenta_template() -> Dict[str, Any]:
    """Get Amprenta-branded Plotly layout template."""
    return {
        "layout": {
            "font": {"family": "Inter, -apple-system, sans-serif", "size": 12, "color": "#E9ECEF"},
            "title": {"font": {"size": 16, "color": "#FFFFFF"}},
            "paper_bgcolor": "rgba(0,0,0,0)",
            "plot_bgcolor": "rgba(30,30,30,0.5)",
            "colorway": CATEGORICAL_PALETTE,
            "xaxis": {"gridcolor": "rgba(255,255,255,0.1)", "zerolinecolor": "rgba(255,255,255,0.2)"},
            "yaxis": {"gridcolor": "rgba(255,255,255,0.1)", "zerolinecolor": "rgba(255,255,255,0.2)"},
            "legend": {"bgcolor": "rgba(0,0,0,0)", "font": {"color": "#E9ECEF"}},
        }
    }


def apply_amprenta_theme(fig: go.Figure) -> go.Figure:
    """Apply Amprenta theme to a Plotly figure."""
    template = get_amprenta_template()
    fig.update_layout(**template["layout"])
    return fig


# =============================================================================
# CHART FACTORIES
# =============================================================================


def create_volcano_plot(
    df: pd.DataFrame,
    fc_col: str,
    pval_col: str,
    fc_thresh: float = 1.0,
    p_thresh: float = 0.05,
) -> go.Figure:
    """Create a volcano plot with threshold lines."""
    fig = px.scatter(
        df,
        x=fc_col,
        y=pval_col,
        color=df[pval_col] < p_thresh,
        color_discrete_map={True: SIGNIFICANCE_COLORS[True], False: SIGNIFICANCE_COLORS[False]},
        hover_data=df.columns,
    )
    fig.add_hline(y=p_thresh, line_dash="dot", line_color="#888")
    fig.add_vline(x=fc_thresh, line_dash="dot", line_color="#888")
    fig.add_vline(x=-fc_thresh, line_dash="dot", line_color="#888")
    fig.update_layout(title="Volcano Plot", xaxis_title="log2 Fold Change", yaxis_title=pval_col)
    return apply_amprenta_theme(fig)


def create_heatmap(
    matrix: pd.DataFrame,
    title: str = "Heatmap",
    colorscale: str = "RdBu",
    cluster: bool = True,
) -> go.Figure:
    """Create a (clustered) heatmap from a matrix DataFrame."""
    data = matrix.copy()
    if cluster and len(data) > 1 and len(data.columns) > 1:
        # Simple clustering via correlation ordering
        row_order = data.corr().mean(axis=1).sort_values().index
        col_order = data.T.corr().mean(axis=1).sort_values().index
        data = data.loc[row_order, col_order]
    fig = px.imshow(
        data,
        color_continuous_scale=colorscale,
        aspect="auto",
        origin="lower",
    )
    fig.update_layout(title=title)
    return apply_amprenta_theme(fig)


def create_scatter(
    df: pd.DataFrame,
    x: str,
    y: str,
    color_by: Optional[str] = None,
    title: str = "Scatter Plot",
) -> go.Figure:
    """Create a basic scatter plot."""
    fig = px.scatter(
        df,
        x=x,
        y=y,
        color=color_by,
        color_discrete_sequence=CATEGORICAL_PALETTE,
        hover_data=df.columns,
        title=title,
    )
    return apply_amprenta_theme(fig)


def create_bar_chart(
    df: pd.DataFrame,
    x: str,
    y: str,
    color_by: Optional[str] = None,
    title: str = "Bar Chart",
    orientation: str = "v",
) -> go.Figure:
    """Create a bar chart."""
    fig = px.bar(
        df,
        x=x if orientation == "v" else y,
        y=y if orientation == "v" else x,
        color=color_by,
        color_discrete_sequence=CATEGORICAL_PALETTE,
        title=title,
        orientation=orientation,
    )
    return apply_amprenta_theme(fig)


# =============================================================================
# DATA HELPERS
# =============================================================================


def list_datasets_for_dropdown(db: Session, limit: int = 200) -> List[Dict[str, str]]:
    """Return list of dataset dicts for dropdowns."""
    from amprenta_rag.database.models import Dataset

    rows = (
        db.query(Dataset)
        .order_by(Dataset.updated_at.desc())
        .limit(limit)
        .all()
    )
    return [
        {
            "id": str(d.id) if d.id is not None else "",
            "name": str(d.name) if d.name is not None else "",
            "omics_type": str(d.omics_type) if d.omics_type is not None else "unknown",
        }
        for d in rows
        if d.id is not None and d.name is not None
    ]


def extract_feature_stats(feature) -> Dict[str, Optional[float]]:
    """Extract common stat fields from feature.external_ids."""
    stats = feature.external_ids or {}
    log2fc = stats.get("log2FoldChange") or stats.get("log2FC") or stats.get("fc")
    pval = stats.get("pvalue") or stats.get("p_value") or stats.get("pval") or stats.get("padj")
    return {"log2fc": log2fc, "pval": pval}


def get_feature_value(feature, default: float = 1.0) -> float:
    """Extract numeric value for matrices; fallback to default."""
    stats = feature.external_ids or {}
    for key in ["log2FoldChange", "log2FC", "fc", "value"]:
        if key in stats:
            try:
                return float(stats[key])
            except Exception:
                continue
    return default

