# Visualization Guide

**Last Updated**: December 7, 2025

This document describes the interactive visualization capabilities of the Amprenta RAG Platform for exploring multi-omics data.

---

## Overview

The platform provides **four interactive visualization types** accessible through the Streamlit dashboard:

1. **Volcano Plot** - Differential expression analysis
2. **Heatmap** - Feature × Dataset intensity matrix
3. **PCA Scatter** - Dataset similarity and clustering
4. **Signature Network** - Feature co-occurrence networks

All visualizations are:
- **Interactive**: Zoom, pan, hover for details
- **Exportable**: Download as CSV or image
- **Real-time**: Query Postgres database on-demand
- **Customizable**: Adjust thresholds, colors, and filters

---

## Accessing Visualizations

### Start the Dashboard

```bash
# From project root
cd scripts/dashboard
streamlit run app.py

# Dashboard opens at: http://localhost:8501
```

### Navigate to Visualizations

1. **Open dashboard** in browser (`http://localhost:8501`)
2. **Select visualization** from sidebar:
   - Volcano Plot
   - Heatmap
   - PCA Scatter
   - Signature Network

---

## 1. Volcano Plot

### Purpose
Visualize differential expression results showing fold-change vs statistical significance.

### Use Cases
- **Identify significant features**: Find features with large fold-changes and low p-values
- **Quality check**: Verify expected patterns in differential expression data
- **Publication figures**: Export publication-ready volcano plots

### Data Requirements

**Dataset must contain**:
- Features with fold-change data
- Statistical significance values (p-values)

**Supported metadata fields** (in `Feature.external_ids`):
- `log2FC`, `log2fc`, `log2_fold_change`, `fold_change` (for X-axis)
- `pvalue`, `p_value`, `adj_pvalue`, `adj_p` (for Y-axis)

### How to Use

1. **Select Dataset**: Choose from dropdown (shows omics type)
2. **Select Feature Type**: 
   - `all` - All features
   - `gene` - Genes only
   - `protein` - Proteins only
   - `metabolite` - Metabolites only
   - `lipid` - Lipids only
3. **Set Thresholds**:
   - **Fold-change threshold**: Absolute log2FC cutoff (default: 1.0 = 2-fold)
   - **p-value threshold**: Significance cutoff (default: 0.05)
4. **Generate Plot**: Click "Generate Plot" button

### Interpretation

**Plot Elements**:
- **X-axis**: log2 Fold Change (negative = down-regulated, positive = up-regulated)
- **Y-axis**: -log10(p-value) (higher = more significant)
- **Red points**: Significant features (pass both thresholds)
- **Blue points**: Non-significant features
- **Dashed lines**: Threshold boundaries

**What to Look For**:
- **Top-right quadrant**: Significantly up-regulated features
- **Top-left quadrant**: Significantly down-regulated features
- **Bottom**: Non-significant features (low p-value, regardless of fold-change)
- **Center**: Small fold-changes (may be significant but biologically uninteresting)

**Example Interpretation**:
```
Red point at (log2FC=2.5, -log10p=4.5):
- Feature is up-regulated 2^2.5 = 5.7-fold
- p-value = 10^-4.5 = 0.00003 (highly significant)
- Conclusion: Strong up-regulation with high confidence
```

### Export Options

- **Download CSV**: Exports data table with columns: `feature`, `log2FC`, `pvalue`, `-log10p`, `significant`
- **Save Image**: Right-click plot → "Save image as..." (PNG format)

---

## 2. Heatmap

### Purpose
Compare feature intensities across multiple datasets with hierarchical clustering.

### Use Cases
- **Compare datasets**: Identify shared and differential features
- **Find patterns**: Discover feature clusters that co-occur
- **Visual comparison**: Quickly assess similarity between datasets

### Data Requirements

**Datasets must contain**:
- Features with intensity values

**Supported metadata fields** (in `Feature.external_ids`):
- `log2FC`, `log2fc`, `value`, `fold_change` (for intensity)
- Fallback: presence/absence (1.0 if feature present, 0.0 if absent)

### How to Use

1. **Select Datasets**: Choose 2+ datasets from multiselect dropdown
   - Default: First 3 datasets selected
   - Supports any number of datasets
2. **Generate Heatmap**: Click "Generate Heatmap" button

### Interpretation

**Plot Elements**:
- **Rows**: Features (automatically clustered)
- **Columns**: Datasets (automatically clustered)
- **Color**: Intensity value (Viridis scale: dark purple = low, bright yellow = high)
- **Clustering**: Similar features/datasets grouped together

**What to Look For**:
- **Horizontal bands**: Features that behave similarly across datasets
- **Vertical bands**: Datasets with similar feature profiles
- **Blocks**: Groups of features that co-occur in specific datasets
- **Outliers**: Individual high/low values (bright/dark cells)

**Clustering Algorithm**:
- **Method**: Average linkage hierarchical clustering
- **Distance**: Euclidean distance
- **Applied to**: Both rows (features) and columns (datasets)

**Example Interpretation**:
```
Bright yellow block in top-left corner:
- Group of features highly expressed in first 2 datasets
- Features are clustered together (similar patterns)
- Datasets are clustered together (similar profiles)
- Conclusion: Shared biological signature between these datasets
```

### Export Options

- **Download CSV**: Exports clustered matrix (rows = features, columns = dataset IDs)
- **Save Image**: Right-click plot → "Save image as..." (PNG format)

---

## 3. PCA Scatter

### Purpose
Visualize dataset similarity using Principal Component Analysis (PCA) dimensionality reduction.

### Use Cases
- **Dataset clustering**: Identify groups of similar datasets
- **Quality control**: Detect outliers or batch effects
- **Explore relationships**: Understand how datasets relate to each other
- **Color by metadata**: Visualize patterns by program, omics type, or disease

### Data Requirements

**Datasets must contain**:
- Features with values (used to build feature matrix)
- At least 2 datasets required for PCA

**Supported metadata fields** (in `Feature.external_ids`):
- `log2FC`, `log2fc`, `value`, `fold_change` (for feature vector)
- Fallback: presence/absence (1.0 if present)

### How to Use

1. **Select Color By**: Choose grouping variable
   - `program` - Color by research program
   - `omics_type` - Color by lipidomics/metabolomics/proteomics/transcriptomics
   - `disease` - Color by disease context
2. **Select Program** (optional): Filter datasets to specific program
3. **Select Datasets**: Choose datasets from multiselect dropdown
   - Default: First 5 datasets selected
   - Minimum: 2 datasets required
4. **Run PCA**: Click "Run PCA" button

### Interpretation

**Plot Elements**:
- **X-axis**: Principal Component 1 (PC1) - largest variance direction
- **Y-axis**: Principal Component 2 (PC2) - second largest variance direction
- **Points**: Individual datasets
- **Colors**: Groups (based on selected color variable)
- **Hover**: Dataset ID and metadata

**What to Look For**:
- **Clusters**: Groups of nearby points indicate similar datasets
- **Outliers**: Isolated points may be unusual or erroneous datasets
- **Color patterns**: Same-color clustering suggests metadata explains variation
- **Axes percentages**: Higher % = more variance explained by that component

**Example Interpretation**:
```
PC1 (45%), PC2 (22%):
- PC1 explains 45% of variance (major pattern)
- PC2 explains 22% of variance (secondary pattern)
- Together: 67% of total variance captured

Cluster of red points (ALS datasets) in top-right:
- ALS datasets are similar to each other
- Separated from blue points (PD datasets)
- Conclusion: ALS signature distinguishable from PD
```

**Statistical Notes**:
- PCA centers and scales data automatically
- Variance explained shown in axis labels
- 2 components shown (can't visualize higher dimensions)

### Export Options

- **Download PCA CSV**: Exports coordinates with columns: `PC1`, `PC2`, `dataset_id`, `color` (grouping variable)
- **Save Image**: Right-click plot → "Save image as..." (PNG format)

---

## 4. Signature Network

### Purpose
Visualize feature co-occurrence within a signature as a network graph.

### Use Cases
- **Understand signatures**: See which features are part of a signature
- **Visual representation**: Network diagram for presentations/publications
- **Feature relationships**: Understand signature composition at a glance

### Data Requirements

**Signature must contain**:
- At least 1 feature
- Features stored in `Signature.features` relationship

### How to Use

1. **Select Signature**: Choose from dropdown (shows signature names)
2. **Build Network**: Click "Build Network" button

### Interpretation

**Plot Elements**:
- **Nodes**: Individual features (genes, proteins, metabolites, lipids)
- **Edges**: Connections showing co-occurrence in signature
- **Labels**: Feature names displayed above nodes
- **Layout**: Grid layout (6 nodes per row)

**What to Look For**:
- **Number of nodes**: Signature size (how many features)
- **Feature names**: Identify which molecules are in the signature
- **Connectivity**: Fully connected network (all features part of same signature)

**Network Structure**:
- **Fully connected**: Every feature connected to every other feature
- **Rationale**: All features co-occur in the same signature by definition
- **Alternative interpretation**: Can be modified to show feature-feature correlations (future enhancement)

**Example Interpretation**:
```
Network with 6 nodes:
- Cer(d18:1/16:0), Cer(d18:1/18:0), Cer(d18:1/24:1)
- SM(d18:1/16:0), SM(d18:1/18:0), HexCer(d18:1/16:0)

Conclusion: 
- Signature contains 6 ceramide/sphingolipid species
- Mix of ceramides, sphingomyelins, and hexosylceramides
- All species co-occur in signature definition
```

### Export Options

- **Save Image**: Right-click plot → "Save image as..." (PNG format)
- **CSV export**: Not applicable (network structure, not tabular data)

---

## General Features

### Interactive Controls

All visualizations support:
- **Zoom**: Scroll wheel or click-drag with box tool
- **Pan**: Click-drag with hand tool
- **Hover**: Mouse over points to see details
- **Reset**: Double-click to reset view
- **Download**: Camera icon to save as PNG

### Refresh Data

All visualization pages have **"Refresh data"** button:
- Reloads dataset/signature list from Postgres
- Use after ingesting new data
- Forces dashboard to query database again

### Performance

**Data Source**: All visualizations query Postgres in real-time
- **Volcano Plot**: ~1-5 seconds (depends on feature count)
- **Heatmap**: ~2-10 seconds (depends on dataset count and features)
- **PCA**: ~3-15 seconds (depends on feature matrix size)
- **Network**: ~1-3 seconds (depends on signature size)

**Optimization Tips**:
- Limit dataset selection for heatmap/PCA (< 50 datasets recommended)
- Filter by feature type for volcano plots
- Large signatures (>100 features) may render slowly in network view

---

## Data Format Requirements

### Feature Metadata Structure

For visualizations to work, features must have metadata in `Feature.external_ids` (JSON field):

**Example Transcriptomics/Proteomics**:
```json
{
  "log2FC": 2.5,
  "pvalue": 0.0001,
  "adj_pvalue": 0.001
}
```

**Example Lipidomics/Metabolomics**:
```json
{
  "fold_change": 3.2,
  "p_value": 0.01,
  "value": 1250.5
}
```

**Fallback Behavior**:
- If statistical fields missing: Volcano plot shows warning
- If intensity fields missing: Heatmap/PCA use presence/absence (1.0 if feature exists)

### Dataset Metadata

Datasets should have:
- `name` - Display name
- `omics_type` - For filtering and coloring
- `disease` - Array of disease contexts (for PCA coloring)
- `programs` - Relationship to Program table (for PCA coloring)

---

## Troubleshooting

### Issue: "No feature statistics found"

**Problem**: Volcano plot shows warning, no data plotted

**Solutions**:
1. Check feature metadata has `log2FC`/`pvalue` fields
2. Verify feature type filter matches your data
3. Ensure dataset has features with statistical results

**Diagnostic**:
```python
from amprenta_rag.database.base import get_db
from amprenta_rag.database.models import Dataset

db = next(get_db())
dataset = db.query(Dataset).filter(Dataset.id == "your-id").first()
for feat in dataset.features[:5]:  # Check first 5 features
    print(f"{feat.name}: {feat.external_ids}")
```

### Issue: "No features found for selected datasets"

**Problem**: Heatmap or PCA shows warning, no plot generated

**Solutions**:
1. Verify selected datasets have features
2. Check datasets are not empty
3. Ensure feature-dataset relationships exist in Postgres

**Diagnostic**:
```python
dataset = db.query(Dataset).filter(Dataset.id == "your-id").first()
print(f"Feature count: {len(dataset.features)}")
```

### Issue: Plot is too slow to generate

**Problem**: Visualization takes >30 seconds to render

**Solutions**:
1. **Reduce dataset count**: Select fewer datasets for heatmap/PCA
2. **Filter by feature type**: Limit to specific omics type for volcano plot
3. **Limit signatures**: Select signatures with <50 features for network view

**Performance Targets**:
- Heatmap: < 20 datasets, < 1000 features per dataset
- PCA: < 50 datasets
- Volcano: Any single dataset (fast)
- Network: < 100 features in signature

### Issue: Wrong features shown in plot

**Problem**: Features don't match expectations

**Solutions**:
1. Check feature type filter (volcano plot)
2. Verify dataset-feature relationships in Postgres
3. Ensure correct dataset selected

---

## Advanced Usage

### Export for Publications

**High-Resolution Images**:
1. Generate plot in dashboard
2. Use Plotly camera icon (top-right of plot)
3. Choose "Download plot as PNG"
4. Resolution: 1200x800px (suitable for presentations)

**For higher resolution**:
- Right-click → "Save image as SVG" (if available)
- Use external tools (Inkscape, Illustrator) to edit SVG

**CSV Data for Custom Plots**:
1. Export CSV from dashboard
2. Import into R, Python, or Excel
3. Create custom plots with ggplot2, matplotlib, etc.

### Programmatic Access

**Generate plots in Python**:

```python
from amprenta_rag.database.base import get_db
from amprenta_rag.database.models import Dataset
import plotly.express as px
import numpy as np
import pandas as pd

# Get dataset
db = next(get_db())
dataset = db.query(Dataset).filter(Dataset.name == "ST004396").first()

# Extract volcano data
rows = []
for feat in dataset.features:
    ext = feat.external_ids or {}
    if "log2FC" in ext and "pvalue" in ext:
        rows.append({
            "feature": feat.name,
            "log2FC": float(ext["log2FC"]),
            "pvalue": float(ext["pvalue"])
        })

# Create plot
df = pd.DataFrame(rows)
df["-log10p"] = -np.log10(df["pvalue"])
fig = px.scatter(df, x="log2FC", y="-log10p", hover_name="feature")
fig.show()
```

### Custom Visualizations

Want to create custom plots? Use these modules:
- `scripts/dashboard/pages/visualizations/` - Copy and modify existing visualizations
- `amprenta_rag.database.models` - Access Postgres data
- `plotly` - Create interactive plots
- `streamlit` - Build dashboard pages

**Example: Add new visualization**:
1. Create `scripts/dashboard/pages/visualizations/my_plot.py`
2. Implement `render()` function
3. Add to dashboard navigation

---

## Future Enhancements

### Planned Visualizations (Q1 2026)

1. **Box Plots**: Compare feature distributions across groups
2. **Time Series**: Track features over time (longitudinal data)
3. **Pathway Maps**: KEGG/Reactome pathway diagrams with highlighted features
4. **Correlation Networks**: Feature-feature correlation networks (not just co-occurrence)
5. **3D PCA**: Three principal components for richer visualization

### Planned Features

- **Save plot configurations**: Bookmark favorite plot settings
- **Batch export**: Export all plots for a program at once
- **Interactive filtering**: Click points to filter/highlight
- **Annotations**: Add text labels and arrows to plots
- **Themes**: Light/dark mode, color-blind friendly palettes

---

## See Also

- [Streamlit Dashboard Guide](STREAMLIT_DASHBOARD_GUIDE.md) - Dashboard overview
- [API Reference](API_REFERENCE.md) - Programmatic data access
- [Usage Examples](USAGE_EXAMPLES.md) - Data ingestion workflows
- [Database Models](../amprenta_rag/database/models.py) - Postgres schema

---

**Questions or Issues?**

See [Troubleshooting Guide](TROUBLESHOOTING.md) or file an issue on GitHub.

