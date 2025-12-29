---
name: High-Dimensional Projector
overview: Add interactive 3D visualization for high-dimensional data using UMAP/t-SNE projections with Plotly, supporting any dataset's feature matrix.
todos:
  - id: proj-batch-1
    content: Core Projection Engine (projection_engine.py)
    status: completed
  - id: proj-batch-2
    content: API Endpoints (projector router)
    status: completed
    dependencies:
      - proj-batch-1
  - id: proj-batch-3
    content: Dashboard Page (projector.py)
    status: completed
    dependencies:
      - proj-batch-2
  - id: proj-batch-4
    content: Tests (unit, API, E2E)
    status: completed
    dependencies:
      - proj-batch-3
---

# High-Dimensional Projector

## Overview

Create an interactive 3D visualization tool for exploring high-dimensional datasets using UMAP and t-SNE dimensionality reduction, with point cloud rendering via Plotly.

## Current State

Existing related code:
- `amprenta_rag/single_cell/scanpy_pipeline.py` - UMAP for single-cell only
- `scripts/dashboard/pages/visualizations/pca.py` - PCA visualization
- `amprenta_rag/analysis/timeseries.py` - K-means clustering

Gap: No general-purpose projector for arbitrary feature matrices (proteomics, metabolomics, compound descriptors, signatures).

## Architecture

```
Input Sources:
- Dataset feature matrices (any omics type)
- Compound descriptor matrices (Morgan FP, RDKit)
- Signature vectors
- Custom uploaded matrices

Pipeline:
[Feature Matrix] → [Preprocessing] → [UMAP/t-SNE] → [3D Plotly] → [Interactive Dashboard]

Features:
- Algorithm selection (UMAP, t-SNE, PCA)
- Parameter tuning (n_neighbors, perplexity, n_components)
- Color by metadata (cluster, label, continuous value)
- Point size by metric
- Hover tooltips with entity details
- Export coordinates
```

## Batch Organization

### Batch 1: Core Projection Engine

**Dependency to add:**
- `umap-learn>=0.5.4` to requirements.txt

Files to create:
- `amprenta_rag/analysis/projection_engine.py`

Features:
- ProjectorEngine class with UMAP, t-SNE, PCA methods
- Input: numpy array or pandas DataFrame
- Output: ProjectionResult with coordinates + metadata
- Caching by input hash + parameters
- Support 2D and 3D projections

Pydantic schemas:
- ProjectionParams: algorithm, n_components, n_neighbors, perplexity, random_state
- ProjectionResult: coordinates, explained_variance, algorithm_used, params

### Batch 2: API Endpoints

Files to modify/create:
- `amprenta_rag/api/routers/projector.py` (new router)
- `amprenta_rag/api/main.py` (register router)

Endpoints:
- POST /api/v1/projector/compute - Compute projection from dataset_id or uploaded matrix
- GET /api/v1/projector/datasets - List datasets available for projection
- POST /api/v1/projector/export - Export projection coordinates as CSV

### Batch 3: Dashboard Page

Files to create:
- `scripts/dashboard/pages/projector.py`

UI Features:
- Dataset selector (dropdown of available datasets)
- Algorithm selector (UMAP, t-SNE, PCA)
- Parameter sliders (n_neighbors, perplexity)
- 2D/3D toggle
- Plotly scatter with:
  - Color by cluster/label/value
  - Hover with entity details
  - Zoom/pan/rotate
- Export button

### Batch 4: Tests

Files to create:
- `amprenta_rag/tests/analysis/test_projection_engine.py`
- `amprenta_rag/tests/api/test_projector_api.py`
- `amprenta_rag/tests/e2e/test_projector_page.py`

Coverage:
- Projection algorithms work correctly
- API endpoints return valid responses
- Dashboard loads and displays projections

## Success Criteria

- UMAP, t-SNE, PCA projections working
- 2D and 3D visualization with Plotly
- Works with any dataset feature matrix
- 20+ tests passing
- Zero skipped tests

## P3 Items - COMPLETED (No Deferral Policy)

1. **Real Dataset Feature Matrix Loading** ✅
   - Implemented: Loads actual features from Dataset model
   - Commit: 35ebab1

2. **Color-by Metadata in UI** ✅
   - Implemented: Color By selector (None, Index)
   - Commit: 35ebab1
