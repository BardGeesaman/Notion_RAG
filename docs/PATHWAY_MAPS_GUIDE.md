# Pathway Maps Guide

## Overview
**Pathway Maps** provide interactive KEGG pathway diagrams with data overlays so you can interpret omics signals in biological context. The MVP is **KEGG-only**: it fetches pathway structure from **KGML**, renders nodes/edges in Cytoscape, and optionally colors gene nodes using expression metrics (e.g., log2 fold-change).

Key capabilities:
- Fetch and cache KEGG pathway structure (nodes, edges, layout coordinates)
- Render pathways with KGML-derived positions or force-directed layouts
- Overlay expression values onto gene nodes with diverging colormaps
- Export pathway snapshots to PNG from the viewer

## Browsing Pathways
The **Browse Pathways** tab supports two discovery paths:
- **Search KEGG**: enter a query (pathway ID fragment or name keyword) and call the KEGG “find pathway” endpoint. Results include a `pathway_id`, pathway name, and organism code.
- **Enriched Pathways (from dataset)**: select a dataset and fetch the top enriched KEGG pathways computed by the backend enrichment routine. This is useful when you want to go directly from a dataset to candidate pathways for exploration.

Selecting a pathway stores the `pathway_id` in `session_state` so you can switch tabs without losing context.

## Viewing Pathway Maps
The **Pathway Map** tab renders the selected pathway using Cytoscape with KEGG-friendly defaults:
- **Node shapes by type**:
  - gene → rectangle
  - compound → ellipse
  - sub-pathway (map links) → hexagon
- **Edge styles**:
  - activation/expression → green solid
  - inhibition/repression → red dashed
  - other → gray solid

Layout options:
- **KEGG Layout** uses the KGML `graphics x/y` coordinates (normalized and scaled) so the map resembles KEGG’s canonical diagram.
- **Force-Directed** (cose) lets you re-arrange the network to reduce overlap or explore connectivity.

## Data Overlay
The **Data Overlay** tab lets you color gene nodes based on dataset expression:
1. Choose a dataset and click **Load Expression** to retrieve a `gene_expression` mapping (gene_symbol → fold change).
2. Configure the visualization range with `vmin`/`vmax` (default -2 to 2 for log2FC).
3. Pick a colormap (e.g., `RdBu_r`, `coolwarm`, `viridis`) and click **Apply Overlay**.

Colors are computed with a diverging normalization centered at 0 (TwoSlopeNorm), so negative values map to one side of the palette and positive values to the other. A legend is shown to help interpret the scale.

## API Reference
All endpoints are under `/api/pathway-maps`:

### 1) Get pathway structure
`GET /api/pathway-maps/structure/{pathway_id}`

### 2) Compute overlay for a pathway
`POST /api/pathway-maps/overlay`

Example request:
```json
{
  "pathway_id": "hsa00010",
  "expression_data": {"TP53": 1.2, "BRCA1": -0.8},
  "colormap": "RdBu_r",
  "vmin": -2.0,
  "vmax": 2.0
}
```

### 3) Search KEGG pathways
`GET /api/pathway-maps/search?query=glycolysis`

### 4) Top enriched pathways for a dataset
`GET /api/pathway-maps/enriched/{dataset_id}`

### 5) Expression values for overlay
`GET /api/pathway-maps/expression/{dataset_id}`

## Technical Notes
- **KGML parsing**: the service parses KGML `<entry>` elements into nodes and `<relation>/<reaction>` elements into edges. Reaction edges are modeled as substrate → product links.
- **Coordinate normalization**: KGML x/y positions are normalized by max values so the renderer can scale consistently.
- **Caching**: pathway structures are cached on disk in `data/pathway_cache/{pathway_id}.json` with a **7-day TTL** to reduce repeated KEGG requests.
- **Rate limiting**: KEGG requests are throttled with a 0.5s delay between calls (structure fetch and search).


