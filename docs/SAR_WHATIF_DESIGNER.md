# SAR What-If Designer

## Overview

The SAR What-If Designer provides an interactive workflow to:

- Validate SMILES (RDKit parsing)
- Render 2D structures (PNG)
- Predict ADMET endpoints (LogP, LogS, hERG) using the existing model registry predictor
- Compare multiple compounds side-by-side
- Generate “scaffold hop” variants via predefined transformations and re-score them

## API Endpoints

All endpoints are under `/api/v1/sar` (existing SAR router).

### Validate SMILES

```bash
curl -X POST "http://localhost:8000/api/v1/sar/validate" \
  -H "Content-Type: application/json" \
  -d '{"smiles":"c1ccccc1"}'
```

### Predict properties (batch)

```bash
curl -X POST "http://localhost:8000/api/v1/sar/predict" \
  -H "Content-Type: application/json" \
  -d '{"smiles_list":["c1ccccc1","CCO"]}'
```

Response includes:

- `smiles`
- `structure_img` (base64 PNG)
- `logp`, `logs` (regression value if model is registered)
- `herg` (probability if model is registered)

### List supported transformations

```bash
curl -X GET "http://localhost:8000/api/v1/sar/transformations"
```

### Scaffold hop

```bash
curl -X POST "http://localhost:8000/api/v1/sar/scaffold-hop" \
  -H "Content-Type: application/json" \
  -d '{"smiles":"c1ccccc1","transformation":"benzene_to_pyridine"}'
```

## Supported Transformations (MVP)

- `benzene_to_pyridine`: Replace phenyl ring with pyridine (one C→N).
- `cyclohexane_to_piperidine`: Replace cyclohexane with piperidine (one C→N).

## Dashboard

The Streamlit page **SAR What-If** includes:

1. **Property Prediction**:
   - Multi-line SMILES input
   - Results table with 2D structure and ADMET endpoints
   - Simple risk heuristic (GREEN/YELLOW/RED) and property bar chart

2. **Scaffold Hopping**:
   - Choose a transformation and generate variants
   - Predict properties for all generated variants


