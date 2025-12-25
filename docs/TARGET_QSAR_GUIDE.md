# Target-Activity QSAR Guide

## Overview
Target-Activity QSAR provides **per-target activity prediction** as a binary classifier:
- **Active** if predicted to have \(IC_{50} < \) a chosen threshold (default: **1000 nM**)
- **Inactive** otherwise

This is intended for quick triage and prioritization, not as a replacement for confirmatory assays.

## Data Sources

### ChEMBL integration (default)
The QSAR pipeline can fetch target activity data from ChEMBL via the REST API:
- Activity filter: **IC50**
- Units normalized to **nM**
- Per-compound aggregation: **median IC50**

### Local BiochemicalResult (fallback)
If you have internal assay data, QSAR can train from the local database table:
- `BiochemicalResult.target`
- `BiochemicalResult.ic50` + `BiochemicalResult.units`
- Joined with `Compound.smiles`

## Training Models

### CLI
Use:
- `python scripts/train_qsar_models.py --targets=EGFR,CDK2 --source=chembl`
- `python scripts/train_qsar_models.py --targets=common --dry-run`

### Minimum data requirements
Training enforces:
- **≥ 100 compounds** per target
- **≥ 20% active ratio** per target (based on the IC50 threshold)

## Making Predictions

### API
`POST /api/qsar/predict`
- Request includes `smiles_list` (max 100) and `targets`
- Response returns per-compound, per-target probabilities and uncertainty

### Dashboard
Use the **Target QSAR** page:
- View available models
- Run predictions for one or more targets
- Inspect model metadata

## Interpreting Results
- **Probability**: predicted likelihood of activity (class 1).
- **Uncertainty (std)**: ensemble disagreement (higher = less stable prediction).
- **In-domain**: applicability domain check (Tanimoto similarity to training centroid).
- **Traffic light (dashboard)**:
  - **GREEN**: probability > 0.7
  - **YELLOW**: 0.3–0.7
  - **RED**: probability < 0.3

## Model Details
- **BootstrapEnsemble**: 5 XGBoost classifiers trained on bootstrap resamples.
- **Calibration**: isotonic regression on a calibration split.
- **Features**: Morgan fingerprint (2048 bits) + 6 RDKit descriptors (2054 dims).



