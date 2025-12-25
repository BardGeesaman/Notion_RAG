## ADMET Uncertainty Quantification (UQ) Validation Report

Status: **Placeholder / TODO populate after ChEMBL training is run** (Phase 2 training deferred due to 4GB download).

### Section 1: Dataset Summary

**Curated from**: ChEMBL SQLite (`chembl_33.db`)  
**Splits**: 70% train / 15% calibration / 15% test  

| Endpoint | Task | Train | Cal | Test | Notes |
|---|---|---:|---:|---:|---|
| hERG | Classification | TODO | TODO | TODO | Label: IC50 < 10uM |
| LogS | Regression | TODO | TODO | TODO | Units normalized to uM where possible |
| LogP | Regression | TODO | TODO | TODO | LogP/LogD values |

**hERG class balance (train)**: TODO: {0: N, 1: N}

### Section 2: Model Performance

Metrics are computed on the **test** split using the ensemble mean and 90% CI.

| Endpoint | AUC | RMSE | R2 | ECE | Brier | Coverage_90 |
|---|---:|---:|---:|---:|---:|---:|
| hERG | TODO | - | - | TODO | TODO | TODO |
| LogS | - | TODO | TODO | - | - | TODO |
| LogP | - | TODO | TODO | - | - | TODO |

**Comparison**: raw ensemble vs calibrated ensemble (classification only)

- hERG raw: AUC=TODO, ECE=TODO, Brier=TODO
- hERG calibrated: AUC=TODO, ECE=TODO, Brier=TODO

### Section 3: Ablation Study

#### 3.1 Calibration Ablation (hERG)

| Variant | AUC | ECE | Brier | Coverage_90 |
|---|---:|---:|---:|---:|
| Ensemble (raw) | TODO | TODO | TODO | TODO |
| Ensemble + isotonic | TODO | TODO | TODO | TODO |
| Ensemble + platt | TODO | TODO | TODO | TODO |

**Conclusion (expected)**: calibration improves **ECE** and often improves **Brier**, while AUC should remain similar.

#### 3.2 Applicability Domain Ablation

| Variant | Endpoint | Coverage_90 | Notes |
|---|---|---:|---|
| No domain widening | hERG/LogS/LogP | TODO | baseline CI width |
| Domain-based CI widening | hERG/LogS/LogP | TODO | OOD samples widened 2× |

**Conclusion (expected)**: widening increases Coverage_90 for OOD samples at the cost of wider intervals.

### Section 4: Uncertainty Quality

- **90% CI coverage**: TODO
- **Reliability diagrams**:
  - View in dashboard: **ADMET Predictor → Calibration**
  - TODO: save static images after models trained
- **OOD detection accuracy**:
  - Similarity threshold: 0.3 (generalized Tanimoto to training centroid)
  - TODO: quantify false positives/false negatives after training

### Section 5: Limitations

- **Label noise**: ChEMBL activities have heterogeneous assay protocols; filtering reduces but does not remove noise.
- **Endpoint definition**: hERG derived from KCNH2 binding IC50 only; other assay types may differ.
- **Applicability domain**: centroid similarity is coarse; future work may use kNN distance or conformal prediction.
- **Known failure modes**:
  - RDKit featurization failures for malformed SMILES
  - Sparse training regimes for rare chemotypes



