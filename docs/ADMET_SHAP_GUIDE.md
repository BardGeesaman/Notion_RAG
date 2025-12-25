# ADMET SHAP Explainability Guide

## Overview

**SHAP** (SHapley Additive exPlanations) is a standard approach for explaining model predictions by attributing a portion of a prediction to each input feature. Conceptually, SHAP answers:

- “Which features pushed this prediction higher?”
- “Which features pushed this prediction lower?”

Explainability matters for ADMET predictions because the output is often used to make **risk and prioritization decisions** (e.g., flagging likely hERG liability). A numeric prediction alone is useful, but an explanation helps you:

- validate whether the model is responding to chemically plausible signals
- compare two compounds and understand why they differ
- decide when uncertainty / out-of-domain signals should override the point estimate

## Feature Types

The ADMET models use a **2054-dimensional** feature vector:

1. **Morgan Fingerprint bits (2048)**
   - Indices 0–2047
   - Names: `MorganBit_0000` … `MorganBit_2047`
   - Each bit indicates presence/absence of a substructure pattern (radius=2 in our feature pipeline).

2. **RDKit Descriptors (6)**
   - Indices 2048–2053
   - Names:
     - `MolWt`
     - `MolLogP`
     - `TPSA`
     - `NumHDonors`
     - `NumHAcceptors`
     - `NumRotatableBonds`

These names are provided by `amprenta_rag/ml/admet/features.py`.

## Interpreting the Waterfall Chart

The Explain tab displays a **waterfall chart** that decomposes a prediction into additive contributions:

- **Base value**
  - The average model output over the training data (or the explainer’s expected value).
  - Think of this as the “starting point” before adding feature contributions.

- **Positive bars (green)**
  - Features that **increase** the prediction relative to the base value.

- **Negative bars (red)**
  - Features that **decrease** the prediction relative to the base value.

- **Other bar**
  - The combined contribution of the remaining features not shown in Top K.
  - This makes the explanation complete without listing all 2054 features.

- **Final output**
  - The model output implied by the explanation:

    \( \text{output} \approx \text{base} + \sum_i \text{SHAP}_i \)

In the UI we show Top K feature contributions plus an “Other” bar so the chart remains readable.

## Traffic Light Integration

The ADMET Predictor uses a traffic-light heuristic to help with rapid triage:

- **hERG (classification probability)**
  - Green: prob < 0.3
  - Yellow: 0.3–0.7
  - Red: > 0.7

- **LogS**
  - Green: > -4
  - Yellow: -6 to -4
  - Red: < -6

- **LogP**
  - Green: 1–3
  - Yellow: 0–1 or 3–5
  - Red: < 0 or > 5

How SHAP relates:

- Use SHAP to understand *why* a prediction landed in a given traffic-light band.
- If SHAP shows that the prediction is driven by a small number of high-impact features, it can be a sign to scrutinize that compound more closely.

**In-domain vs out-of-domain**

- If `in_domain=False`, treat explanations and point estimates with extra caution.
- Out-of-domain indicates low similarity to the training centroid, which can mean:
  - the model hasn’t seen similar chemistry
  - uncertainty intervals may be widened

## Global vs Local Importance

**Local importance**

- Explains *this specific compound*.
- The Explain tab is local: it shows top features for a single SMILES.

**Global importance**

- Answers: “Which features matter most across many compounds?”
- In the MVP Global Importance tab we accumulate explanations across multiple runs and compute a simple aggregate score (mean absolute SHAP magnitude over collected explanations).

## Limitations

- **Morgan bit names are not human-readable**
  - `MorganBit_0123` indicates “some substructure pattern,” but not the pattern itself.
  - Future improvement: decode bits to substructures (e.g., RDKit bitInfo → SMARTS fragments) and display them.

- **Ensemble averaging smooths individual model quirks**
  - We average SHAP values across ensemble members for stability.
  - This can hide disagreement among models, but provides a more robust explanation overall.



