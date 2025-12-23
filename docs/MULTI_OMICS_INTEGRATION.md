# Multi-Omics Integration (Latent Factors / MOFA-style)

## Overview

This module supports multi-omics latent factor analysis (MOFA-style) by:

1. Creating a `MultiOmicsExperiment` that references multiple datasets (`dataset_ids`) and a cross-view `sample_mapping`
2. Extracting per-dataset feature matrices from `Dataset.file_paths[0]`
3. Aligning sample columns across views using `sample_mapping`
4. Running MOFA via `mofapy2` and saving an HDF5 model
5. Parsing the HDF5 and persisting:
   - `LatentFactor` (per factor)
   - `FactorLoading` (per feature per factor per view)
   - `FactorScore` (per sample per factor)

## Sample Mapping Format (row-oriented JSON)

The recommended `sample_mapping` is a list of rows, one per unified sample:

```json
[
  {"sample_id": "S1", "transcriptomics": "RNA_A", "proteomics": "PROT_A"},
  {"sample_id": "S2", "transcriptomics": "RNA_B", "proteomics": "PROT_B"}
]
```

- `sample_id`: unified sample identifier used across all views
- Each additional key is a view/omics name (must match the view keys used in `dataset_ids`)
- Values are the source sample column names in that view’s matrix file

## Dataset Requirements

Each dataset referenced by the experiment must have:

- `Dataset.file_paths[0]` set to a readable CSV/TSV file on disk
- A **feature identifier column** (first column by default, or one of: `feature`, `feature_id`, `gene`, `gene_id`, `id`, `name`)
- One or more **numeric sample columns** (counts/abundances/etc.)

Example CSV:

```text
feature,SampleA,SampleB
TP53,10.2,11.1
BRCA1,5.0,4.8
```

Internally, the extractor returns a DataFrame shaped **features × samples**.

## MOFA2 Installation

To run MOFA locally, install:

```bash
pip install mofapy2 h5py
```

Notes:

- The runner imports `mofapy2` inside `run_mofa()` and raises a clear ImportError if it’s not installed.
- The result parser requires `h5py` to read the output HDF5.

## API Usage Examples

### Create an experiment

```bash
curl -X POST "http://localhost:8000/api/multi-omics/experiments" \
  -H "Content-Type: application/json" \
  -d '{
    "name": "MOFA Example",
    "dataset_ids": [
      {"omics_type":"transcriptomics", "dataset_id":"<uuid>"},
      {"omics_type":"proteomics", "dataset_id":"<uuid>"}
    ],
    "sample_mapping": [
      {"sample_id":"S1","transcriptomics":"RNA_A","proteomics":"PROT_A"},
      {"sample_id":"S2","transcriptomics":"RNA_B","proteomics":"PROT_B"}
    ],
    "n_factors": 10,
    "convergence_mode": "fast"
  }'
```

### Run an experiment

```bash
curl -X POST "http://localhost:8000/api/multi-omics/experiments/<experiment_id>/run"
```

### List factors (variance per factor, if available)

```bash
curl -X GET "http://localhost:8000/api/multi-omics/experiments/<experiment_id>/factors"
```

### Fetch top loadings for a factor

```bash
curl -X GET "http://localhost:8000/api/multi-omics/experiments/<experiment_id>/factors/0/loadings?limit=50"
```

### Fetch scores for a factor

```bash
curl -X GET "http://localhost:8000/api/multi-omics/experiments/<experiment_id>/factors/0/scores"
```


