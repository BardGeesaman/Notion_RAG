# CRISPR Screen Analysis (MAGeCK)

## Overview

The CRISPR screen module provides an MVP pipeline to:

1. Store screen metadata (`CRISPRScreen`) linked to a `Dataset`
2. Run MAGeCK (`mageck test`) on a count matrix
3. Parse gene-level results and persist them as `CRISPRResult`
4. Explore results in the dashboard (table + volcano plot)

## Count Matrix Format

The analysis expects a MAGeCK-style count matrix with:

- Required columns:
  - `sgRNA` (guide identifier / sequence ID)
  - `Gene` (target gene symbol)
- Sample columns:
  - One or more sample columns containing integer read counts

Example (TSV):

```text
sgRNA    Gene     control    treatment
sg1      TP53     120        80
sg2      TP53     140        90
sg3      BRCA1    200        260
```

Important:

- `CRISPRScreen.control_label` and `CRISPRScreen.treatment_label` must match sample column names exactly.
- The analysis uses `Dataset.file_paths[0]` as the path to this count matrix.

## MAGeCK Installation

The runner is docker-first with a local fallback:

1. Docker (recommended): `biocontainers/mageck`
2. Local binary fallback: `mageck` available on PATH

The runner will:

- Copy the input count matrix into the run output folder
- Run `mageck test`
- Expect `{output_prefix}.gene_summary.txt` to be created

## Confidence / Hit Calling

The MVP hit definition is:

- `is_hit = (fdr < 0.05)`

Results are ranked by ascending FDR (NULL FDR sorts last).

## API Usage Examples

### Create a screen

```bash
curl -X POST "http://localhost:8000/api/crispr/screens" \
  -H "Content-Type: application/json" \
  -d '{
    "name": "My CRISPR Screen",
    "dataset_id": "<dataset_uuid>",
    "library_type": "Brunello",
    "cell_line": "A375",
    "control_label": "control",
    "treatment_label": "treatment"
  }'
```

### Run analysis

```bash
curl -X POST "http://localhost:8000/api/crispr/screens/<screen_id>/analyze" \
  -H "Content-Type: application/json" \
  -d '{"method":"test"}'
```

### Fetch results

```bash
curl -X GET "http://localhost:8000/api/crispr/screens/<screen_id>/results?fdr_threshold=0.05&limit=100"
```

### Fetch volcano plot points

```bash
curl -X GET "http://localhost:8000/api/crispr/screens/<screen_id>/volcano-data"
```


