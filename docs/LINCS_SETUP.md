## Overview

This repository supports **LINCS / Connectivity Map (CMap)** workflows:

- Start from an internal **Signature** (`signatures`).
- Ingest LINCS Level 5 signatures into **LINCSSignature** (`lincs_signatures`).
- Compute **ConnectivityScore** (`connectivity_scores`) between an internal signature and each LINCS signature.
  - **Reversals**: negative scores (candidate perturbagens that may reverse the signature)
  - **Mimics**: positive scores (perturbagens that may mimic the signature)

The Dashboard page **Connectivity Map** surfaces the compute + top reversals/mimics views.

---

## Prerequisites

### 1) Entrez IDs must exist on gene Features

Connectivity scoring expects internal gene features to have an Entrez ID in:

- `Feature.external_ids.entrez_id` (integer)

If this field is missing, scoring will fail with:

> No gene Entrez IDs available for signature (missing Feature.external_ids.entrez_id)

### 2) Migrations applied

Ensure Alembic migrations are applied through the LINCS schema migration:

- `alembic/versions/d4e6f7a8b9c0_add_lincs.py`

---

## Populating Entrez IDs (gene Features)

The helper module `[amprenta_rag/connectivity/gene_mapper.py](../amprenta_rag/connectivity/gene_mapper.py)` downloads and caches NCBI `gene_info.gz` and can map symbols → Entrez IDs.

### Example (Python snippet)

```python
from amprenta_rag.connectivity.gene_mapper import map_symbols_to_entrez
from amprenta_rag.database.models import Feature
from amprenta_rag.database.session import db_session

symbols = ["TP53", "EGFR", "MAPK1"]
mapping = map_symbols_to_entrez(symbols)  # {"TP53": 7157, ...}

with db_session() as db:
    genes = (
        db.query(Feature)
        .filter(Feature.feature_type == "gene")
        .filter(Feature.name.in_(symbols))
        .all()
    )
    for g in genes:
        entrez = mapping.get(g.name)
        if not entrez:
            continue
        ext = dict(g.external_ids or {})
        ext["entrez_id"] = int(entrez)
        g.external_ids = ext
        db.add(g)
```

---

## Environment variables

### LINCS download / storage

- `LINCS_LEVEL5_URL`: URL to the full Level 5 dataset (GCT or GCTX).
- `LINCS_LEVEL5_SUBSET_URL`: URL to a small subset dataset (recommended for dev/CI).
  - Precedence: if `LINCS_LEVEL5_URL` is set, it is used; otherwise the subset URL is used.
- `LINCS_STORE_BASE`: where downloaded datasets should be stored (default: `data/lincs/`).

### Cache

- `LINCS_CACHE_DIR`: cache directory for NCBI gene info and related artifacts (default: `data/lincs/`).

### Optional metadata sidecar for small test datasets

If you ingest a small `.gct` file and want to preserve `pert_iname`, `pert_id`, and `cell_id`, the parser supports an optional sidecar JSON:

- `{file}.meta.json`
- Format: `{ "sig_id": {"pert_iname": "...", "pert_id": "...", "cell_id": "..."}, ... }`

---

## Recommended test dataset strategy

Do not hardcode protected URLs in code or docs. For development, prefer setting:

- `LINCS_LEVEL5_SUBSET_URL` to a small `.gct` or `.gctx` you can access

Optionally, store a corresponding `{file}.meta.json` next to it to include perturbation metadata.

---

## Usage (API)

### 1) Ingest LINCS data

Trigger ingestion (downloads + parses + inserts LINCS signatures):

- `POST /api/connectivity/ingest`

Body:

```json
{"max_compounds": 1000}
```

### 2) Compute connectivity scores for an internal signature

- `POST /api/connectivity/compute/{signature_id}`

### 3) Read results

- Top reversals (negative scores): `GET /api/connectivity/signatures/{signature_id}/top-reversals?n=50`
- Top mimics (positive scores): `GET /api/connectivity/signatures/{signature_id}/top-mimics?n=50`

### 4) Dashboard

Open the Dashboard page: **Connectivity Map**


