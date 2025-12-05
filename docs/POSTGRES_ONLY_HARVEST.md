# Postgres-Only Repository Harvest

## Overview

The repository harvest script now supports **Postgres-only** operation by default. Notion is completely optional and only used when explicitly requested with the `--create-notion` flag.

## Migration Status

✅ **Postgres is the default** - Datasets are created in Postgres by default  
✅ **Notion is optional** - Only used if `--create-notion` is passed  
⚠️ **Ingestion** - Currently requires Notion page ID (Postgres-only ingestion coming soon)

## Usage

### Postgres-Only (Recommended)

```bash
# Create dataset in Postgres only (no Notion queries)
python scripts/harvest_repository_study.py \
    --study-id ST004168 \
    --repository MW

# With ingestion (requires --create-notion for now)
python scripts/harvest_repository_study.py \
    --study-id ST004168 \
    --repository MW \
    --ingest \
    --create-notion
```

### With Notion Sync (Optional)

```bash
# Create dataset in Postgres AND Notion
python scripts/harvest_repository_study.py \
    --study-id ST004168 \
    --repository MW \
    --create-notion
```

## Why Notion is Still Accessed

The script will query Notion **only** if:
1. You pass `--create-notion` flag (explicitly requesting Notion)
2. The script needs to check for existing Notion pages before creating new ones

## To Avoid Notion Completely

**Don't pass `--create-notion`** - The script will work entirely with Postgres:

```bash
# Postgres-only - no Notion queries
python scripts/harvest_repository_study.py \
    --study-id ST004168 \
    --repository MW
```

## Current Limitations

1. **Ingestion requires Notion**: The `--ingest` flag currently requires a Notion page ID. Postgres-only ingestion is planned.
2. **Feature linking**: Features are linked to Postgres, but embedding/RAG still uses Notion-based pipeline.

## Next Steps

- [ ] Implement Postgres-only ingestion pipeline
- [ ] Remove Notion dependency from embedding/RAG pipeline
- [ ] Add Postgres dataset search by external IDs

