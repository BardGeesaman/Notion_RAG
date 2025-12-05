# Postgres Migration Guide

**Last Updated**: 2025-12-04  
**Status**: Postgres Integration Complete

## Overview

All omics ingestion pipelines (Lipidomics, Metabolomics, Proteomics, Transcriptomics) now support Postgres as the primary database. This guide explains the migration, architecture, and how to use the new system.

## Architecture

### Previous Architecture (Notion-Only)
```
File → Notion Page → Pinecone Embedding
```

### New Architecture (Postgres-First)
```
File → Postgres Dataset → Optional Notion Sync → Pinecone Embedding
```

## Key Features

### 1. Postgres as Single Source of Truth
- All datasets stored in Postgres first
- Fast bulk ingestion (no Notion rate limits)
- Direct database access for queries
- Relational integrity

### 2. Optional Notion Sync
- Notion can be used as documentation layer
- Controlled by `ENABLE_NOTION_SYNC` config
- Links Postgres datasets to Notion pages

### 3. Postgres-Aware Embeddings
- Embeddings use Postgres metadata
- Fallback to Notion if Postgres unavailable
- Consistent metadata across all omics types

## Configuration

### Required Settings

Add to your `.env` file:

```bash
# Postgres Database
POSTGRES_HOST=localhost
POSTGRES_PORT=5432
POSTGRES_DB=amprenta_rag
POSTGRES_USER=your_user
POSTGRES_PASSWORD=your_password

# Enable Postgres as source of truth
USE_POSTGRES_AS_SOT=true

# Optional: Enable Notion sync (for documentation)
ENABLE_NOTION_SYNC=false
ENABLE_DUAL_WRITE=false
```

### Configuration Options

- `USE_POSTGRES_AS_SOT`: Use Postgres as primary database (default: `true`)
- `ENABLE_NOTION_SYNC`: Sync to Notion for documentation (default: `false`)
- `ENABLE_DUAL_WRITE`: Write to both Postgres and Notion during transition (default: `false`)

## Usage

### Ingestion Examples

#### Metabolomics
```python
from amprenta_rag.ingestion.metabolomics.ingestion import ingest_metabolomics_file

# Postgres-only ingestion
dataset_id = ingest_metabolomics_file(
    file_path="data/metabolomics_sample.csv",
    create_page=False,  # No Notion page needed
)
```

#### Proteomics
```python
from amprenta_rag.ingestion.proteomics.ingestion import ingest_proteomics_file

dataset_id = ingest_proteomics_file(
    file_path="data/proteomics_sample.csv",
    create_page=False,
)
```

#### Transcriptomics
```python
from amprenta_rag.ingestion.transcriptomics.ingestion import ingest_transcriptomics_file

dataset_id = ingest_transcriptomics_file(
    file_path="data/transcriptomics_sample.csv",
    create_page=False,
)
```

#### Lipidomics
```python
from amprenta_rag.ingestion.lipidomics.ingestion import ingest_lipidomics_file

dataset_id = ingest_lipidomics_file(
    file_path="data/lipidomics_sample.csv",
    create_page=False,
)
```

### With Notion Sync

If you want to also create Notion pages for documentation:

```python
# Set in .env: ENABLE_NOTION_SYNC=true

dataset_id = ingest_metabolomics_file(
    file_path="data/metabolomics_sample.csv",
    create_page=True,  # Creates Notion page and links to Postgres
)
```

### Linking to Existing Notion Pages

```python
dataset_id = ingest_metabolomics_file(
    file_path="data/metabolomics_sample.csv",
    notion_page_id="existing-notion-page-id",
    create_page=False,
)
```

## Migration Path

### Option 1: Fresh Start (Recommended)
- Start using Postgres-only ingestion
- No migration needed
- Faster, cleaner architecture

### Option 2: Gradual Migration
- Enable `ENABLE_DUAL_WRITE=true`
- Both Postgres and Notion updated
- Switch to Postgres-only when ready

### Option 3: Notion-First
- Keep using Notion as primary
- Set `USE_POSTGRES_AS_SOT=false`
- Postgres integration disabled

## Benefits

### Performance
- **10-100x faster** bulk ingestion
- No API rate limits
- Direct database queries

### Scalability
- Handle thousands of datasets
- Efficient bulk operations
- Better query performance

### Consistency
- Unified architecture across omics
- Relational integrity
- Better data quality

## API Access

### FastAPI Endpoints

All datasets accessible via REST API:

```bash
# List all datasets
curl http://localhost:8000/api/v1/datasets

# Get specific dataset
curl http://localhost:8000/api/v1/datasets/{dataset_id}

# Filter by omics type
curl http://localhost:8000/api/v1/datasets?omics_type=Metabolomics
```

### Dashboard

Access the Streamlit dashboard:
```bash
streamlit run scripts/run_dashboard.py
```

Browse datasets at: http://localhost:8501

## Troubleshooting

### Postgres Connection Issues

**Error**: `Postgres creation failed (required)`

**Solution**:
1. Verify Postgres is running
2. Check connection credentials in `.env`
3. Ensure database exists
4. Verify user permissions

### Notion Sync Issues

**Error**: `Notion sync skipped (error)`

**Solution**:
- Notion sync is optional
- If `USE_POSTGRES_AS_SOT=true`, Notion failures won't block ingestion
- Check `ENABLE_NOTION_SYNC` setting if you need Notion pages

### Feature Linking Issues

Feature linking currently uses Notion:
- Requires `ENABLE_NOTION_SYNC=true` or `ENABLE_FEATURE_LINKING=false`
- Future: Feature linking will use Postgres directly

## Next Steps

1. **Feature Linking in Postgres**: Link features to datasets in Postgres
2. **Program/Experiment Linking**: Link datasets to programs/experiments
3. **Enhanced Dashboard**: More visualizations and filters
4. **Performance Optimization**: Database indexing, caching

See `docs/NEXT_STEPS.md` for full roadmap.

## Support

- Configuration: `amprenta_rag/config.py`
- Postgres Integration: `amprenta_rag/ingestion/postgres_integration.py`
- Dashboard: `scripts/run_dashboard.py`
- API: `amprenta_rag/api/`

