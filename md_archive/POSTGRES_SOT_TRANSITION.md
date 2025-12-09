# Phase 6: Transition to Postgres as Source of Truth

**Status**: Complete  
**Roadmap**: `context/UNIFIED_STRATEGIC_ROADMAP.md` Section 3.1, Phase 6

## Overview

Phase 6 completes the transition to Postgres as the source of truth for structured data. This phase provides the infrastructure and integration points needed for ingestion pipelines to use Postgres as primary, with Notion as an optional documentation layer.

## Strategy

### Architecture Principles

1. **Postgres is Source of Truth**: All new structured data goes to Postgres first
2. **Notion is Documentation Layer**: Optional human-readable overlay (ELN, dashboards)
3. **Pinecone uses Postgres IDs**: Metadata references Postgres entities
4. **Dual-Write During Transition**: Write to both systems until fully transitioned
5. **Backward Compatible**: Existing Notion-only pipelines continue to work

### Implementation Approach

Rather than rewriting all ingestion pipelines, Phase 6 provides:

1. **Postgres Integration Utilities** - Functions to create/update entities in Postgres
2. **Postgres-Aware Embedding** - Embedding functions that use Postgres metadata
3. **Configuration Flags** - Toggle Postgres as SoT via environment variables
4. **Dual-Write Manager** - Write to both Postgres and Notion during transition

## Components

### 1. Postgres Integration (`amprenta_rag/ingestion/postgres_integration.py`)

Provides utilities for ingestion pipelines:

- `create_or_update_dataset_in_postgres()` - Create/update datasets in Postgres
- `embed_dataset_with_postgres_metadata()` - Embed with Postgres IDs

### 2. Configuration

New environment variables:

```bash
# Enable Postgres as Source of Truth
USE_POSTGRES_AS_SOT=true

# Enable dual-write (write to both Postgres and Notion)
ENABLE_DUAL_WRITE=true
```

### 3. Dual-Write Manager (`amprenta_rag/migration/dual_write.py`)

Already created in Phase 4 - provides simultaneous writes to both systems.

## Integration into Ingestion Pipelines

### Current Flow (Notion-Only)

```
File → Notion Page → Feature Linking → Signature Scoring → Pinecone Embedding
```

### New Flow (Postgres as SoT)

```
File → Postgres Dataset → Dual-Write to Notion (optional) → Feature Linking → 
Signature Scoring → Pinecone Embedding (with Postgres IDs)
```

## Migration Path

### Step 1: Enable Postgres Integration

Set environment variables:
```bash
USE_POSTGRES_AS_SOT=true
ENABLE_DUAL_WRITE=true
```

### Step 2: Update Ingestion Pipelines

Ingestion pipelines should:

1. Create/update entity in Postgres first
2. Optionally create/update Notion page (dual-write)
3. Use Postgres ID for feature linking
4. Use Postgres metadata for Pinecone embedding

### Step 3: Update Embedding

Use `embed_dataset_with_postgres_metadata()` instead of Notion-only embedding functions.

### Step 4: Update RAG Queries

Use hybrid chunk collection (already implemented in Phase 5):
- Prefer Postgres structured data
- Fallback to Notion narrative

## Example Integration

### Before (Notion-Only)

```python
# Create Notion page
page_id = create_omics_dataset_page(...)

# Embed with Notion ID
embed_lipidomics_dataset(page_id, dataset_name, species)
```

### After (Postgres as SoT)

```python
from amprenta_rag.ingestion.postgres_integration import (
    create_or_update_dataset_in_postgres,
    embed_dataset_with_postgres_metadata,
)

# Create in Postgres
dataset = create_or_update_dataset_in_postgres(
    name=dataset_name,
    omics_type=OmicsType.LIPIDOMICS,
    file_paths=[file_path],
)

# Optional: Dual-write to Notion
if config.pipeline.enable_dual_write:
    notion_page_id = create_omics_dataset_page(...)
    dataset.notion_page_id = notion_page_id
    db.commit()

# Embed with Postgres ID
embed_dataset_with_postgres_metadata(
    dataset_id=dataset.id,
    dataset_name=dataset.name,
    species_or_features=species,
    omics_type=OmicsType.LIPIDOMICS,
    notion_page_id=dataset.notion_page_id,  # Optional
)
```

## Configuration

### Environment Variables

```bash
# Postgres as Source of Truth
USE_POSTGRES_AS_SOT=true          # Enable Postgres as primary
ENABLE_DUAL_WRITE=true             # Also write to Notion (transition)

# Postgres Connection (from Phase 2)
POSTGRES_URL=postgresql://...      # Or individual components
```

### Pipeline Configuration

The `PipelineConfig` now includes:
- `use_postgres_as_sot`: Use Postgres as source of truth
- `enable_dual_write`: Write to both Postgres and Notion

## Testing

### Test Postgres Integration

1. Set `USE_POSTGRES_AS_SOT=true`
2. Ingest a dataset
3. Verify dataset created in Postgres
4. Verify metadata uses Postgres ID
5. Verify RAG queries work with Postgres data

### Test Dual-Write

1. Set `USE_POSTGRES_AS_SOT=true` and `ENABLE_DUAL_WRITE=true`
2. Ingest a dataset
3. Verify dataset in both Postgres and Notion
4. Verify metadata includes both IDs

## Rollout Plan

1. **Phase 6.1**: Infrastructure ready (✅ Complete)
2. **Phase 6.2**: Update one pipeline as pilot (e.g., lipidomics)
3. **Phase 6.3**: Roll out to all pipelines
4. **Phase 6.4**: Disable Notion writes (Notion becomes read-only)
5. **Phase 6.5**: Final validation

## Current Status

✅ **Infrastructure Complete**:
- Postgres integration utilities created
- Postgres-aware embedding functions
- Configuration flags added
- Dual-write manager available

⏳ **Next Steps**:
- Update individual ingestion pipelines to use Postgres integration
- Test with real data ingestion
- Gradually roll out to all pipelines

## Benefits

1. **Scalability**: Postgres handles large datasets better than Notion
2. **Reliability**: Postgres is more robust for structured data
3. **Performance**: Faster queries and relationships
4. **Flexibility**: Easy to add APIs, frontends, etc.
5. **Separation**: Structured data in Postgres, narratives in Notion

## Notes

- Existing test data in Notion/Pinecone is not migrated
- New data will be ingested into Postgres
- Notion remains for ELN/documentation
- Backward compatible with existing Notion-only workflows

