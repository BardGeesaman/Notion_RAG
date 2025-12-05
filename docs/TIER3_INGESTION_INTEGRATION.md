# TIER 3 Postgres Integration into Ingestion Pipelines

**Status**: Pilot Integration Complete  
**Last Updated**: 2025-01-XX

## Overview

Postgres integration has been added to the lipidomics ingestion pipeline as a pilot. This allows datasets to be created in Postgres alongside (or instead of) Notion, enabling the transition to Postgres as source of truth.

## Integration Points

### 1. Dataset Creation

The `ingest_lipidomics_file()` function now:

1. **Creates Notion page** (if enabled/required)
2. **Creates Postgres dataset** (if `USE_POSTGRES_AS_SOT=true`)
3. **Links both** via `notion_page_id` field in Postgres

### 2. Embedding

The embedding process now:

1. **Uses Postgres metadata** (if Postgres dataset exists)
2. **Falls back to Notion** (if Postgres not enabled or fails)
3. **Maintains backward compatibility**

## Configuration

### Enable Postgres Integration

Add to `.env` file:

```bash
# Enable Postgres as Source of Truth
USE_POSTGRES_AS_SOT=true

# Enable dual-write (write to both Postgres and Notion)
ENABLE_DUAL_WRITE=true
```

### Current Behavior

- **`USE_POSTGRES_AS_SOT=false`** (default):
  - Creates dataset in Notion only
  - Uses Notion-based embedding
  - Backward compatible with existing workflow

- **`USE_POSTGRES_AS_SOT=true`**:
  - Creates dataset in Postgres
  - Optionally creates in Notion (if `ENABLE_DUAL_WRITE=true`)
  - Uses Postgres metadata for embedding
  - Links Postgres and Notion via `notion_page_id`

## Integration Flow

### Current Flow (Notion-Only)

```
File → Extract Species → Notion Page → Feature Linking → 
Signature Scoring → Notion Embedding → Pinecone
```

### New Flow (Postgres as SoT)

```
File → Extract Species → Postgres Dataset → (Optional: Notion Page) →
Feature Linking → Signature Scoring → Postgres Embedding → Pinecone
```

## Implementation Details

### Dataset Creation

```python
# In ingest_lipidomics_file()
if cfg.pipeline.use_postgres_as_sot:
    postgres_dataset = create_or_update_dataset_in_postgres(
        name=dataset_name,
        omics_type=OmicsType.LIPIDOMICS,
        file_paths=[file_path],
        description=...,
        notion_page_id=page_id,  # Links to Notion
    )
```

### Embedding

```python
# Uses Postgres metadata if available
if cfg.pipeline.use_postgres_as_sot and postgres_dataset:
    embed_dataset_with_postgres_metadata(
        dataset_id=postgres_dataset.id,
        dataset_name=dataset_name,
        species_or_features=list(species_set),
        omics_type=OmicsType.LIPIDOMICS,
        signature_matches=signature_matches,
        notion_page_id=page_id,
    )
```

## Testing

### Test Postgres Integration

1. **Enable Postgres**:
   ```bash
   # In .env
   USE_POSTGRES_AS_SOT=true
   ENABLE_DUAL_WRITE=true
   ```

2. **Ingest a dataset**:
   ```bash
   python scripts/ingest_lipidomics.py --file path/to/file.csv --create-page
   ```

3. **Verify in Postgres**:
   ```bash
   python scripts/validate_postgres_setup.py
   # Or query directly
   psql amprenta_rag -c "SELECT name, omics_type FROM datasets;"
   ```

4. **Verify in Notion**:
   - Check Experimental Data Assets database
   - Verify page exists and is linked

5. **Verify in Pinecone**:
   - Check metadata includes Postgres ID
   - Verify chunks are embedded

## Rollout Plan

### Phase 1: Pilot (Current)
- ✅ Lipidomics ingestion integrated
- ✅ Postgres creation working
- ✅ Embedding with Postgres metadata
- ✅ Backward compatible

### Phase 2: Expand to Other Omics
- ⏳ Metabolomics ingestion
- ⏳ Proteomics ingestion
- ⏳ Transcriptomics ingestion

### Phase 3: Full Migration
- ⏳ All pipelines use Postgres
- ⏳ Notion becomes read-only documentation layer
- ⏳ RAG queries use Postgres as primary

## Benefits

1. **Scalability**: Postgres handles large datasets better
2. **Performance**: Faster queries and relationships
3. **Flexibility**: Easy to add APIs, frontends, etc.
4. **Separation**: Structured data in Postgres, narratives in Notion
5. **Backward Compatible**: Existing workflows continue to work

## Notes

- Integration is **non-blocking**: If Postgres fails, Notion workflow continues
- **Idempotent**: Re-running ingestion updates existing datasets
- **Dual-write**: Can write to both systems during transition
- **Graceful degradation**: Falls back to Notion if Postgres unavailable

## Next Steps

1. **Test with real data**: Ingest actual lipidomics datasets
2. **Expand to other omics**: Add Postgres integration to other pipelines
3. **Update RAG queries**: Use Postgres as primary source
4. **Monitor performance**: Track ingestion times and success rates

