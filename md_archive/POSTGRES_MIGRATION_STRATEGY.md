# Postgres Migration Strategy - Performance Optimization

**Status**: Active Migration  
**Goal**: Make Postgres primary, Notion optional documentation layer  
**Reason**: Improve performance for bulk data ingestion

## Current State

- ✅ Postgres infrastructure built and tested
- ✅ Postgres integration added to lipidomics ingestion (pilot)
- ⚠️ Notion still primary for most operations
- ⚠️ Bulk ingestion performance limited by Notion API rate limits

## Migration Goals

1. **Postgres becomes primary** for all structured data
2. **Notion becomes optional** documentation layer (can be disabled)
3. **Performance**: 10-100x faster bulk ingestion
4. **Scalability**: Handle thousands of datasets without rate limits
5. **Backward compatibility**: Existing Notion data remains accessible

## Migration Phases

### Phase 1: Make Postgres Primary (Current)

**Goal**: Update ingestion pipelines to use Postgres first, Notion optionally

**Changes Needed**:
1. ✅ Update `USE_POSTGRES_AS_SOT` default to `true`
2. ✅ Update all ingestion pipelines to create Postgres records first
3. ✅ Make Notion creation optional (controlled by flag)
4. ✅ Update RAG queries to use Postgres as primary source

**Status**: In Progress

### Phase 2: Update All Ingestion Pipelines

**Pipelines to Update**:
- ✅ Lipidomics (done)
- ⏳ Metabolomics
- ⏳ Proteomics
- ⏳ Transcriptomics
- ⏳ Signatures
- ⏳ Features
- ⏳ Programs/Experiments

### Phase 3: Update RAG Queries

**Changes**:
- Use Postgres for structured metadata
- Use Notion only for narrative text (if available)
- Fallback to Postgres if Notion unavailable

### Phase 4: Make Notion Optional

**Changes**:
- Add `ENABLE_NOTION_SYNC` flag (default: false)
- Only create/update Notion if flag enabled
- Notion becomes read-only documentation layer

## Configuration

### New Configuration Flags

```bash
# .env file

# Postgres as Source of Truth (REQUIRED for performance)
USE_POSTGRES_AS_SOT=true

# Optional: Sync to Notion (for documentation)
ENABLE_NOTION_SYNC=false  # Set to true if you want Notion pages

# Dual-write mode (transition period only)
ENABLE_DUAL_WRITE=false  # Only needed during migration
```

### Performance Impact

**Before (Notion Primary)**:
- ~1-2 seconds per dataset (Notion API rate limits)
- Bulk ingestion: ~30-60 datasets/minute
- Rate limit errors common

**After (Postgres Primary)**:
- ~10-50ms per dataset (Postgres local)
- Bulk ingestion: ~1000+ datasets/minute
- No rate limits

## Implementation Plan

### Step 1: Update Default Configuration

Change default `USE_POSTGRES_AS_SOT` to `true` in `config.py`

### Step 2: Update Ingestion Pipelines

For each pipeline:
1. Create Postgres record first
2. Optionally create Notion page (if `ENABLE_NOTION_SYNC=true`)
3. Link Postgres and Notion via `notion_page_id`

### Step 3: Update RAG Builders

1. Use Postgres as primary metadata source
2. Fallback to Notion only for narrative text
3. Update Pinecone metadata to use Postgres IDs

### Step 4: Update Query Functions

1. Query Postgres first
2. Optionally enrich with Notion narrative
3. Return Postgres-structured data

## Migration Checklist

- [x] Postgres infrastructure built
- [x] Lipidomics ingestion updated (pilot)
- [ ] Update default config to use Postgres
- [ ] Update all ingestion pipelines
- [ ] Update RAG builders
- [ ] Update query functions
- [ ] Test bulk ingestion performance
- [ ] Document new workflow
- [ ] Update CLI scripts

## Rollback Plan

If issues arise:
1. Set `USE_POSTGRES_AS_SOT=false` in `.env`
2. System reverts to Notion-primary mode
3. Postgres data remains (can be re-synced later)

## Testing

### Performance Test

```bash
# Test bulk ingestion
python scripts/batch_ingest_lipidomics.py --dir /path/to/many/files

# Measure time
time python scripts/batch_ingest_lipidomics.py --dir /path/to/many/files
```

### Functional Test

```bash
# Ingest single dataset
python scripts/ingest_lipidomics.py --file test.csv --create-page

# Verify in Postgres
psql amprenta_rag -c "SELECT name, omics_type FROM datasets;"

# Verify in Notion (if enabled)
# Check Experimental Data Assets database
```

## Next Steps

1. **Update default config** - Make Postgres primary by default
2. **Update all pipelines** - Migrate remaining ingestion modules
3. **Update RAG** - Use Postgres metadata
4. **Test performance** - Verify 10-100x improvement
5. **Document** - Update user guides

