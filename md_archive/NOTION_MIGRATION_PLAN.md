# Complete Notion Migration Plan

## Goal
Remove all Notion dependencies to eliminate slow API calls and enable fast, scalable Postgres-only operation.

## Current Notion Dependencies

### Critical (Blocks Postgres-only Operation)
1. **Dataset Ingestion** (`dataset_ingestion.py`)
   - Requires Notion page ID as input
   - Fetches dataset metadata from Notion
   - Updates Notion with embedding metadata

2. **Omics Ingestion Pipelines** (all omics types)
   - Optional Notion sync still makes API calls
   - Fetches Notion pages even when not needed

3. **Repository Harvest**
   - ✅ Already fixed to be Postgres-only by default
   - Still has optional Notion creation

### Non-Critical (Optional Features)
4. **Feature Linking**
   - Has Postgres linking already
   - Still has Notion linking as fallback

5. **Signature Matching**
   - Fetches signatures from Notion
   - Should fetch from Postgres instead

6. **Embedding Metadata Updates**
   - Updates Notion after Pinecone upsert
   - Can be removed or made optional

## Migration Steps

### Phase 1: Core Dataset Ingestion (Priority 1) ⚡
- [ ] Create `ingest_dataset_from_postgres(dataset_id: UUID)` function
- [ ] Extract dataset metadata from Postgres instead of Notion
- [ ] Build text content from Postgres fields (description, file paths, etc.)
- [ ] Remove Notion page fetching
- [ ] Make embedding metadata updates optional

### Phase 2: Repository Ingestion (Priority 1) ⚡
- [ ] Update harvest script to use Postgres-only ingestion
- [ ] Remove `--create-notion` requirement from ingestion
- [ ] Test end-to-end repository harvest + ingestion

### Phase 3: Omics Ingestion Cleanup (Priority 2)
- [ ] Remove Notion sync from all omics ingestion pipelines
- [ ] Make Notion completely optional (disabled by default)
- [ ] Update all ingestion scripts

### Phase 4: Feature & Signature Systems (Priority 2)
- [ ] Update signature matching to use Postgres
- [ ] Remove Notion feature linking
- [ ] Update all feature-related queries

### Phase 5: Cleanup (Priority 3)
- [ ] Remove unused Notion client imports
- [ ] Update configuration to disable Notion by default
- [ ] Remove Notion from documentation (or mark as deprecated)

## Implementation Details

### New Function: `ingest_dataset_from_postgres()`

```python
def ingest_dataset_from_postgres(
    dataset_id: UUID,
    force: bool = False,
    update_notion: bool = False  # Optional, disabled by default
) -> None:
    """
    Ingest a dataset from Postgres into Pinecone.
    
    No Notion dependencies - works entirely from Postgres data.
    """
    # 1. Fetch dataset from Postgres
    # 2. Build text content from description, file paths, external IDs
    # 3. Extract metadata from Postgres fields
    # 4. Chunk and embed
    # 5. Upsert to Pinecone
    # 6. Optionally update Notion if requested
```

### Configuration Changes

```python
# In config.py
class PipelineConfig:
    use_postgres_as_sot: bool = True
    enable_notion_sync: bool = False  # Disabled by default
    notion_optional: bool = True  # Notion is always optional
```

## Testing Plan

1. **Unit Tests**
   - Test Postgres dataset ingestion without Notion
   - Test repository harvest + ingestion end-to-end
   - Verify no Notion API calls when disabled

2. **Integration Tests**
   - Full repository study import
   - Feature linking from Postgres
   - Embedding and RAG query

3. **Performance Tests**
   - Compare ingestion speed (Notion vs Postgres-only)
   - Measure API call reduction

## Rollout Strategy

1. ✅ **Completed**: Repository harvest Postgres-only
2. **Next**: Core dataset ingestion (this session)
3. **Then**: Repository ingestion integration
4. **Finally**: Cleanup and documentation

## Breaking Changes

- `ingest_dataset(page_id)` → `ingest_dataset_from_postgres(dataset_id)`
- Notion page IDs no longer required for ingestion
- Notion sync must be explicitly enabled if desired

## Rollback Plan

- Keep Notion functions but mark as deprecated
- Maintain dual-write capability during transition
- Can re-enable Notion sync via config if needed

