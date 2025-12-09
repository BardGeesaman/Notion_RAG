# Priority 4 Progress: Program/Experiment Linking

**Date**: 2025-12-04  
**Status**: 🚧 In Progress

## Summary

Implementing Postgres-based linking for datasets to programs and experiments, completing the relationship tracking infrastructure.

## Tasks Completed

### ✅ 1. Created Postgres Program/Experiment Linking Module

**File**: `amprenta_rag/ingestion/postgres_program_experiment_linking.py`

**Functions Created**:
- `find_program_by_notion_id()` - Find Program by Notion page ID
- `find_experiment_by_notion_id()` - Find Experiment by Notion page ID
- `convert_notion_ids_to_postgres_uuids()` - Convert Notion IDs to Postgres UUIDs
- `link_dataset_to_programs_in_postgres()` - Link dataset to programs
- `link_dataset_to_experiments_in_postgres()` - Link dataset to experiments
- `link_dataset_to_programs_and_experiments_in_postgres()` - Combined linking function

**Features**:
- ✅ Automatic Notion ID to Postgres UUID conversion
- ✅ Idempotent linking (checks for existing links)
- ✅ Transaction-safe operations
- ✅ Comprehensive error handling and logging
- ✅ Supports both direct UUIDs and Notion page IDs

## Architecture

### ID Mapping Strategy

The system supports two approaches:

1. **Direct UUIDs**: If you already have Postgres UUIDs
2. **Notion Page IDs**: Automatically converted to Postgres UUIDs

```python
# Option 1: Direct UUIDs
link_dataset_to_programs_and_experiments_in_postgres(
    dataset_id=dataset_uuid,
    program_ids=[program_uuid1, program_uuid2],
    experiment_ids=[experiment_uuid1],
)

# Option 2: Notion Page IDs (auto-converted)
link_dataset_to_programs_and_experiments_in_postgres(
    dataset_id=dataset_uuid,
    notion_program_ids=["notion-page-id-1", "notion-page-id-2"],
    notion_experiment_ids=["notion-page-id-3"],
)
```

### Database Schema

The linking uses existing association tables:
- `program_dataset_assoc` - Links programs to datasets
- `experiment_dataset_assoc` - Links experiments to datasets

Both tables support many-to-many relationships and are already defined in the schema.

## Next Steps

### Remaining Tasks

1. **Integrate into Ingestion Pipelines** 🚧
   - Update all omics ingestion modules
   - Replace TODO comments with actual Postgres linking
   - Support both Notion IDs (for backward compatibility) and Postgres UUIDs

2. **Update Harvest Script** 🚧
   - Add program/experiment linking for repository studies
   - Support metadata-based program/experiment detection

3. **Update Dashboard** 🚧
   - Show Postgres-linked programs/experiments
   - Display relationship counts
   - Add navigation between related entities

4. **Testing** 🚧
   - Test with real datasets
   - Verify linking works correctly
   - Test Notion ID conversion

## Integration Points

### Ingestion Pipelines

Current state:
- All pipelines have TODO comments: `# TODO: Link to Programs/Experiments in Postgres`
- Currently only link to Notion

Planned update:
```python
# After Postgres dataset creation
if cfg.pipeline.use_postgres_as_sot and postgres_dataset:
    from amprenta_rag.ingestion.postgres_program_experiment_linking import (
        link_dataset_to_programs_and_experiments_in_postgres,
    )
    
    link_dataset_to_programs_and_experiments_in_postgres(
        dataset_id=postgres_dataset.id,
        notion_program_ids=program_ids,  # Auto-converted
        notion_experiment_ids=experiment_ids,  # Auto-converted
    )
```

### Harvest Script

Current state:
- Creates Postgres datasets
- No program/experiment linking yet

Planned update:
- Extract program/experiment information from repository metadata
- Link to existing programs/experiments or create new ones

## Testing Checklist

- [ ] Test finding programs by Notion ID
- [ ] Test finding experiments by Notion ID
- [ ] Test ID conversion (Notion → Postgres)
- [ ] Test linking datasets to programs
- [ ] Test linking datasets to experiments
- [ ] Test idempotency (linking twice doesn't create duplicates)
- [ ] Test with real datasets from ingestion pipelines
- [ ] Test with repository-harvested datasets
- [ ] Verify dashboard displays relationships correctly

## Files

**Created**:
- `amprenta_rag/ingestion/postgres_program_experiment_linking.py`
- `docs/PRIORITY_4_PROGRESS.md` (this file)

**To Update**:
- `amprenta_rag/ingestion/metabolomics/ingestion.py`
- `amprenta_rag/ingestion/proteomics/ingestion.py`
- `amprenta_rag/ingestion/transcriptomics/ingestion.py`
- `amprenta_rag/ingestion/lipidomics/ingestion.py`
- `scripts/harvest_repository_study.py`
- `scripts/run_dashboard.py`

## See Also

- `docs/NEXT_STEPS.md` - Roadmap and priorities
- `docs/POSTGRES_MIGRATION_GUIDE.md` - Postgres architecture details

