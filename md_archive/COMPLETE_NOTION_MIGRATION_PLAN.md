# Complete Notion Migration Plan - All Ingestion Types

## Goal
Migrate ALL ingestion pipelines away from Notion to Postgres-only operation:
- ✅ Datasets (completed)
- 🔄 Experiments (in progress)
- 📧 Emails
- 📚 Zotero/Literature

## Strategy

### Phase 1: Experiments (Postgres Model Exists) ✅
- Experiment model already in Postgres
- Create `ingest_experiment_from_postgres(experiment_id: UUID)`
- Build text from Postgres fields (name, description, type, disease, matrix, etc.)
- Link to programs/datasets from Postgres relationships
- Embed to Pinecone with Postgres metadata

### Phase 2: Emails & Zotero (New Models Needed)
- Create Postgres models for Email and Literature items
- Store email/literature metadata in Postgres
- Ingest content directly to Pinecone
- Link features to Postgres

### Phase 3: Cleanup
- Remove all Notion API calls
- Update configuration to disable Notion
- Remove unused Notion client code

## Implementation Details

### Experiments
- **Source**: Postgres Experiment model
- **Fields**: name, description, type, disease, matrix, model_systems
- **Relationships**: programs, datasets (already linked)
- **Content**: Build from description + relationships

### Emails
- **Options**:
  1. Create Email model in Postgres (store metadata)
  2. Ingest directly to Pinecone (no Postgres storage)
- **Recommendation**: Option 2 for simplicity (emails are transient)

### Zotero/Literature
- **Options**:
  1. Create Literature model in Postgres (store metadata)
  2. Ingest directly to Pinecone (no Postgres storage)
- **Recommendation**: Option 1 (literature is important for linking)

## Migration Steps

1. ✅ Dataset ingestion - DONE
2. 🔄 Experiment ingestion - NEXT
3. ⏳ Email ingestion
4. ⏳ Zotero/Literature ingestion
5. ⏳ Remove Notion dependencies
6. ⏳ Update all scripts

