# Migration Strategy - TIER 3 Architecture Evolution

## Overview

The migration from Notion-as-database to Postgres-as-database is a **forward-looking transition**, not a migration of existing test data. The strategy is:

1. **New data going forward** is ingested into both Notion and Postgres (dual-write)
2. **Existing test data in Notion/Pinecone is not migrated** - it will be reingested as needed
3. **Postgres becomes the source of truth** for all new structured data
4. **Notion remains a documentation/knowledge layer** (ELN, SOPs, dashboards)

## Dual-Write Period

During the transition, all new ingestion pipelines should use `DualWriteManager` to write to both systems simultaneously:

```python
from amprenta_rag.migration.dual_write import DualWriteManager
from amprenta_rag.database.base import get_db

# In your ingestion pipeline
db = next(get_db())
dual_write = DualWriteManager(db, enable_notion=True, enable_postgres=True)

# Create a program
program_id = dual_write.create_program({
    "name": "ALS Research Program",
    "description": "...",
    "disease": ["ALS"]
})
```

## Migration Phases

### Phase 1-3: Foundation ✅
- Domain models
- Postgres schema
- FastAPI service layer

### Phase 4: Dual-Write Capability ✅
- `DualWriteManager` for simultaneous writes
- Transition utilities

### Phase 5: RAG Integration
- Update RAG builders to use Postgres IDs
- Hybrid Postgres + Notion metadata

### Phase 6: Postgres as Source of Truth
- Switch ingestion pipelines to Postgres-first
- Notion becomes read-only documentation layer

## Key Principles

1. **No data export/migration** - Existing test data stays in Notion
2. **Re-ingest as needed** - Real data will be reingested into Postgres
3. **Dual-write during transition** - New data goes to both systems
4. **Postgres is authoritative** - Once stable, Postgres is source of truth
5. **Notion for humans** - Remains for ELN, dashboards, documentation

## What Gets Reingested

- Real experimental datasets
- Production signatures
- Actual programs and experiments
- Real chemistry/screening data

These will be ingested through the normal pipelines with dual-write enabled.

## What Stays in Notion

- Test data (not migrated)
- ELN entries (documentation)
- Dashboards and views
- Human-readable documentation

## Implementation Status

- ✅ Dual-write manager created
- ⏳ Integration into ingestion pipelines (Phase 6)
- ⏳ RAG builder updates (Phase 5)

