# TIER 3: Architecture Evolution - Postgres + FastAPI

**Status**: In Progress  
**Roadmap**: `context/UNIFIED_STRATEGIC_ROADMAP.md` Section 3.1

## Overview

This is a **major architectural change** to migrate from Notion-as-database to Postgres-as-database, with FastAPI as the service layer. The goal is to improve scalability, performance, and enable multi-user support.

## Architecture Principles

1. **Postgres becomes system of record** for structured data
2. **SQLite remains internal detail** or stepping stone (chemistry data)
3. **Notion remains human-friendly overlay** (ELN, SOPs, dashboards)
4. **Pinecone remains semantic index** (unchanged)
5. **Keep ID semantics stable** across systems

## Implementation Phases

### ✅ Phase 1: Domain Model Extraction (COMPLETE)

- Created unified domain models in `amprenta_rag/models/domain.py`
- Models: `Program`, `Experiment`, `Dataset`, `Feature`, `Signature`
- Type-safe enums: `OmicsType`, `FeatureType`, `SignatureDirection`
- Created tracking system: `scripts/check_implementation_status.py`

### ✅ Phase 2: Postgres Schema Design (COMPLETE)

- ✅ Added Postgres configuration to `amprenta_rag/config.py`
- ✅ Created SQLAlchemy models in `amprenta_rag/database/models.py`
- ✅ Set up database base configuration in `amprenta_rag/database/base.py`
- ✅ Set up Alembic migrations (`alembic/` directory)
- ✅ Created helper scripts for migrations
- ✅ Added UUID default generation to models

### ⏳ Phase 3: FastAPI Service Layer (PENDING)

- Set up FastAPI application structure
- Implement core API endpoints:
  - Programs
  - Experiments
  - Datasets
  - Features
  - Signatures
  - Compounds
  - Screening
- Add authentication (basic, later OIDC)

### ⏳ Phase 4: Migration Utilities (PENDING)

- Export Notion + SQLite data
- Bootstrap Postgres with initial data
- Implement dual-write capability (transition phase)

### ⏳ Phase 5: RAG Integration with Postgres (PENDING)

- Update RAG builders for Postgres + Notion hybrid
- DB for structured details, Notion for narratives
- Update metadata to use Postgres IDs

### ⏳ Phase 6: Transition to Postgres SoT (PENDING)

- Pivot ingestion to Postgres as primary
- Keep Notion as documentation layer
- Final testing and validation

## Configuration

Add these environment variables to your `.env` file:

```bash
# Postgres Database (optional - for TIER 3)
POSTGRES_URL=postgresql://user:password@localhost:5432/amprenta_rag
# OR individual components:
POSTGRES_HOST=localhost
POSTGRES_PORT=5432
POSTGRES_DB=amprenta_rag
POSTGRES_USER=postgres
POSTGRES_PASSWORD=your_password
POSTGRES_ECHO=false  # Set to true for SQL logging
```

## Database Schema

The Postgres schema includes:

- `programs` - Research programs/therapeutic areas
- `experiments` - Experimental studies/assays
- `datasets` - Experimental datasets (omics data)
- `features` - Biological features (genes, proteins, metabolites, lipids)
- `signatures` - Multi-omics signature definitions
- `signature_components` - Individual components of signatures

All tables include:
- UUID primary keys
- `notion_page_id` for migration support (dual-write)
- `external_ids` JSON field for flexibility
- Timestamps (`created_at`, `updated_at`)

## Migration Strategy

1. **Dual-write phase**: Write to both Notion and Postgres
2. **Gradual migration**: Move ingestion pipelines one by one
3. **Notion remains read-only mirror** during transition
4. **Final cutover**: Postgres becomes source of truth

## Files Created

- `amprenta_rag/models/domain.py` - Domain models (dataclasses)
- `amprenta_rag/database/base.py` - Database configuration
- `amprenta_rag/database/models.py` - SQLAlchemy models
- `scripts/check_implementation_status.py` - Tracking system

## Next Steps

1. Set up Alembic migrations
2. Create initial migration script
3. Begin Phase 3 (FastAPI service layer)

## Estimated Effort

- Phase 1: ✅ Complete (3-4 days)
- Phase 2: 🚧 In Progress (5-7 days)
- Phase 3: ⏳ Pending (6-8 days)
- Phase 4: ⏳ Pending (6-8 days)
- Phase 5: ⏳ Pending (4-5 days)
- Phase 6: ⏳ Pending (1-2 days)

**Total**: 25-34 days

