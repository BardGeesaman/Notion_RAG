# Notion Removal Implementation Instructions

**For:** Implementor Agent  
**From:** Architect Agent  
**Date:** 2025-01-XX  
**Status:** Ready for Implementation

## Overview

This document provides step-by-step instructions for removing all Notion dependencies from the codebase. Postgres is now the Source of Truth (100% migration complete), and all Notion code should be removed.

## Scope

- **Files with Notion references:** 114 files found
- **Total Notion imports/references:** 367 matches across 78 files
- **Estimated effort:** 4-5 weeks of systematic work

## Implementation Phases

### Phase 1: Delete Notion-Specific Modules

**Objective:** Remove entire modules that are Notion-only and have no Postgres equivalents.

**Files to Delete:**
1. `amprenta_rag/clients/notion_client.py` - Notion API client
2. `amprenta_rag/ingestion/notion_pages.py` - Notion page creation/management
3. `amprenta_rag/ingestion/dataset_notion_utils.py` - Dataset Notion operations
4. `amprenta_rag/ingestion/signature_notion_crud.py` - Signature Notion CRUD (already re-exports, safe to delete)
5. `amprenta_rag/analysis/pathway_notion_integration.py` - Pathway Notion integration
6. `amprenta_rag/chemistry/notion_integration.py` - Chemistry Notion integration
7. `amprenta_rag/migration/dual_write.py` - Dual-write manager (no longer needed)

**Before deletion:** Ensure no critical functionality depends on these modules. All have Postgres equivalents.

### Phase 2: Remove Notion Imports and Calls

**Objective:** Remove all imports and function calls to Notion modules from remaining files.

**Key Import Patterns to Remove:**
- `from amprenta_rag.clients.notion_client import notion_headers, get_page_text`
- `from amprenta_rag.ingestion.notion_pages import ensure_literature_page, create_rag_chunk_page, extract_page_content`
- `from amprenta_rag.ingestion.dataset_notion_utils import fetch_dataset_page, update_dataset_embedding_metadata, update_dataset_scientific_metadata`
- `from amprenta_rag.analysis.pathway_notion_integration import create_or_update_pathway_page`
- `from amprenta_rag.chemistry.notion_integration import create_compound_feature_page`
- `from amprenta_rag.migration.dual_write import DualWriteManager`

**Files Requiring Import Removal:**
1. `amprenta_rag/query/rag/query.py` - Remove `fetch_dataset_page`, `extract_page_content`
2. `amprenta_rag/ingestion/dataset_ingestion.py` - Remove all Notion utils imports
3. All omics embedding modules (`metabolomics/embedding.py`, `proteomics/embedding.py`, etc.) - Remove `update_dataset_embedding_metadata` calls
4. All omics ingestion modules - Remove Notion page creation calls
5. `amprenta_rag/query/cross_omics/helpers.py` - Remove `fetch_notion_page`, `notion_headers`
6. `amprenta_rag/query/cross_omics/program_summary.py` - Remove Notion fetching
7. `amprenta_rag/query/cross_omics/signature_summary.py` - Remove Notion fetching
8. `amprenta_rag/query/cross_omics/dataset_summary.py` - Remove Notion fetching
9. `amprenta_rag/rag/postgres_resolver.py` - Remove Notion ID lookups
10. `amprenta_rag/analysis/enrichment.py` - Remove pathway Notion integration

**Action Required:**
- Remove import statements
- Remove function calls to Notion APIs
- Replace with Postgres queries where data is needed
- For metadata updates (embedding IDs, etc.), simply remove - data is in Pinecone/Postgres

### Phase 3: Refactor Cross-Omics Functions to Use Postgres

**Objective:** Update cross-omics summary functions to accept Postgres UUIDs instead of Notion page IDs.

**Files to Refactor:**

1. **`amprenta_rag/query/cross_omics/program_summary.py`**
   - Change function signature: `cross_omics_program_summary(program_id: UUID)` instead of `program_page_id: str`
   - Fetch program from Postgres: `db.query(ProgramModel).filter(ProgramModel.id == program_id).first()`
   - Extract related experiments/datasets from Postgres relationships: `program.experiments`, `program.datasets`
   - Build context from Postgres model fields instead of Notion properties
   - Remove all `fetch_notion_page()` calls

2. **`amprenta_rag/query/cross_omics/signature_summary.py`**
   - Change function signature: `cross_omics_signature_summary(signature_id: UUID)`
   - Fetch signature from Postgres
   - Extract components from `signature.components` relationship
   - Extract related datasets from `signature.datasets` relationship

3. **`amprenta_rag/query/cross_omics/dataset_summary.py`**
   - Change function signature: `cross_omics_dataset_summary(dataset_id: UUID)`
   - Fetch dataset from Postgres
   - Extract features from `dataset.features` relationship
   - Extract signatures from `dataset.signatures` relationship
   - Extract experiments from `dataset.experiments` relationship

4. **`amprenta_rag/query/cross_omics/helpers.py`**
   - Remove `fetch_notion_page()` function
   - Add `fetch_postgres_entity()` helper that takes UUID and entity type
   - Update `extract_relation_ids()` to work with Postgres relationship objects
   - Remove all Notion client imports

**Database Session:** These functions will need a database session parameter. Add `db: Session` parameter to all cross-omics summary functions.

### Phase 4: Update RAG Query Functions

**Objective:** Remove Notion dependencies from RAG query pipeline.

**Files to Update:**

1. **`amprenta_rag/query/rag/query.py`**
   - Remove `fetch_dataset_page()` import and calls
   - Remove `extract_page_content()` import and calls
   - For `signature_similarity_query()`: Change parameter from `dataset_page_id: str` to `dataset_id: UUID`
   - Fetch dataset from Postgres instead of Notion
   - Extract features from Postgres `dataset.features` relationship
   - Remove Notion page content extraction

2. **`amprenta_rag/query/rag/chunk_collection.py`**
   - Remove `get_page_text()` calls (from notion_client)
   - Use Pinecone metadata or Postgres data for chunk text
   - Verify `collect_chunks()` works without Notion

3. **`amprenta_rag/rag/hybrid_chunk_collection.py`**
   - Remove hybrid Postgres + Notion logic
   - Use Postgres-only chunk collection
   - Remove `get_page_text()` calls

### Phase 5: Update Ingestion Modules

**Objective:** Remove Notion sync from all ingestion pipelines.

**Files to Update:**

1. **Omics Ingestion Modules:**
   - `amprenta_rag/ingestion/lipidomics/ingestion.py`
   - `amprenta_rag/ingestion/metabolomics/ingestion.py`
   - `amprenta_rag/ingestion/proteomics/ingestion.py`
   - `amprenta_rag/ingestion/transcriptomics/ingestion.py`
   
   **Changes:**
   - Remove `create_*_dataset_page()` functions (Notion page creation)
   - Remove `notion_page_id` parameters from ingestion functions
   - Remove Notion sync code paths
   - Ensure all ingestion uses Postgres only

2. **`amprenta_rag/ingestion/dataset_ingestion.py`**
   - Remove `notion_page_id` parameter requirement
   - Change to accept `dataset_id: UUID` instead
   - Fetch dataset from Postgres
   - Remove `fetch_dataset_page()` calls
   - Remove `update_dataset_embedding_metadata()` calls
   - Remove `update_dataset_scientific_metadata()` calls

3. **`amprenta_rag/ingestion/experiments_ingestion.py`**
   - Remove Notion dependencies
   - Use Postgres experiment model

4. **`amprenta_rag/ingestion/email_ingestion.py`**
   - Remove Notion email DB operations
   - Use Postgres for email storage (if needed)

5. **`amprenta_rag/ingestion/zotero/ingestion.py`**
   - Remove Notion literature page creation
   - Use Postgres for literature storage

6. **Omics Embedding Modules:**
   - `amprenta_rag/ingestion/metabolomics/embedding.py`
   - `amprenta_rag/ingestion/proteomics/embedding.py`
   - `amprenta_rag/ingestion/lipidomics/embedding.py`
   - `amprenta_rag/ingestion/transcriptomics/embedding.py`
   
   **Changes:**
   - Remove `update_dataset_embedding_metadata()` import and calls
   - Embedding metadata is stored in Pinecone, no need to update Notion

### Phase 6: Update Database Models

**Objective:** Remove `notion_page_id` columns from all models.

**File:** `amprenta_rag/database/models.py`

**Changes:**
- Remove `notion_page_id` column from:
  - `Program` model (line ~103)
  - `Experiment` model (line ~139)
  - `Feature` model (line ~162)
  - `Dataset` model (line ~203)
  - `Signature` model (line ~257)

**Migration Required:**
- Create Alembic migration file: `alembic/versions/XXXX_remove_notion_page_id_columns.py`
- Migration should:
  - Drop `notion_page_id` columns from all tables
  - Handle any foreign key constraints if they exist
  - Test migration on development database first

### Phase 7: Update Configuration

**Objective:** Remove all Notion configuration.

**File:** `amprenta_rag/config.py`

**Changes:**
- Remove all Notion database IDs:
  - `NOTION_EMAIL_DB_ID`
  - `NOTION_RAG_DB_ID`
  - `NOTION_EXP_DATA_DB_ID`
  - `NOTION_METABOLITE_FEATURES_DB_ID`
  - `NOTION_PROTEIN_FEATURES_DB_ID`
  - `NOTION_GENE_FEATURES_DB_ID`
  - `NOTION_SIGNATURE_DB_ID`
  - `NOTION_PROGRAMS_DB_ID`
  - `NOTION_EXPERIMENTS_DB_ID`
  - All other Notion DB IDs

- Make Notion API key completely optional:
  - Change `NOTION_API_KEY` to optional (no error if missing)
  - Remove `ENABLE_NOTION_SYNC` check (no longer needed)

- Remove Notion config class if it exists:
  - Remove `NotionConfig` dataclass
  - Remove Notion base URL configuration

### Phase 8: Update API Schemas

**Objective:** Remove `notion_page_id` fields from API response schemas.

**File:** `amprenta_rag/api/schemas.py`

**Changes:**
- Remove `notion_page_id` field from:
  - `Program` response schema
  - `Experiment` response schema
  - `Dataset` response schema
  - `Feature` response schema
  - `Signature` response schema

### Phase 9: Update Scripts

**Objective:** Remove Notion options from all scripts.

**Files to Update:**

1. **`scripts/harvest_repository_study.py`**
   - Remove `--create-notion` option
   - Remove `find_existing_dataset_page()` function (Notion search)
   - Use Postgres dataset lookup instead: `db.query(DatasetModel).filter(DatasetModel.external_ids.contains({repository: study_id})).first()`

2. **`scripts/harvest_mw_studies.py`**
   - Remove `create_notion` parameter
   - Remove Notion page creation

3. **`scripts/list_cross_omics_test_ids.py`**
   - Update to list Postgres UUIDs instead of Notion page IDs
   - Query Postgres: `db.query(ProgramModel).all()` etc.

4. **All ingestion scripts** (`scripts/ingest_*.py`)
   - Remove `--create-notion` options
   - Remove `notion_page_id` parameters
   - Use Postgres UUIDs only

### Phase 10: Update Reporting and Analysis

**Files to Update:**

1. **`amprenta_rag/reporting/evidence_report.py`**
   - Remove `write_evidence_report_to_notion()` function
   - Keep Postgres-based report generation only

2. **`amprenta_rag/analysis/dataset_comparison.py`**
   - Remove Notion dataset fetching
   - Use Postgres dataset queries

3. **`amprenta_rag/analysis/program_signature_maps.py`**
   - Remove Notion program/signature fetching
   - Use Postgres queries

4. **`amprenta_rag/analysis/enrichment.py`**
   - Remove pathway Notion integration imports
   - Remove calls to `create_or_update_pathway_page()`

### Phase 11: Update Maintenance Scripts

**Files to Update:**

1. **`amprenta_rag/maintenance/verify.py`**
   - Remove Notion verification checks
   - Add Postgres verification checks

2. **`amprenta_rag/maintenance/zotero_universe.py`**
   - Remove Notion literature DB operations
   - Use Postgres for literature tracking

### Phase 12: Update Tests

**Files to Update:**

1. **`amprenta_rag/tests/ingestion/test_cross_omics_helpers.py`**
   - Update to test Postgres-based helpers
   - Remove Notion mocking
   - Add Postgres test fixtures

2. **All test files with Notion mocks**
   - Remove Notion API mocks
   - Add Postgres test fixtures
   - Update test data to use Postgres UUIDs

### Phase 13: Clean Up Imports Globally

**Objective:** Remove all remaining Notion imports.

**Action:**
- Search entire codebase for: `notion_client`, `notion_pages`, `dataset_notion_utils`, `pathway_notion_integration`, `chemistry.notion_integration`, `dual_write`
- Remove all import statements
- Remove all function calls to these modules

## Implementation Order

1. **Week 1:** Phases 1-2 (Delete modules, remove imports)
2. **Week 2:** Phases 3-4 (Refactor cross-omics, update RAG)
3. **Week 3:** Phases 5-6 (Update ingestion, database models)
4. **Week 4:** Phases 7-9 (Config, API, scripts)
5. **Week 5:** Phases 10-13 (Reporting, tests, cleanup)

## Verification Checklist

After each phase, verify:
- [ ] No import errors
- [ ] No undefined function calls
- [ ] Tests pass (update tests as needed)
- [ ] No Notion API calls in logs
- [ ] System works with Postgres only

## Key Principles

1. **Postgres is Source of Truth:** All data operations use Postgres
2. **Pinecone for Embeddings:** Embedding metadata is in Pinecone, not Notion
3. **Remove, Don't Replace:** For optional features (like Notion metadata updates), simply remove them
4. **Test Thoroughly:** After each phase, test affected functionality
5. **Update Documentation:** Update any docs that reference Notion

## Questions for Architect

If you encounter issues or need clarification:
1. What should replace Notion page ID lookups? → Use Postgres UUIDs
2. How to handle chunk text retrieval? → Use Pinecone metadata or Postgres data
3. What about backward compatibility? → Not needed - migration is complete
4. Should we keep Notion code commented? → No, delete it completely

## Success Criteria

- Zero Notion API calls in codebase
- All functions work with Postgres UUIDs only
- No Notion dependencies in requirements
- All tests pass without Notion
- Documentation reflects Postgres-only architecture
- System operates completely independently of Notion

