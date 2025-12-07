# Complete Notion Removal Plan - Agent Delegation

## Overview

Remove all Notion dependencies from the codebase. Postgres is now the Source of Truth (100% migration complete).

**Scope:**
- Files with Notion references: 114 files found
- Total Notion imports/references: 367 matches across 78 files
- Estimated effort: 4-5 weeks

## Agent Roles

- **Implementor**: Makes code changes, removes imports, refactors functions
- **Tester**: Designs and writes tests, validates functionality  
- **Reviewer**: Reviews code changes for correctness, safety, style
- **Automator**: Handles file deletion, cleanup, script generation

---

## Phase 1: Delete Notion-Specific Modules

### Task 1.1: Verify Dependencies (Automator)
**Agent:** Automator  
**Action:** Search codebase for imports of Notion modules to verify they're safe to delete  
**Files:** All Python files  
**Output:** List of files that import each Notion module  
**Acceptance:** Confirmed list of files that will need import removal

### Task 1.2: Delete Notion Modules (Automator)
**Agent:** Automator  
**Action:** Delete the following files:
1. `amprenta_rag/clients/notion_client.py`
2. `amprenta_rag/ingestion/notion_pages.py`
3. `amprenta_rag/ingestion/dataset_notion_utils.py`
4. `amprenta_rag/ingestion/signature_notion_crud.py`
5. `amprenta_rag/analysis/pathway_notion_integration.py`
6. `amprenta_rag/chemistry/notion_integration.py`
7. `amprenta_rag/migration/dual_write.py`
**Acceptance:** Files deleted, no broken imports remain

### Task 1.3: Review Deletion (Reviewer)
**Agent:** Reviewer  
**Action:** Review that deleted modules had no critical functionality  
**Acceptance:** Confirmed safe to delete

---

## Phase 2: Remove Notion Imports and Calls

### Task 2.1: Remove Imports from RAG Query Module (Implementor)
**Agent:** Implementor  
**Files:** `amprenta_rag/query/rag/query.py`  
**Action:** Remove `fetch_dataset_page`, `extract_page_content` imports and calls  
**Review:** Reviewer  
**Test:** Tester - Verify RAG queries still work  
**Acceptance:** No Notion imports, code compiles

### Task 2.2: Remove Imports from Dataset Ingestion (Implementor)
**Agent:** Implementor  
**Files:** `amprenta_rag/ingestion/dataset_ingestion.py`  
**Action:** Remove all Notion utils imports (`fetch_dataset_page`, `update_dataset_embedding_metadata`, `update_dataset_scientific_metadata`)  
**Review:** Reviewer  
**Test:** Tester - Verify dataset ingestion still works  
**Acceptance:** No Notion imports, code compiles

### Task 2.3: Remove Embedding Metadata Updates from Metabolomics (Implementor)
**Agent:** Implementor  
**Files:** `amprenta_rag/ingestion/metabolomics/embedding.py`  
**Action:** Remove `update_dataset_embedding_metadata` import and call  
**Review:** Reviewer  
**Test:** Tester - Verify embedding still works  
**Acceptance:** No Notion calls, embedding succeeds

### Task 2.4: Remove Embedding Metadata Updates from Proteomics (Implementor)
**Agent:** Implementor  
**Files:** `amprenta_rag/ingestion/proteomics/embedding.py`  
**Action:** Remove `update_dataset_embedding_metadata` import and call  
**Review:** Reviewer  
**Test:** Tester - Verify embedding still works  
**Acceptance:** No Notion calls, embedding succeeds

### Task 2.5: Remove Embedding Metadata Updates from Lipidomics (Implementor)
**Agent:** Implementor  
**Files:** `amprenta_rag/ingestion/lipidomics/embedding.py`  
**Action:** Remove `update_dataset_embedding_metadata` import and call  
**Review:** Reviewer  
**Test:** Tester - Verify embedding still works  
**Acceptance:** No Notion calls, embedding succeeds

### Task 2.6: Remove Embedding Metadata Updates from Transcriptomics (Implementor)
**Agent:** Implementor  
**Files:** `amprenta_rag/ingestion/transcriptomics/embedding.py`  
**Action:** Remove `update_dataset_embedding_metadata` import and call  
**Review:** Reviewer  
**Test:** Tester - Verify embedding still works  
**Acceptance:** No Notion calls, embedding succeeds

### Task 2.7: Remove Notion Calls from Cross-Omics Helpers (Implementor)
**Agent:** Implementor  
**Files:** `amprenta_rag/query/cross_omics/helpers.py`  
**Action:** Remove `fetch_notion_page`, `notion_headers` imports and functions  
**Review:** Reviewer  
**Test:** Tester - Verify cross-omics helpers work  
**Acceptance:** No Notion imports, code compiles

### Task 2.8: Remove Notion Calls from Program Summary (Implementor)
**Agent:** Implementor  
**Files:** `amprenta_rag/query/cross_omics/program_summary.py`  
**Action:** Remove Notion fetching calls (will be replaced in Phase 3)  
**Review:** Reviewer  
**Acceptance:** No Notion calls remain

### Task 2.9: Remove Notion Calls from Signature Summary (Implementor)
**Agent:** Implementor  
**Files:** `amprenta_rag/query/cross_omics/signature_summary.py`  
**Action:** Remove Notion fetching calls (will be replaced in Phase 3)  
**Review:** Reviewer  
**Acceptance:** No Notion calls remain

### Task 2.10: Remove Notion Calls from Dataset Summary (Implementor)
**Agent:** Implementor  
**Files:** `amprenta_rag/query/cross_omics/dataset_summary.py`  
**Action:** Remove Notion fetching calls (will be replaced in Phase 3)  
**Review:** Reviewer  
**Acceptance:** No Notion calls remain

### Task 2.11: Remove Notion Calls from Postgres Resolver (Implementor)
**Agent:** Implementor  
**Files:** `amprenta_rag/rag/postgres_resolver.py`  
**Action:** Remove Notion ID lookup functions  
**Review:** Reviewer  
**Test:** Tester - Verify Postgres resolution still works  
**Acceptance:** No Notion calls remain

### Task 2.12: Remove Notion Calls from Enrichment Analysis (Implementor)
**Agent:** Implementor  
**Files:** `amprenta_rag/analysis/enrichment.py`  
**Action:** Remove pathway Notion integration imports and calls  
**Review:** Reviewer  
**Test:** Tester - Verify enrichment analysis still works  
**Acceptance:** No Notion calls remain

---

## Phase 3: Refactor Cross-Omics Functions to Use Postgres

### Task 3.1: Refactor Cross-Omics Helpers (Implementor)
**Agent:** Implementor  
**Files:** `amprenta_rag/query/cross_omics/helpers.py`  
**Action:** 
- Remove `fetch_notion_page()` function
- Add `fetch_postgres_entity(uuid: UUID, entity_type: str, db: Session)` helper
- Update `extract_relation_ids()` to work with Postgres relationship objects
**Review:** Reviewer  
**Test:** Tester - Write tests for new Postgres helpers  
**Acceptance:** Helpers work with Postgres models

### Task 3.2: Refactor Program Summary (Implementor)
**Agent:** Implementor  
**Files:** `amprenta_rag/query/cross_omics/program_summary.py`  
**Action:** 
- Change signature: `cross_omics_program_summary(program_id: UUID, db: Session)`
- Fetch program from Postgres: `db.query(ProgramModel).filter(ProgramModel.id == program_id).first()`
- Extract related experiments/datasets from Postgres relationships
- Build context from Postgres model fields
**Review:** Reviewer  
**Test:** Tester - Test program summary with Postgres UUIDs  
**Acceptance:** Function works with Postgres UUIDs

### Task 3.3: Refactor Signature Summary (Implementor)
**Agent:** Implementor  
**Files:** `amprenta_rag/query/cross_omics/signature_summary.py`  
**Action:** 
- Change signature: `cross_omics_signature_summary(signature_id: UUID, db: Session)`
- Fetch signature from Postgres
- Extract components from `signature.components` relationship
- Extract related datasets from `signature.datasets` relationship
**Review:** Reviewer  
**Test:** Tester - Test signature summary with Postgres UUIDs  
**Acceptance:** Function works with Postgres UUIDs

### Task 3.4: Refactor Dataset Summary (Implementor)
**Agent:** Implementor  
**Files:** `amprenta_rag/query/cross_omics/dataset_summary.py`  
**Action:** 
- Change signature: `cross_omics_dataset_summary(dataset_id: UUID, db: Session)`
- Fetch dataset from Postgres
- Extract features from `dataset.features` relationship
- Extract signatures from `dataset.signatures` relationship
- Extract experiments from `dataset.experiments` relationship
**Review:** Reviewer  
**Test:** Tester - Test dataset summary with Postgres UUIDs  
**Acceptance:** Function works with Postgres UUIDs

---

## Phase 4: Update RAG Query Functions

### Task 4.1: Update RAG Query Function (Implementor)
**Agent:** Implementor  
**Files:** `amprenta_rag/query/rag/query.py`  
**Action:** 
- Remove `fetch_dataset_page()` import and calls
- Remove `extract_page_content()` import and calls
- For `signature_similarity_query()`: Change parameter from `dataset_page_id: str` to `dataset_id: UUID`
- Fetch dataset from Postgres instead of Notion
- Extract features from Postgres `dataset.features` relationship
**Review:** Reviewer  
**Test:** Tester - Test signature similarity query with Postgres UUIDs  
**Acceptance:** Function works with Postgres UUIDs

### Task 4.2: Update Chunk Collection (Implementor)
**Agent:** Implementor  
**Files:** `amprenta_rag/query/rag/chunk_collection.py`  
**Action:** 
- Remove `get_page_text()` calls (from notion_client)
- Use Pinecone metadata or Postgres data for chunk text
**Review:** Reviewer  
**Test:** Tester - Verify `collect_chunks()` works without Notion  
**Acceptance:** Chunk collection works without Notion

### Task 4.3: Update Hybrid Chunk Collection (Implementor)
**Agent:** Implementor  
**Files:** `amprenta_rag/rag/hybrid_chunk_collection.py`  
**Action:** 
- Remove hybrid Postgres + Notion logic
- Use Postgres-only chunk collection
- Remove `get_page_text()` calls
**Review:** Reviewer  
**Test:** Tester - Verify hybrid collection works Postgres-only  
**Acceptance:** Hybrid collection works Postgres-only

---

## Phase 5: Update Ingestion Modules

### Task 5.1: Remove Notion from Lipidomics Ingestion (Implementor)
**Agent:** Implementor  
**Files:** `amprenta_rag/ingestion/lipidomics/ingestion.py`  
**Action:** 
- Remove `create_lipidomics_dataset_page()` function
- Remove `notion_page_id` parameters
- Remove Notion sync code paths
**Review:** Reviewer  
**Test:** Tester - Verify lipidomics ingestion works Postgres-only  
**Acceptance:** Ingestion works Postgres-only

### Task 5.2: Remove Notion from Metabolomics Ingestion (Implementor)
**Agent:** Implementor  
**Files:** `amprenta_rag/ingestion/metabolomics/ingestion.py`  
**Action:** 
- Remove `create_metabolomics_dataset_page()` function
- Remove `notion_page_id` parameters
- Remove Notion sync code paths
**Review:** Reviewer  
**Test:** Tester - Verify metabolomics ingestion works Postgres-only  
**Acceptance:** Ingestion works Postgres-only

### Task 5.3: Remove Notion from Proteomics Ingestion (Implementor)
**Agent:** Implementor  
**Files:** `amprenta_rag/ingestion/proteomics/ingestion.py`  
**Action:** 
- Remove `create_proteomics_dataset_page()` function
- Remove `notion_page_id` parameters
- Remove Notion sync code paths
**Review:** Reviewer  
**Test:** Tester - Verify proteomics ingestion works Postgres-only  
**Acceptance:** Ingestion works Postgres-only

### Task 5.4: Remove Notion from Transcriptomics Ingestion (Implementor)
**Agent:** Implementor  
**Files:** `amprenta_rag/ingestion/transcriptomics/ingestion.py`  
**Action:** 
- Remove `create_transcriptomics_dataset_page()` function
- Remove `notion_page_id` parameters
- Remove Notion sync code paths
**Review:** Reviewer  
**Test:** Tester - Verify transcriptomics ingestion works Postgres-only  
**Acceptance:** Ingestion works Postgres-only

### Task 5.5: Refactor Dataset Ingestion (Implementor)
**Agent:** Implementor  
**Files:** `amprenta_rag/ingestion/dataset_ingestion.py`  
**Action:** 
- Remove `notion_page_id` parameter requirement
- Change to accept `dataset_id: UUID` instead
- Fetch dataset from Postgres
- Remove all Notion function calls
**Review:** Reviewer  
**Test:** Tester - Verify dataset ingestion works with Postgres UUIDs  
**Acceptance:** Ingestion works with Postgres UUIDs

### Task 5.6: Remove Notion from Other Ingestion Modules (Implementor)
**Agent:** Implementor  
**Files:** 
- `amprenta_rag/ingestion/experiments_ingestion.py`
- `amprenta_rag/ingestion/email_ingestion.py`
- `amprenta_rag/ingestion/zotero/ingestion.py`
- `amprenta_rag/ingestion/signature_ingestion.py`
- `amprenta_rag/ingestion/signature_matching/signature_loader.py`
**Action:** Remove Notion dependencies from each module  
**Review:** Reviewer  
**Test:** Tester - Verify each module works Postgres-only  
**Acceptance:** All modules work Postgres-only

---

## Phase 6: Update Database Models

### Task 6.1: Remove notion_page_id from Models (Implementor)
**Agent:** Implementor  
**Files:** `amprenta_rag/database/models.py`  
**Action:** Remove `notion_page_id` column from:
- `Program` model
- `Experiment` model
- `Feature` model
- `Dataset` model
- `Signature` model
**Review:** Reviewer - Verify no breaking changes  
**Acceptance:** Models updated, no notion_page_id columns

### Task 6.2: Create Alembic Migration (Implementor)
**Agent:** Implementor  
**Files:** `alembic/versions/XXXX_remove_notion_page_id_columns.py`  
**Action:** Create migration to drop `notion_page_id` columns from all tables  
**Review:** Reviewer - Verify migration safety  
**Test:** Tester - Test migration on development database  
**Acceptance:** Migration works correctly

---

## Phase 7: Update Configuration

### Task 7.1: Remove Notion Config (Implementor)
**Agent:** Implementor  
**Files:** `amprenta_rag/config.py`  
**Action:** 
- Remove all Notion database IDs
- Make `NOTION_API_KEY` optional (no error if missing)
- Remove `ENABLE_NOTION_SYNC` check
- Remove `NotionConfig` dataclass if exists
**Review:** Reviewer  
**Test:** Tester - Verify config loads without Notion settings  
**Acceptance:** Config works without Notion

---

## Phase 8: Update API Schemas

### Task 8.1: Remove notion_page_id from Schemas (Implementor)
**Agent:** Implementor  
**Files:** `amprenta_rag/api/schemas.py`  
**Action:** Remove `notion_page_id` field from:
- `Program` response schema
- `Experiment` response schema
- `Dataset` response schema
- `Feature` response schema
- `Signature` response schema
**Review:** Reviewer  
**Test:** Tester - Verify API responses don't include notion_page_id  
**Acceptance:** Schemas updated, API works

---

## Phase 9: Update Scripts

### Task 9.1: Update Harvest Repository Script (Implementor)
**Agent:** Implementor  
**Files:** `scripts/harvest_repository_study.py`  
**Action:** 
- Remove `--create-notion` option
- Remove `find_existing_dataset_page()` function
- Use Postgres dataset lookup instead
**Review:** Reviewer  
**Test:** Tester - Verify harvest script works Postgres-only  
**Acceptance:** Script works Postgres-only

### Task 9.2: Update Other Harvest Scripts (Implementor)
**Agent:** Implementor  
**Files:** `scripts/harvest_mw_studies.py`, `scripts/list_cross_omics_test_ids.py`  
**Action:** Remove Notion options, use Postgres queries  
**Review:** Reviewer  
**Test:** Tester - Verify scripts work Postgres-only  
**Acceptance:** Scripts work Postgres-only

### Task 9.3: Update Ingestion Scripts (Automator)
**Agent:** Automator  
**Files:** All `scripts/ingest_*.py` files  
**Action:** Remove `--create-notion` options and `notion_page_id` parameters  
**Review:** Reviewer  
**Acceptance:** Scripts updated

---

## Phase 10: Update Reporting and Analysis

### Task 10.1: Update Evidence Report (Implementor)
**Agent:** Implementor  
**Files:** `amprenta_rag/reporting/evidence_report.py`  
**Action:** Remove `write_evidence_report_to_notion()` function  
**Review:** Reviewer  
**Test:** Tester - Verify report generation still works  
**Acceptance:** Reports generate Postgres-only

### Task 10.2: Update Analysis Modules (Implementor)
**Agent:** Implementor  
**Files:** 
- `amprenta_rag/analysis/dataset_comparison.py`
- `amprenta_rag/analysis/program_signature_maps.py`
**Action:** Remove Notion fetching, use Postgres queries  
**Review:** Reviewer  
**Test:** Tester - Verify analysis modules work Postgres-only  
**Acceptance:** Analysis works Postgres-only

---

## Phase 11: Update Maintenance Scripts

### Task 11.1: Update Maintenance Scripts (Implementor)
**Agent:** Implementor  
**Files:** 
- `amprenta_rag/maintenance/verify.py`
- `amprenta_rag/maintenance/zotero_universe.py`
**Action:** Remove Notion operations, add Postgres checks  
**Review:** Reviewer  
**Test:** Tester - Verify maintenance scripts work  
**Acceptance:** Maintenance scripts work Postgres-only

---

## Phase 12: Update Tests

### Task 12.1: Design Test Strategy (Tester)
**Agent:** Tester  
**Action:** Design test strategy for Postgres-based functionality  
**Output:** Test plan document  
**Acceptance:** Comprehensive test plan created

### Task 12.2: Update Cross-Omics Tests (Tester)
**Agent:** Tester  
**Files:** `amprenta_rag/tests/ingestion/test_cross_omics_helpers.py`  
**Action:** Update tests to use Postgres fixtures, remove Notion mocks  
**Review:** Reviewer  
**Acceptance:** Tests pass with Postgres

### Task 12.3: Update All Tests with Notion Mocks (Tester)
**Agent:** Tester  
**Files:** All test files  
**Action:** Remove Notion API mocks, add Postgres test fixtures  
**Review:** Reviewer  
**Acceptance:** All tests pass without Notion

---

## Phase 13: Clean Up Imports Globally

### Task 13.1: Global Import Cleanup (Automator)
**Agent:** Automator  
**Action:** Search entire codebase for remaining Notion imports and remove them  
**Files:** All Python files  
**Review:** Reviewer  
**Acceptance:** Zero Notion imports remain

### Task 13.2: Final Verification (Automator)
**Agent:** Automator  
**Action:** Generate report of all Notion references remaining (should be zero)  
**Output:** Verification report  
**Acceptance:** Zero Notion references confirmed

---

## Implementation Order

1. **Week 1:** Phases 1-2 (Delete modules, remove imports)
2. **Week 2:** Phases 3-4 (Refactor cross-omics, update RAG)
3. **Week 3:** Phases 5-6 (Update ingestion, database models)
4. **Week 4:** Phases 7-9 (Config, API, scripts)
5. **Week 5:** Phases 10-13 (Reporting, tests, cleanup)

## Verification Checklist

After each task, verify:
- [ ] No import errors
- [ ] No undefined function calls
- [ ] Tests pass (update tests as needed)
- [ ] No Notion API calls in logs
- [ ] System works with Postgres only

## Success Criteria

- Zero Notion API calls in codebase
- All functions work with Postgres UUIDs only
- No Notion dependencies in requirements
- All tests pass without Notion
- Documentation reflects Postgres-only architecture
- System operates completely independently of Notion

