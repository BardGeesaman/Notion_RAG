# Amprenta RAG System – Engineering Roadmap

This file tracks the **current state** and **next steps** for the Amprenta RAG system:
- Notion + Zotero + Pinecone + OpenAI
- Ingestion, metadata, query, classification
- Refactors, tests, logging, and features

Whenever we update the code and create a new zip, this file should be kept roughly in sync so both humans and AI can see the state at a glance.

---

## 1. Goals

- Support Amprenta’s ceramide/sphingolipid neurodegeneration program with a **high-quality, explainable RAG system**.
- Keep the codebase **modular, testable, and comprehensible** to non-software scientists.
- Make it easy to add:
  - New data sources (Experiments, Clinical, etc.)
  - New metadata fields (lipid signatures, phenotypes, model systems)
  - New downstream tools (dashboards, analyses)

---

## 2. Current State (High-Level)

### ✅ Ingestion

- **Zotero → Notion → Pinecone → RAG**:
  - Full PDFs + Zotero notes ingested.
  - Attachments & notes chunked and embedded.
  - Idempotent with md5 (attachments) and hash (notes).
  - Snapshot HTML attachments correctly skipped.

- **Email → Notion → Pinecone**:
  - Emails/notes can be ingested and chunked.
  - Used mainly for free-form insights, meeting notes, etc.

### ✅ Metadata

- For literature chunks, metadata now includes:
  - `doc_id`, `doc_source`, `doc_type`
  - `year`, `journal`, `importance`
  - `diseases`, `targets` (molecular targets), `modality`, `stage`
  - `model_systems`, `biomarker_role`
  - `lipids_raw` (lipid species mentions)
  - `lipid_signatures`, `phenotype_axes`, `matrix`, `treatment_arms`
  - `zotero_item_key`, `notion_chunk_page_id`, legacy tags, etc.

### ✅ Query / RAG

- `rag_query.py`:
  - Semantic search over Pinecone.
  - Filters: `--source-type`, `--disease`, `--target`, `--lipid`, `--signature`, `--tag`.
  - Uses the correct Pinecone namespace.
  - Synthesizes answers from top-ranked chunks with context display.

### ✅ Classifier

- `classify_literature_metadata.py`:
  - Uses OpenAI to classify each paper into:
    - Disease, Targets, Modality, Stage, Model Systems, Biomarker Role.
    - Lipid species (raw), lipid signatures, phenotype axes, matrix, treatment arms.
  - **Targets** now restricted to **molecular targets only** (e.g., SPTLC1, DEGS1, CERS2, ORMDL3).
  - Lipid species (e.g., ceramide, sphingomyelin) live in `Lipid Species (raw)` instead of Targets.
  - Has a `--force` mode to reclassify even if metadata already exists.

---

## 3. Refactor & Cleanup Status

### ✅ Completed / In-Progress

- [x] Introduced `amprenta_rag` as a package.
- [x] Split Zotero API into `ingestion/zotero_api.py`.
- [x] Split text extraction + boilerplate logic into `ingestion/text_extraction.py`.
- [x] Introduced `ingestion/zotero_ingest.py` as the orchestrator.
- [x] Added `sanitize_metadata` and centralized Pinecone metadata cleanup.
- [x] Fixed namespace issues in query code.
- [x] Fixed classifier outputs and Targets vs Lipids separation.
- [x] Confirmed ingestion + query still work after initial refactors.

### ✅ Completed Refactoring Phases 1-5 (December 2025)

**Phase 1-2: Ingestion Pipeline Cleanup**
- [x] Moved Notion-specific helpers to `ingestion/notion_pages.py`:
  - `ensure_literature_page`, `update_literature_page`, `create_rag_chunk_page`
  - `fetch_not_embedded_emails`, `extract_page_content`, `update_email_page`
- [x] Moved semantic metadata helpers to `ingestion/metadata_semantic.py`:
  - `get_literature_semantic_metadata`, `get_email_semantic_metadata`
- [x] Consolidated Pinecone idempotency in `ingestion/pinecone_utils.py`:
  - `attachment_already_ingested`, `note_already_ingested` (already there)
- [x] Modularized `email_ingestion.py` using shared helpers
- [x] Removed ~200 lines of duplicate code

**Phase 3: Query Engine Refactoring**
- [x] Split query engine into:
  - `query/pinecone_query.py` (embed + raw Pinecone search)
  - `query/rag_engine.py` (filtering + answer synthesis)
  - `query/rag_query_engine.py` (compatibility wrapper)
- [x] Maintained 100% backward compatibility

**Phase 4: Logging & Error Handling**
- [x] Normalized logging with consistent prefixes ([INGEST][ZOTERO], [NOTION], [PINECONE], [RAG])
- [x] Replaced all `print()` statements with structured logger calls
- [x] Added error logging for all external API calls with proper exception handling

**Phase 5: Code Hygiene**
- [x] Removed unused imports
- [x] Verified consistent formatting
- [x] Config now uses environment variables (secure)

### Optional Future Enhancements

- [ ] Remove legacy `zotero_item.py` file (marked as LEGACY, functionality migrated)
- [ ] Consider splitting classifier module further:
  - `metadata/notion_utils.py`
  - `metadata/classifier_base.py`
  - `metadata/classify_literature.py`

---

## 4. Testing Plan (Lightweight but Useful)

### ✅ Already sketched / partially implemented

- `tests/test_sanitize_metadata.py`:
  - Ensures no `None` or unsupported types are sent to Pinecone.

- `tests/test_metadata_classifier.py`:
  - Ensures Targets ≠ Lipids:
    - Targets do **not** include "ceramide", "sphingomyelin".
    - Lipid terms end up in `Lipid Species (raw)`.

- `tests/test_build_meta_filter.py`:
  - Tests `build_meta_filter()` builds correct Pinecone filter structure:
    - `diseases`, `targets`, `lipid_signatures`
    - `$or` across `lipids` and `lipids_raw`.

### ⏭️ Planned tests

- [ ] Add tests for:
  - `_chunk_text` behavior.
  - Idempotency logic (attachment_already_ingested / note_already_ingested) with mocked Pinecone responses.
  - Query filtering by `--source-type` (Email vs Literature) with synthetic metadata.

### Manual smoke tests (keep doing these after major changes)

- [ ] Run:

  ```bash
  python scripts/ingest_collection.py --collection-key 3RGXZTAY --parent-type Literature
  