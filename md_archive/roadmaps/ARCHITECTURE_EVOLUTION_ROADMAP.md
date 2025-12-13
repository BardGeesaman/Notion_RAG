# Architecture Evolution: Postgres + FastAPI + Frontend

**Status**: 📋 Roadmap Item - Architecture Evolution

**Priority**: 🔥🔥 **MEDIUM-HIGH** (Future Evolution)

**Timeline**: After core features are stable (Weeks 17+)

This document defines the evolution path from the current Notion + SQLite + Pinecone architecture to a more robust Postgres + FastAPI + Frontend architecture.

---

## 0. CORE PRINCIPLES

1. **Postgres becomes the system of record** for structured data:
   - Programs, Experiments, Datasets, Features, Signatures, Compounds, Campaigns.

2. **SQLite remains an internal detail or stepping stone**:
   - May continue to be used for rapid chem/screening prototyping or as a cache.
   - But primary schema should eventually be implemented in Postgres.

3. **Notion remains a human-friendly overlay:**
   - ELN-style experiment pages.
   - SOPs, design docs, decisions.
   - High-level dashboards that can be fed by Postgres (via exports or API integration).

4. **Pinecone remains a semantic index:**
   - Embeddings are generated from text derived from Postgres + Notion, not raw DB rows.

5. Keep **ID semantics stable**:
   - Use stable IDs that can be stored in all systems (Postgres, Notion, RAG metadata).

---

## 1. EXTRACT & FORMALIZE THE DOMAIN MODEL

### 1.1 Domain Models (Pydantic/Dataclasses)

Create a set of Python domain models that represent the core entities, independent of storage:

**Location**: `amprenta_rag/models/domain.py`

**Core Entities**:
- `Program` - Research programs
- `Experiment` - Experimental designs
- `Dataset` - Experimental Data Assets (lipidomics, metabolomics, proteomics, transcriptomics)
- `Feature` - Genes, Proteins, Metabolites, Lipids, Compounds
- `Signature` - Multi-omics signatures
- `SignatureComponent` - Individual signature components
- `Compound` - Chemical compounds
- `ScreeningCampaign` - HTS campaigns
- `BiochemicalResult` - Biochemical assay results

**Requirements**:
- Not depend on Notion or SQLAlchemy directly
- Encode relationships via IDs (e.g., `program_id`, `experiment_ids`, `feature_ids`)
- Provide a stable abstraction layer for ingestion modules

**Goal**: All ingestion logic should operate on these models, not directly on Notion/SQLite.

**Estimated Effort**: Medium (3-4 days)

---

## 2. DESIGN THE POSTGRES SCHEMA

### 2.1 Schema Design

Based on the domain models, design a normalized Postgres schema using SQLAlchemy models.

**Location**: `amprenta_api/db/models.py`

**Tables** (high-level):
- `programs` - Research programs
- `experiments` - Experimental designs
- `datasets` - Experimental Data Assets
- `features` - Unified features table with `feature_type` field (gene/protein/metabolite/lipid/compound)
- `signatures` - Multi-omics signatures
- `signature_components` - Signature components
- `compounds` - Chemical compounds (from SQLite migration)
- `screening_campaigns` - HTS campaigns
- `hts_results` - HTS screening results
- `biochemical_results` - Biochemical assay results
- `compound_program` - Junction table (compound ↔ program)
- `dataset_feature` - Junction table (dataset ↔ feature)
- `dataset_signature` - Junction table (dataset ↔ signature, with scores)

**Key Requirements**:
- Use Alembic for migrations
- Maintain **external IDs** (e.g., `notion_page_id`, `sqlite_compound_id`, `inchi_key`) as columns
- Support migration & syncing between systems
- Indexes on foreign keys and frequently queried fields

**Estimated Effort**: High (5-7 days)

---

## 3. DESIGN THE FASTAPI SERVICE LAYER

### 3.1 API Structure

Create a new app: `amprenta_api/`

**Stack**:
- FastAPI + SQLAlchemy + Pydantic models
- Clear API modules organized by domain

### 3.2 Core API Endpoints

#### Programs:
- `GET /programs` - List all programs
- `GET /programs/{program_id}` - Get program details
- `POST /programs` - Create program (internal use)
- `GET /programs/{program_id}/experiments` - Get program experiments
- `GET /programs/{program_id}/datasets` - Get program datasets
- `GET /programs/{program_id}/signatures` - Get program signatures

#### Experiments:
- `GET /experiments` - List experiments
- `GET /experiments/{experiment_id}` - Get experiment details
- `POST /experiments` - Create experiment (from ingestion/ELN sync)
- `GET /experiments/{experiment_id}/datasets` - Get experiment datasets

#### Datasets:
- `GET /datasets` - List datasets (with filtering)
- `GET /datasets/{dataset_id}` - Get dataset details
- `GET /datasets/{dataset_id}/features` - Get dataset features
- `GET /datasets/{dataset_id}/signatures` - Get dataset signatures

#### Features:
- `GET /features?type=gene&symbol=TP53` - Query features
- `GET /features/{feature_id}` - Get feature details

#### Signatures:
- `GET /signatures` - List signatures
- `GET /signatures/{signature_id}` - Get signature details
- `POST /signatures` - Create signature (from signature ingestion)
- `POST /signatures/{signature_id}/score` - Score signature against dataset

#### Compounds:
- `GET /compounds` - List compounds
- `GET /compounds/{compound_id}` - Get compound details
- `GET /compounds/{compound_id}/experiments` - Get compound experiments
- `GET /compounds/{compound_id}/datasets` - Get compound datasets
- `GET /compounds/{compound_id}/signatures` - Get compound signatures

#### Screening:
- `GET /screening/campaigns` - List campaigns
- `GET /screening/campaigns/{campaign_id}` - Get campaign details
- `GET /screening/campaigns/{campaign_id}/hits` - Get campaign hits
- `POST /screening/campaigns/{campaign_id}/ingest` - Ingest campaign data (future)

**API Design Principles**:
- Thin API layer - model domain, enforce consistency
- Call underlying business logic (ingestion, scoring, etc.)
- JSON-based responses
- Pagination for large lists
- Cacheable where appropriate

**Estimated Effort**: High (6-8 days)

---

## 4. MIGRATION PLAN: NOTION + SQLITE → POSTGRES

### 4.1 Stage 0 — Introspection & Export

**Utilities**:
- Export Notion entities (Programs, Experiments, Datasets, Features, Signatures) to JSON/CSV
- Export all SQLite tables (compounds, hts_results, etc.) to CSV/Parquet
- Store exports in `/data/export/{date}/`

**Scripts**:
- `scripts/export_notion_data.py`
- `scripts/export_sqlite_data.py`

**Estimated Effort**: Low (2 days)

---

### 4.2 Stage 1 — Bootstrap Postgres

**Tasks**:
- Apply initial Alembic migration to create schema
- Load Notion + SQLite exports into Postgres
- Map Notion `page_id` → appropriate ID fields in Postgres
- Map SQLite `compound_id` → `compounds` table

**Migration Scripts**:
- `scripts/migrate_notion_to_postgres.py`
- `scripts/migrate_sqlite_to_postgres.py`

**Estimated Effort**: Medium (4-5 days)

---

### 4.3 Stage 2 — Dual Write (Optional Transition Phase)

**Configuration**:
- `WRITE_TO_NOTION = True/False`
- `WRITE_TO_POSTGRES = True/False`

**Behavior**:
- Ingestion scripts continue writing to Notion (for prototype compatibility)
- **AND** start writing to Postgres as the system of record
- Allows flipping behavior gradually

**Estimated Effort**: Medium (3-4 days)

---

### 4.4 Stage 3 — Pivot to Postgres as SoT

**Tasks**:
- Stop writing ingestion results to Notion (or make it optional)
- Keep Notion as:
  - read-only mirror for curated entities
  - documentation layer
  - high-level views
- Continue writing ingestion outputs only to Postgres

**Estimated Effort**: Low (1-2 days)

---

## 5. RAG INTEGRATION WITH POSTGRES

### 5.1 Text Building Strategy

Once Postgres is live, RAG should build text from:
- **Postgres queries**: datasets, experiments, signatures, features, compounds (structured details, counts, metrics)
- **Notion content**: ELN entries, SOPs, decisions (narratives, methods, interpretation)

### 5.2 RAG Builder Updates

**Location**: `amprenta_rag/ingestion/rag_builders.py`

**Changes**:
- Replace current "Notion-only" text extraction with DB + Notion hybrid
- DB for structured details (counts, metrics, features)
- Notion for narratives (methods, decisions, interpretation)

### 5.3 Metadata Requirements

Metadata written to Pinecone must include:
- `dataset_id` (Postgres ID)
- `experiment_id` (Postgres ID)
- `program_id` (Postgres ID)
- `feature_id` / `signature_id` (Postgres ID)
- `source_type` and `omics_type`
- `notion_page_id` (for reference)

**All RAG reasoning should be updated to rely on these stable IDs.**

**Estimated Effort**: Medium (4-5 days)

---

## 6. FRONT-END ARCHITECTURE (FOR FUTURE)

### 6.1 Recommended Stack

- **Frontend**: Next.js + React + TypeScript
- **Backend**: FastAPI service
- **Auth**: Basic token-based (later OIDC with Auth0/Cognito/etc.)

### 6.2 API Design for Frontend

Design API responses to be:
- JSON-based
- Simple and cacheable
- Paginated where needed (e.g., lists of datasets, compounds)
- RESTful conventions

**Note**: Cursor doesn't build the frontend directly but should design APIs with a frontend in mind.

**Estimated Effort**: N/A (Frontend is separate work)

---

## 7. CODING GUIDELINES FOR FUTURE WORK

### 7.1 Architecture Principles

- Keep **domain models** separated from storage specifics
- Use clear interfaces for:
  - `storage_notion`
  - `storage_postgres`
  - `storage_sqlite`
- Maintain **idempotent** ingestion behavior
- Avoid deep Notion coupling in new code:
  - treat Notion as a *client* or a *view*, not as SoT

### 7.2 Testing Requirements

Add tests for:
- Migration scripts
- API endpoints
- RAG integration
- Domain model validation

---

## 8. COORDINATION WITH NOTION AGENT

### 8.1 Schema Change Process

When schema changes are needed:
1. Cursor produces well-structured markdown instructions (like a "Notion schema spec")
2. User pastes those instructions into the Notion agent
3. User confirms schema creation
4. Cursor builds schema detection utilities to check for expected structure

### 8.2 Schema Detection

Cursor should:
- Never assume a new DB or property exists until the user confirms
- Build schema detection utilities to check for expected Notion structure before performing operations

---

## 9. IMPLEMENTATION PHASES

### Phase 1: Domain Model Extraction (Week 17)
- Extract and formalize domain models
- Create Pydantic/dataclass models
- Update ingestion modules to use domain models

**Estimated Effort**: 3-4 days

### Phase 2: Postgres Schema Design (Week 18)
- Design normalized Postgres schema
- Create SQLAlchemy models
- Set up Alembic migrations

**Estimated Effort**: 5-7 days

### Phase 3: FastAPI Service Layer (Week 19-20)
- Implement core API endpoints
- Set up FastAPI application structure
- Add authentication (basic)

**Estimated Effort**: 6-8 days

### Phase 4: Migration Utilities (Week 21)
- Export Notion + SQLite data
- Bootstrap Postgres with initial data
- Implement dual-write capability

**Estimated Effort**: 6-8 days

### Phase 5: RAG Integration (Week 22)
- Update RAG builders for Postgres + Notion hybrid
- Update metadata to use Postgres IDs
- Test RAG queries with new architecture

**Estimated Effort**: 4-5 days

### Phase 6: Transition to Postgres SoT (Week 23)
- Pivot ingestion to Postgres as primary
- Keep Notion as documentation layer
- Final testing and validation

**Estimated Effort**: 1-2 days

**Total Estimated Effort**: 25-34 days (Weeks 17-23)

---

## 10. SUMMARY

### Long-Term Architecture

- **Postgres**: Full relational SoT for all structured data
- **FastAPI**: Service layer for all internal tools & frontends
- **Notion**: Curated narrative + ELN + SOPs + high-level dashboards
- **SQLite**: Optional "chemistry cache" or early-stage prototype DB
- **Pinecone**: Vector store for embeddings & RAG

### Cursor's Responsibilities

1. Capture the domain model
2. Design & implement the Postgres schema and FastAPI
3. Build migration utilities from Notion + SQLite → Postgres
4. Adapt ingestion scripts to target Postgres
5. Keep Notion as a documentation/view layer
6. Maintain RAG integration using Postgres IDs as anchor

---

**Last Updated**: 2025-12-04
**Status**: Future Evolution (After core features stable)

