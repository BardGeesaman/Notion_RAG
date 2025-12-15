# Session Memory

This file stores persistent session memory for the Architect.
It should be updated at natural breakpoints in work sessions to support continuity across machines and over time.

---

## 1. Project Summary (High-Level)

*A short, evolving description of what the project is and is not.*

* **Purpose:** Amprenta Multi-Omics RAG System - A comprehensive knowledge management platform with RAG (Retrieval-Augmented Generation) capabilities for ceramide/sphingolipid neurodegeneration research, extensible to any omics domain.

* **Scope:**
  - Multi-omics data ingestion (lipidomics, metabolomics, proteomics, transcriptomics)
  - Multi-omics signature management with automatic discovery and scoring
  - Semantic search across all data sources with RAG
  - Evidence-based report generation using cross-omics reasoning
  - Pattern discovery and statistical analysis
  - Integration with public repositories (Metabolomics Workbench, GEO, PRIDE, MetaboLights)
  - Chemistry and HTS (High-Throughput Screening) data management

* **Out of Scope:**
  - Real-time data streaming (batch processing only)
  - Direct wet lab integration
  - Clinical trial management
  - Patient data management (research data only)

* **Key Constraints:**
  - **Postgres** as primary system of record (Notion 100% REMOVED per Chairman directive)
  - **Pinecone** for vector storage (semantic index)
  - **Python 3.10+** codebase
  - API rate limits (OpenAI, Pinecone, public repositories)

* **Guiding Principles:**
  - **Idempotency**: All operations safe to re-run
  - **Non-blocking**: Errors don't stop ingestion
  - **Schema resilience**: Graceful property handling
  - **Comprehensive logging**: Clear prefixes ([INGEST], [RAG], etc.)
  - **Postgres as canonical source**: All truth comes from Postgres/SQLite (chemistry)

---

## 2. Current Architecture Overview

*A brief explanation of how the system is currently structured.*

* **Main Components:**
  - **Postgres**: Primary system of record for ALL structured data (Programs, Experiments, Datasets, Features, Signatures, Compounds, HTS, and all relationships).
  - **Pinecone**: Vector index for RAG with OpenAI embeddings.
  - **OpenAI**: Embedding generation (text-embedding-ada-002) + LLM reasoning (GPT-4o/Claude 3.5).
  - **FastAPI**: REST API service layer with 15+ CRUD endpoints.
  - **Streamlit**: Web dashboard (40+ pages) for data exploration and analysis.
  - **Python Scripts**: CLI tools for ingestion and analysis.
  - ~~**SQLite**: Removed (Dec 2025) - Chemistry/HTS migrated to Postgres.~~

* **How They Interact:**
  1. Data files → Ingestion pipelines → Parse and normalize features → **Postgres**
  2. Entities (Experiments, Datasets) → **Postgres** → Establish relationships
  3. Dataset summaries → OpenAI (embeddings) → **Pinecone** (vector storage)
  4. Signatures → Match against dataset features → Score and writeback to **Postgres**
  5. RAG queries → **Pinecone** (semantic search) + **Postgres** (structured search) → OpenAI (reasoning)
  6. **Streamlit Dashboard** provides human-friendly interface and data visualization

* **Known Limitations:**
  - **Jupyter Integration**: Currently lacks deep data exploration/editing capabilities for scientists.
  - **Public Repository Ingestion**: Semi-automated (harvesting implemented, but full ingestion requires manual review).
  - **Visualization**: Needs advanced network graphs (Cytoscape) and genome browsing (IGV).

* **Planned Architectural Changes:**
  - **Jupyter Integration**: Embedded environment for custom analysis.
  - **AWS Cloud Deployment**: ECS, RDS Aurora, ElastiCache.
  - **Advanced Visualization**: Cytoscape.js and IGV.js integration.
  - **Bioinformatics Pipeline**: Nextflow/Snakemake integration.

* **Dependencies to Monitor:**
  - Pinecone index capacity and query performance
  - OpenAI API costs and embedding model changes
  - External mapping services (KEGG, Reactome, UniProt APIs)

---

## 3. Active Tasks

*A living list of work in progress (Synced with ROADMAP.md).*

| Task | Owner (Agent) | Status | Notes |
| :--- | :--- | :--- | :--- |
| **JupyterHub Integration** | Implementor | **✓ Complete** | All 5 phases done (2025-12-15) |
| — Phase 1: API Client Library | Implementor | **✓ Complete** | amprenta-client with 6 resource clients |
| — Phase 2: Write Endpoints | Implementor | **✓ Complete** | Annotation endpoints operational |
| — Phase 3: JupyterHub Deploy | Implementor | **✓ Complete** | DockerSpawner + RDKit stack |
| — Phase 4: SSO Integration | Implementor | **✓ Complete** | JWT TokenAuthenticator |
| — Phase 5: Templates + Launch | Implementor | **✓ Complete** | Templates + context passing |
| **SAR/Voila Test Coverage** | Implementor | **Next Up** | Unit tests for sar_data.py, rgroup.py, notebook_utils.py + API tests for /api/v1/sar/* (3-4 hours) |
| **Innovator-Approved Features** | Implementor | **Pending** | 7 features from INNOVATION_LOG (HTS QC, MOA, Reports, etc.) |
| **AWS Deployment Architecture** | Architect | **Pending** | ECS/RDS/ElastiCache design |
| **Advanced Viz Suite (Cytoscape/IGV)** | Implementor | **Pending** | Network graphs and genome browser |

---

## 4. Completed Tasks (Recent)

*A reverse-chronological log of what has been done recently.*

* [2025-12-15] – **JupyterHub Integration Complete - All 5 Phases Done**:
  - Phase 1: API Client Library (amprenta-client with 6 resource clients) ✓
  - Phase 2: Write Endpoints (annotation endpoints operational) ✓
  - Phase 3: JupyterHub Deployment (DockerSpawner + RDKit stack) ✓
  - Phase 4: SSO Integration (JWT TokenAuthenticator) ✓
  - Phase 5: Templates + Launch (5 templates with context passing) ✓
  - SAR Delta Explorer fully working with best-effort RDKit parsing and Voila rendering
  - Notebook content recovery from Docker container after debugging session
  - Created idempotent seed script (`scripts/seed_sar_data.py`) for reproducibility
  
* [2025-12-15] – **Voila Notebook Suite Audit & Standardization (Complete)**:
  - Audited all 10 Voila notebooks for consistency and robustness
  - Created shared `deploy/jupyterhub/templates/notebook_utils.py` with:
    - RDKit warning suppression (best-effort parsing + log filtering)
    - API client initialization with fallback
    - Demo mode banner for offline/fallback scenarios
    - Safe SMILES parsing utilities
  - Updated all notebooks to use shared utilities (API fallback + demo mode)
  - Updated singleuser Dockerfile to copy all required modules: client/, chemistry/, analysis/pathway/
  - Updated requirements.txt: rdkit-pypi, numpy<2.0, ipywidgets
  - Commits: defa2f8, 0e85518

* [2025-12-15] – **SAR Delta Explorer + JupyterHub/Voila Debugging (Partial)**:
  - JupyterHub “User not allowed” (403) was resolved by ensuring the running hub container loads `/etc/jupyterhub/jupyterhub_config.py` and rebuilding/recreating the `jupyterhub` service so `c.Authenticator.allowed_users={"scientist"}` takes effect.
  - Singleuser image was iterated to include missing runtime pieces for SAR notebooks: `voila`, `ipywidgets`, `amprenta_rag/logging_utils.py`, and chemistry helpers (`amprenta_rag/chemistry/rgroup.py`, `sar_analysis.py`), plus Voila server extension config.
  - **Remaining blocker**: Voila-under-Hub still renders a blank page due to persistent JupyterLab federated-extension `CORE_OUTPUT` shared-module version mismatch (browser console shows unsatisfied @jupyterlab/* versions). Multiple rebuilds/pins did not change CORE_OUTPUT behavior yet.
  - SAR cliffs endpoint still returns empty at defaults (similarity threshold too strict relative to available/selected pair; needs deterministic/unique high-similarity pair in seeded data and/or API behavior adjustment).

* [2025-12-15] – **SAR Delta Explorer + SAR API restored (Complete, portable)**:
  - **Root causes of “hard” debugging**:
    - **Voila rendering** failures were due to **JupyterLab federated-extension shared-module version mismatches** (`_JUPYTERLAB.CORE_OUTPUT` errors in browser console). Symptoms: blank page / widgets never appear even though kernels start.
    - **Cliffs returning empty** was due to **RDKit SMILES parsing/kekulization failures** in seeded aromatics → **0 valid fingerprints** → no pairs, independent of thresholds.
  - **Fixes shipped**:
    - Restored missing SAR API endpoints by adding:
      - `amprenta_rag/api/routers/sar.py`
      - `amprenta_rag/api/services/sar_data.py`
      - and wiring into `amprenta_rag/api/main.py` with prefix `/api/v1/sar`
    - Fixed SAR targets query to use SQLAlchemy `func/distinct/desc` (don’t use `db.func` / `db.text` from a session).
    - Implemented best-effort RDKit parsing in SAR cliffs detection (fallback `sanitize=False` + sanitize minus kekulize) so demo SMILES still fingerprint.
    - Hardened `deploy/jupyterhub/templates/sar_delta_explorer.ipynb` with the same best-effort parsing for grids/MCS inputs.
    - Pinned a known-good Voila stack in `deploy/jupyterhub/singleuser/Dockerfile`: `voila==0.4.6` + `nbclient==0.7.4` to avoid Lab shared-module regressions.
  - **Portable next-session quickstart**:
    - Restored an idempotent SAR seed script: `scripts/seed_sar_data.py` (supports `--reset`, auto-adds repo root to `sys.path` so it runs as `python scripts/seed_sar_data.py` from any workstation).
    - See `docs/VOILA_SAR_QUICKSTART.md` for exact commands to start API, seed SAR data, build singleuser image, and smoke-test Voila on any workstation.
  - **Debugging discipline (memorialized)**:
    - When debugging multi-system issues (Voila/JupyterHub/API/seed data), change **one axis at a time** and run a deterministic smoke test after each change.
    - See `docs/DEBUG_PROTOCOL.md` for the reusable “avoid combinatorial failures” checklist.
* [2025-12-14] – **JupyterHub Phase 3 Complete**: DockerSpawner working locally with RDKit + amprenta-client.
* [2025-12-13] – **JupyterHub Plan Approved**: 5-phase implementation plan documented in ROADMAP.md (commit aff313a).
* [2025-12-13] – **Innovator Agent Created**: New strategic feature ideation agent with charter and communication flow.
* [2025-12-13] – **7 New Features Approved**: HTS QC, Signature Explainability, MOA Inference, Cross-Omics Pathway, Narrative Reports, Data Quality Watcher, Protocol Diff.
* [2025-12-13] – **Code Quality Complete**: All P1-P3 cleanup done (64 db_session fixes, bare excepts, SQLite migration, model splits).
* [2025-12-13] – **Test Coverage Added**: 24 new tests (auth module: 15, models: 9).
* [2025-12-13] – **SQLite Fully Removed**: Chemistry/HTS migrated to PostgreSQL. database.py, schema.py deleted.
* [2025-12-13] – **Operational Maturity**: Achieved 57 passing E2E tests, full dashboard coverage, and implemented optimistic locking for concurrent editing.
* [2025-12-13] – **Feature Completion**: Delivered Genomic Variant Tracking, Scientific Q&A, Advanced UI/UX (Themes, Global Search, Shortcuts), and massive roadmap update (40+ features completed).
* [2025-12-13] – **Roadmap Expansion**: Added Jupyter Integration, AWS Deployment, and Advanced Visualization to roadmap.
* [2025-12-07] – **COMPLETE ROADMAP IMPLEMENTATION**: Delivered all Tiers 1-5 (15+ major features).
* [2025-12-07] – **Notion 100% Removed**: Migrated all modules to Postgres per Chairman directive.
* [2025-12-07] – **Feature Caching System**: 10-100x speedup with DatasetFeatureCache.
* [2025-12-07] – **Batch Ingestion Framework**: 4x speedup with auto-detection.
* [2025-12-07] – **Enhanced Cross-Omics Reasoning**: Disease/matrix/model context.
* [2025-12-07] – **Chemistry & HTS Integration**: SQLite system with 7 tables.
* [2025-12-07] – **Pathway Analysis ID Mapping**: UniProt, KEGG, Reactome integration complete.

*(Older entries archived)*

---

## 5. Decisions Log (Critical Choices)

*A list of explicit decisions made during development.*

### Decision 1: Postgres as Canonical Source of Truth
* **Decision:** Replace Notion with Postgres as the primary system of record.
* **Date:** 2025-12-07
* **Reasoning:** Performance, scalability, and "Chairman's directive" to remove Notion dependence.
* **Impact:** All pipelines write to Postgres. Notion API removed from codebase.

### Decision 2: SQLite Removed - Chemistry/HTS in Postgres
* **Decision:** Migrated chemistry and screening data from SQLite to Postgres.
* **Date:** 2025-12-13 (Migration complete)
* **Original Date:** 2025-12-02 (SQLite introduced)
* **Reasoning:** Unified architecture simplifies operations; Postgres handles HTS volume fine.
* **Impact:** chemistry/database.py and chemistry/schema.py deleted. All data in Postgres.

### Decision 3: Multi-Omics Signature Architecture
* **Decision:** Use feature type inference and modality tagging for cross-omics signatures.
* **Date:** 2025-12-01
* **Reasoning:** Enables single signature to span genes, proteins, metabolites, and lipids.
* **Impact:** Scoring engine handles cross-omics matching natively.

### Decision 4: Optimistic Locking for Concurrency
* **Decision:** Use SQLAlchemy versioning (optimistic locking) for concurrent editing safety.
* **Date:** 2025-12-13
* **Reasoning:** Prevents data overwrites in multi-user environment without complex locking servers.
* **Impact:** Added `version_id` to models and conflict resolution UI in dashboard.

### Decision 5: DockerSpawner for JupyterHub
* **Decision:** Use DockerSpawner with per-user containers instead of LocalProcessSpawner
* **Date:** 2025-12-14
* **Reasoning:** Proper user isolation, no system user requirements, works with DummyAuthenticator
* **Key Learnings:**
  - hub_connect_ip required for container-to-hub connectivity
  - rdkit-pypi (pip) installs in 30 seconds vs conda (60+ minutes on ARM)
  - jupyter/* base images need JUPYTERHUB_SINGLEUSER_APP env var
  - Empty __init__.py files avoid importing entire module trees
  - numpy<2.0 required for rdkit-pypi compatibility
  - External network with external: true for pre-created networks

---

## 6. Open Questions

*Unresolved issues requiring future design or clarification.*

### Question 1: Jupyter Integration Strategy
* **Question:** How to deeply integrate Jupyter? (Embedded vs. Export vs. JupyterHub)
* **Status:** RESOLVED - JupyterHub chosen
* **Decision:** Full JupyterHub deployment with API client library for two-way data flow. 5-phase plan approved (see ROADMAP.md).

### Question 2: Public Repository Automation
* **Question:** Level of automation for public repo ingestion?
* **Status:** Partially Resolved (Harvesting implemented, Manual Import workflow chosen for now).

### Question 3: Multi-Tenant Architecture
* **Question:** How to segregate data for multi-company support?
* **Status:** Planned (Row-level security vs. separate schemas).

---

## 7. Risks & Unknowns

*Potential pitfalls or areas lacking clarity.*

### Risk 1: Pinecone Index Capacity
* **Risk:** Free tier limits vs. growing public dataset collection.
* **Mitigation:** Monitoring vector counts, potential upgrade to paid tier or self-hosted vector DB (pgvector) in AWS phase.

### Risk 2: AWS Deployment Complexity
* **Risk:** Moving from local/single-server to ECS/RDS adds DevOps overhead.
* **Mitigation:** Use IaC (Terraform/Pulumi) and automate via CI/CD.

### Risk 3: Data Quality from Public Repositories
* **Risk:** Inconsistent metadata in GEO/MetaboLights.
* **Mitigation:** Implemented curation workflows and "Draft" status for imports.

---

## 8. Future Ideas / Backlog

*> **NOTE:** See `NEXT_STEPS.md` for the official, prioritized roadmap.*

* **Tier 1 (Immediate)**: Jupyter Integration, Advanced Viz (Cytoscape/IGV).
* **Tier 2 (Strategic)**: AWS Deployment, Bioinformatics Pipeline Orchestrator.
* **Tier 3 (Future)**: Bayesian Inference, GNN for ADMET, Generative Chemistry.

---

## 9. Session Notes

*A running notes section for the Architect to write contextual thoughts.*

### Notes from 2025-12-15

**JUPYTERHUB INTEGRATION COMPLETE - ALL 5 PHASES DONE**

* **Phase 1-5 Complete**: API client library, write endpoints, JupyterHub deployment, SSO, and templates all delivered and working.
* **SAR Delta Explorer Working**: Fixed RDKit SMILES parsing with best-effort sanitization (fallback `sanitize=False` + selective sanitization without kekulization).
* **Voila Rendering Fixed**: Pinned working stack (`voila==0.4.6`, `nbclient==0.7.4`) resolved JupyterLab federated-extension version mismatches.
* **Debugging Session Recovery**: Notebook content was lost during debugging, recovered from Docker container (`jupyter-scientist`), applied RDKit warning suppression.
* **Reproducibility**: Created idempotent seed script (`scripts/seed_sar_data.py`) with `--reset` option for clean data setup.
* **Key Learning**: When debugging complex integrations (JupyterHub + Voila + RDKit), change ONE variable at a time (see memory 12221705).

**VOILA NOTEBOOK AUDIT & STANDARDIZATION**

* **Shared Utilities Created**: `notebook_utils.py` centralizes RDKit suppression, API client initialization, demo mode banner, and safe SMILES parsing.
* **Consistent API Fallback**: All 10 Voila notebooks now gracefully handle API unavailability with demo mode and offline data.
* **Demo Mode Banner**: Visual indicator when notebooks run without API connectivity (useful for demos/testing).
* **Test Coverage Gap Identified**: New SAR functionality (sar_data.py, rgroup.py, notebook_utils.py) and API endpoints (/api/v1/sar/*) need unit tests and API tests.

### Notes from 2025-12-13 (Evening Session)

**JUPYTERHUB PLAN & CODE QUALITY COMPLETION**

*   **JupyterHub Decision**: After evaluating options (JupyterLite, API client only, JupyterHub), Chairman approved full JupyterHub deployment for maximum capability.
*   **5-Phase Plan**: API Client (5d) → Write Endpoints (3d) → JupyterHub Deploy (7d) → SSO (3d) → Templates (2d) = ~3 weeks total.
*   **Innovator Agent**: Created new agent for strategic feature ideation. Reports directly to Chairman, not part of standard workflow.
*   **7 New Features Approved**: HTS QC, Signature Explainability, MOA Inference, Cross-Omics Pathway, Narrative Reports, Data Quality Watcher, Protocol Diff.
*   **Code Quality 100% Complete**: All db_session fixes (64 instances), bare excepts (10), SQLite removal, model splits done.
*   **SQLite Fully Removed**: Chemistry/HTS now in Postgres. Deleted database.py, schema.py, email_cleanup.py.
*   **Test Coverage Added**: 24 new tests for auth and models modules.

### Notes from 2025-12-13 (Earlier)

**OPERATIONAL MATURITY & ROADMAP EXPANSION**

*   **E2E Testing Milestone**: Reached 57 passing Playwright tests covering all core dashboard pages, chemistry flows, RAG, and admin features.
*   **Feature Completeness**: Marked 40+ features as complete including User Feedback, Teams/Projects, In-App Guidance, and Genomic Variant Tracking.
*   **Architecture Update**: Confirmed Notion removal is final. Focus shifted to "Scientist's Cockpit" (Jupyter/Viz) and Infrastructure (AWS).
*   **Governance**: Implemented "Concurrent Editing Safety" via optimistic locking (verified in code, test pending).
*   **Roadmap Refinement**: Clarified status of Multi-Company (Planned) vs Concurrent Editing (Implemented).

### Notes from 2025-12-07

**MASSIVE SESSION - COMPLETED ALL ROADMAP TIERS 1-5**

*   **Multi-Agent System Reactivated**: Full six-agent coordination.
*   **Notion Removal**: 100% complete.
*   **Performance**: Feature caching (100x speedup) and batch ingestion (4x speedup) implemented.
*   **Architecture**: Transitioned to Postgres-centric model.

---

## 10. Continuity Summary (Auto-Generated by Architect)

*To be produced automatically by the Architect at the end of each session.*

**Last Updated:** 2025-12-15

### Summary

The system has reached **production maturity** with **40+ features**, **57 passing E2E tests**, and **fully unified Postgres architecture** (SQLite removed, Notion removed). Code quality optimization is 100% complete.

**JupyterHub integration is COMPLETE** - all 5 phases delivered (2025-12-15). The SAR Delta Explorer is fully working with Voila dashboards operational. All 10 Voila notebooks have been audited and standardized with shared utilities, API fallback, and demo mode.

### Current State

*   **System Status**: Production-Ready. Code quality at 9/10.
*   **Architecture**: Unified Postgres (no SQLite, no Notion), FastAPI, Streamlit (40+ pages), JupyterHub operational.
*   **Test Coverage**: High (24 new auth/model tests, 57 passing E2E tests). Gap identified for new SAR/Voila functionality.
*   **JupyterHub**: All 5 phases complete (API client, write endpoints, deployment, SSO, templates).
*   **Notebook Suite**: 10 Voila notebooks standardized with shared utilities (notebook_utils.py).
*   **Next Focus**: SAR/Voila test coverage (3-4 hours) - Option B for next session.

### JupyterHub Status (ALL COMPLETE)

| Phase | Duration | Status |
|-------|----------|--------|
| 1. API Client Library | 5 days | ✓ Complete |
| 2. Write Endpoints | 3 days | ✓ Complete |
| 3. JupyterHub Deploy | 7 days | ✓ Complete |
| 4. SSO Integration | 3 days | ✓ Complete |
| 5. Templates + Launch | 2 days | ✓ Complete |

### Notebook Audit (COMPLETE)

| Task | Status |
|------|--------|
| 10 Voila notebooks audited | ✓ Complete |
| Shared notebook_utils.py created | ✓ Complete |
| API fallback + demo mode | ✓ Complete |
| Dockerfile updated (modules) | ✓ Complete |
| requirements.txt updated | ✓ Complete |

### Next Session Recommendation

**Option B: Write Tests for SAR + notebook_utils (3-4 hours)**

Unit tests needed:
- `amprenta_rag/api/services/sar_data.py` (cliffs detection, activity grouping)
- `amprenta_rag/chemistry/rgroup.py` (MCS, R-group decomposition)
- `deploy/jupyterhub/templates/notebook_utils.py` (RDKit suppression, safe parsing)

API tests needed:
- `/api/v1/sar/targets` (list targets with activity counts)
- `/api/v1/sar/cliffs` (activity cliffs detection)
- `/api/v1/sar/series` (compound series identification)

Commits this session: defa2f8, 0e85518

### Key Files to Review
*   `docs/ROADMAP.md`: The canonical roadmap (updated this session).
*   `docs/INNOVATION_LOG.md`: Innovator feature tracking.
*   `agents/innovator.md`: New Innovator agent charter.

---

## 11. Resume Instructions (For User)

1. Open `agents/session-memory.md`.
2. Read the **Continuity Summary**.
3. Check `NEXT_STEPS.md` for the latest task status.
4. In your IDE/agent environment, say:

   ```text
   Architect:
   Please rehydrate your context from agents/session-memory.md and generate a fresh plan for next steps.
   ```
