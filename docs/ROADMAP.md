# ROADMAP (Single Source of Truth)

**Last Updated**: 2025-12 (December 2025)

Simple status legend:
- ✅ DONE
- ⏳ IN PROGRESS / NEXT UP
- ❌ NOT STARTED / FUTURE

---

## ✅ Current State Summary (December 2025)

- ✅ **PostgreSQL is the sole system of record** (no Notion, no SQLite runtime)
- ✅ **Notion removed completely**
- ✅ **Streamlit dashboard production-ready** (40+ pages, modular core/chemistry splits, auth sessions)
- ✅ **FastAPI REST API operational**
- ✅ **Chemistry & HTS integrated with PostgreSQL**
- ✅ **Public repository ingestion complete** (GEO, PRIDE, MetaboLights, Metabolomics Workbench)
- ✅ **RAG operational** (Pinecone-backed retrieval + LLM synthesis)

---

## ✅ DONE (Key Capabilities)

### Platform & Architecture
- ✅ Postgres-only architecture (Programs/Experiments/Datasets/Features/Signatures/Compounds/HTS)
- ✅ Notion removed (no runtime dependencies)
- ✅ Streamlit dashboard: modular page architecture + sessions/auth
- ✅ FastAPI service layer operational

### RAG Maturity
- ✅ Hybrid retrieval (dense + structured search patterns)
- ✅ Citations/source attribution in RAG answers
- ✅ Evaluation and quality controls (per `NEXT_STEPS.md`)

### Chemistry & HTS
- ✅ Compound registration (dedupe, properties, corporate IDs)
- ✅ Structure search + SAR tooling
- ✅ HTS campaign ingestion + querying (Postgres-backed)

### Discovery / Repositories
- ✅ Automated discovery workflow + scheduled harvesting
- ✅ GEO/PRIDE/MetaboLights/MW ingestion

---

## ⏳ NEXT UP (Prioritized)

### Approved Features (2025-12-13)
- ⏳ HTS QC & Triage Assistant (3-5 days)
- ⏳ Signature Match Explainability (2-4 days)
- ⏳ MOA Inference from HTS + Multi-Omics (5-8 days)
- ⏳ Cross-Omics Pathway Analysis (5-8 days)
- ⏳ One-Click Narrative Reports (2-3 days)
- ⏳ Data Quality Watcher (2-3 days)
- ⏳ Protocol Version Diff & Deviations Audit (2-4 days)

### 1) JupyterHub Integration (TOP PRIORITY - ~3 weeks)

*Note: Requires AWS infrastructure (see Section 2).*

**Phase 1 — API Client Library (5 days)**:
- Create `amprenta-client` Python package
- Typed wrappers for all CRUD endpoints
- Auth handling (API keys/JWT)
- pip-installable

**Phase 2 — Write Endpoints (3 days, parallel with Phase 1)**:
- POST /datasets/{id}/annotations
- POST /experiments/{id}/analysis_results
- POST /signatures/{id}/annotations
- POST /compounds/{id}/annotations

**Phase 3 — JupyterHub Deployment (7 days)**:
- ✅ Docker/K8s deployment (local DockerSpawner validated)
- ✅ Persistent user workspaces
- ✅ Pre-installed packages (pandas, seaborn, rdkit, amprenta-client)

**Phase 4 — SSO Integration (3 days)**:
- Shared auth with Streamlit
- Token/key management

**Phase 5 — Templates + Launch (2 days)**:
- "Open in Jupyter" button on Dataset/Experiment pages
- Notebook templates per omics type

### 2) AWS Deployment / Infrastructure Hardening
- ❌ IaC (Terraform/Pulumi)
- ❌ ECS/RDS/ElastiCache deployment architecture
- ❌ CI/CD pipeline (GitHub Actions → deploy)

### 3) Advanced Visualization Suite
- ❌ Cytoscape.js networks (relationships/lineage)
- ❌ IGV.js genome/variant browser
- ❌ Ag-Grid for large-table browsing/editing

### 4) Testing Closure (Implemented-but-Untested)
- ⏳ Playwright E2E: **Concurrent Editing / optimistic locking** flow

### 5) Multi-Company Support (Multi-Tenancy)
- ❌ Company model + admin roles
- ❌ Data segregation (RLS vs schema separation)
- ❌ Company-scoped settings/branding

---

## ❌ FUTURE (Backlog)

- ❌ Bioinformatics pipeline orchestrator (Nextflow/Snakemake)
- ❌ Data Version Control (DVC)
- ❌ **Pinecone → pgvector migration** - Self-host vector search in Postgres for reduced costs and latency
- ❌ Bayesian inference & optimization workflows
- ❌ ML/AI model registry + predictive ADMET
- ❌ OOP Refactoring Review

---

## Notes

- The full, detailed product checklist remains in `NEXT_STEPS.md`.
- Historical roadmap documents are archived under `md_archive/roadmaps/`.

