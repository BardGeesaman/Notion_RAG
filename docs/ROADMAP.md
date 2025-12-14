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

### 1) Jupyter Notebook Integration (Highest Priority)

**Phase 1 — Export to Notebook (Current)**:
- ⏳ **Export to Notebook** button on dataset pages
- ⏳ Generate a `.ipynb` pre-loaded with dataset data
- ⏳ Templates per omics type
- ⏳ Works with any local Jupyter

**Phase 2 — JupyterHub (Future, requires AWS)**:
- ❌ Deploy JupyterHub for multi-user support
- ❌ Persistent workspaces per scientist
- ❌ Write-back results to Postgres

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

---

## Notes

- The full, detailed product checklist remains in `NEXT_STEPS.md`.
- Historical roadmap documents are archived under `md_archive/roadmaps/`.

