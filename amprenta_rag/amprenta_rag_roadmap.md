# Amprenta Platform – Engineering Roadmap (Postgres-Only)

**Last Updated**: 2025-12 (December 2025)

This file tracks the **engineering state** and **next steps** for the Amprenta platform.

---

## 1. Goals

- Support Amprenta’s research programs with a **high-quality, explainable platform** for multi-omics + chemistry + RAG.
- Keep the codebase **modular, testable, and comprehensible**.
- Maintain a clear separation between:
  - **System of record**: PostgreSQL
  - **Semantic retrieval**: Pinecone
  - **UI**: Streamlit dashboard
  - **API**: FastAPI

---

## 2. Current State (December 2025)

### ✅ Data Layer
- ✅ PostgreSQL is the **sole system of record** (no Notion, no SQLite runtime).
- ✅ SQLAlchemy models + Alembic migrations.

### ✅ API Layer (FastAPI)
- ✅ Operational REST API (core CRUD).
- ✅ Authentication supported.

### ✅ UI Layer (Streamlit)
- ✅ Production-ready dashboard with **40+ pages**.
- ✅ Modular page organization (core + chemistry splits).
- ✅ Auth sessions.

### ✅ Chemistry & HTS
- ✅ Integrated in PostgreSQL (compound registration, screening ingestion, SAR workflows).

### ✅ Public Repository Ingestion
- ✅ GEO (transcriptomics)
- ✅ PRIDE (proteomics)
- ✅ MetaboLights (metabolomics)
- ✅ Metabolomics Workbench (metabolomics/lipidomics)

### ✅ RAG
- ✅ Pinecone-backed retrieval + LLM synthesis.
- ✅ RAG UI pages in dashboard.

---

## 3. DONE (Engineering Milestones)

- ✅ DONE: Postgres-only architecture (Notion removed)
- ✅ DONE: Streamlit dashboard modularization + production readiness
- ✅ DONE: Auth system + sessions
- ✅ DONE: FastAPI operational
- ✅ DONE: Chemistry & HTS integrated in Postgres
- ✅ DONE: Public repository ingestion complete (GEO/PRIDE/MetaboLights/MW)

---

## 4. NEXT UP (Engineering Priorities)

### 4.1 Jupyter Integration (Highest Priority)
- Embed notebook workflows for dataset exploration.
- Start read-only; add write-back only after safety model is clear.

### 4.2 Deployment / AWS Hardening
- ECS/RDS/ElastiCache reference architecture
- Terraform/Pulumi IaC
- CI/CD pipeline

### 4.3 Visualization Upgrades
- Cytoscape.js for relationship graphs
- IGV.js for variant visualization
- Ag-Grid for large tables

### 4.4 Testing Closure
- Add E2E test for **Concurrent Editing / optimistic locking** flow.

---

## 5. Pending / Future Work

- Multi-company support (multi-tenancy)
- Nextflow/Snakemake orchestrator
- DVC integration
- Bayesian inference / ML model registry

---

## Historical Note

Older iterations referenced Notion and SQLite prototypes. These are **deprecated** and not part of the production architecture.

