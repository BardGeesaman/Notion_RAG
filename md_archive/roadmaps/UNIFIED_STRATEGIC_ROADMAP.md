# Unified Strategic Roadmap - Amprenta Multi-Omics Platform (Postgres-Only)

**Last Updated**: 2025-12 (December 2025)

This roadmap is the **single source of truth** for what is **done**, what is **next**, and what remains **future work**.

---

## ✅ Current State (December 2025)

### **Architecture**
- ✅ **PostgreSQL is the sole system of record** for all entities (Programs, Experiments, Datasets, Features, Signatures, Compounds, HTS).
- ✅ **Notion has been completely removed** (no ingestion, storage, or UI dependencies).
- ✅ **No SQLite in the runtime data path** (historical/prototype only).
- ✅ **Streamlit Dashboard is production-ready**:
  - Modular page architecture (core + chemistry splits)
  - Auth system with sessions
  - 40+ working pages
- ✅ **FastAPI REST API is operational** (CRUD + auth/permissions where applicable).
- ✅ **Chemistry & HTS are integrated with PostgreSQL**.
- ✅ **Public repository ingestion is complete** (GEO, PRIDE, MetaboLights, Metabolomics Workbench).

### **Canonical Roadmap References**
- `NEXT_STEPS.md` (product roadmap / feature checklist)
- `docs/USER_GUIDE.md` (UI user guide)
- `docs/NARRATIVE_WALKTHROUGH.md` (contextual narrative walkthrough)

---

## ✅ DONE (Major Milestones)

### **Platform & Architecture**
- ✅ DONE: Postgres-only architecture (all core entities)
- ✅ DONE: Notion removal
- ✅ DONE: FastAPI operational
- ✅ DONE: Streamlit dashboard production-ready and modularized

### **Data Ingestion & Operations**
- ✅ DONE: Batch ingestion framework + operational tooling
- ✅ DONE: Repository harvesting + manual review/import workflow
- ✅ DONE: Public repository ingestion (GEO/PRIDE/MetaboLights/MW)

### **Chemistry & HTS**
- ✅ DONE: Compound registration + dedupe + property calculation
- ✅ DONE: HTS campaign ingestion + querying
- ✅ DONE: Postgres-backed chemistry + HTS integration

### **RAG & Retrieval**
- ✅ DONE: Pinecone-backed semantic retrieval + LLM synthesis
- ✅ DONE: Hybrid retrieval and citation support (see `NEXT_STEPS.md` for the detailed checklist)

---

## 🎯 NEXT UP (Realistic Priorities)

### **1) Jupyter Notebook Integration (Highest Priority)**
**Goal**: scientists can launch a notebook from a dataset for custom analysis and visualization.

- Deliverables:
  - Launch notebook from any dataset (read-only first)
  - Notebook templates for common workflows
  - Optional controlled write-back in a later phase

### **2) Infrastructure Hardening / AWS Deployment**
**Goal**: reproducible, scalable deployment.

- Deliverables:
  - IaC (Terraform/Pulumi)
  - ECS/RDS/ElastiCache reference architecture
  - CI/CD pipeline

### **3) Advanced Visualization Suite**
**Goal**: stronger interactivity and performance on large datasets.

- Deliverables:
  - Cytoscape.js networks (relationships/lineage)
  - IGV.js variant browser
  - Ag-Grid tables for large-scale browsing/editing

### **4) E2E Coverage Closure (Implemented-but-Untested)**
**Goal**: eliminate the remaining E2E gaps.

- Deliverables:
  - Add Playwright E2E for **Concurrent Editing / optimistic locking** flow

### **5) Multi-Company Support (Multi-Tenancy)**
**Goal**: company-level segregation and admin roles.

- Deliverables:
  - Company model + per-company settings
  - Data segregation strategy (RLS or schema separation)
  - UI/admin role model

---

## 🧭 Future (Backlog)

- Bayesian inference & optimization workflows
- ML/AI model registry + predictive ADMET
- Bioinformatics pipeline orchestrator (Nextflow/Snakemake)
- Data Version Control (DVC)

---

## Historical Note (for context)

Earlier prototypes referenced Notion/SQLite to accelerate iteration. These are **no longer used** in the production architecture.

