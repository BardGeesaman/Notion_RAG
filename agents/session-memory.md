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
  - **Postgres**: Primary system of record for all structured data (Programs, Experiments, Datasets, Features, Signatures, and all relationships).
  - **SQLite**: Chemistry and screening data (7 tables: compounds, libraries, campaigns, results).
  - **Pinecone**: Vector index for RAG with OpenAI embeddings.
  - **OpenAI**: Embedding generation (text-embedding-ada-002) + LLM reasoning (GPT-4o/Claude 3.5).
  - **FastAPI**: REST API service layer.
  - **Streamlit**: Web dashboard (30+ pages) for data exploration and analysis.
  - **Python Scripts**: CLI tools for ingestion and analysis.

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

*A living list of work in progress (Synced with NEXT_STEPS.md).*

| Task | Owner (Agent) | Status | Notes |
| :--- | :--- | :--- | :--- |
| **Jupyter Integration Design** | Architect | **Pending** | Highest priority: Data Viz & Editing |
| **AWS Deployment Architecture** | Architect | **Pending** | ECS/RDS/ElastiCache design |
| **Advanced Viz Suite (Cytoscape/IGV)** | Implementor | **Pending** | Network graphs and genome browser |
| **E2E Test: Concurrent Editing** | Tester | **Pending** | Gap identified in roadmap |
| **Public Repo Expansion** | Architect | **Planned** | Add S2/OpenAlex integration |

---

## 4. Completed Tasks (Recent)

*A reverse-chronological log of what has been done recently.*

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

### Decision 2: SQLite for Chemistry/HTS Data
* **Decision:** Use SQLite for chemistry and screening data (high volume).
* **Date:** 2025-12-02
* **Reasoning:** HTS campaigns can have 100k-1M molecules; relational structure needed but simpler than full Postgres for prototype.
* **Impact:** Chemistry data stays in SQLite (now migrating to Postgres for unification).

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

---

## 6. Open Questions

*Unresolved issues requiring future design or clarification.*

### Question 1: Jupyter Integration Strategy
* **Question:** How to deeply integrate Jupyter? (Embedded vs. Export vs. JupyterHub)
* **Status:** Open (High Priority)
* **Proposed:** Embedded JupyterLite or Kernel Gateway for tight integration.

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

### Notes from 2025-12-13

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

**Last Updated:** 2025-12-13

### Summary

The system has reached a high level of operational maturity with **40+ features completed** and **57 passing E2E tests**. The architecture is now fully **Postgres-centric** (Notion removed), with a robust Streamlit dashboard and RAG capabilities. Key recent additions include **Genomic Variant Tracking**, **Scientific Q&A**, and **Optimistic Locking** for concurrent editing.

The roadmap has been updated to prioritize **Jupyter Notebook Integration** (for deep analysis), **Advanced Visualization** (Cytoscape/IGV), and **AWS Cloud Deployment** (for scalability). Multi-company support and automated reporting are planned but pending UI implementation.

### Current State

*   **System Status**: Production-Ready (Core & RAG), expanding into Advanced Analytics.
*   **Test Coverage**: High (Core flows covered), gap identified in Concurrent Editing testing.
*   **Architecture**: Stable Postgres/FastAPI/Streamlit stack.
*   **Next Focus**: Integration (Jupyter, AWS) and Visualization.

### Recommended Next Steps

1.  **Design Jupyter Integration**: Decide on architecture (embedded vs. hub) and begin implementation.
2.  **AWS Architecture**: Design ECS/RDS infrastructure and create IaC scripts.
3.  **Advanced Visualization**: Integrate Cytoscape.js for network graphs.
4.  **Test Coverage**: Add E2E test for Concurrent Editing (Optimistic Locking).

### Key Files to Review
*   `NEXT_STEPS.md`: The canonical roadmap.
*   `amprenta_rag/tests/dashboard/test_e2e_platform.py`: The E2E test suite.
*   `amprenta_rag/analysis/`: Analytics modules (enrichment, stats).

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
