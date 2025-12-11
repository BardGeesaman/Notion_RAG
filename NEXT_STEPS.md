# Next Steps - Implementation Roadmap

**Last Updated**: 2025-12-10  
**Status**: Pathway analysis complete, Postgres-only architecture

---

## ✅ COMPLETED FEATURES

### 1. Pathway Analysis ID Mapping ✅ COMPLETE

**Status**: Fully implemented and tested (December 2025)

**Implemented Features**:
- ✅ UniProt mapping for proteins (`map_protein_to_uniprot`)
- ✅ Gene symbol → KEGG ID mapping (`map_gene_to_kegg`)
- ✅ Metabolite name → KEGG compound ID mapping (`map_metabolite_to_kegg`)
- ✅ Protein → KEGG mapping via UniProt (`map_protein_to_kegg`)
- ✅ Gene/Protein → Reactome mapping (`map_gene_to_reactome`, `map_protein_to_reactome`)
- ✅ Batch mapping utility (`batch_map_features_to_pathway_ids`)
- ✅ Caching and rate limiting for API calls
- ✅ Full KEGG/Reactome API integration
- ✅ Fisher's exact test with Benjamini-Hochberg FDR correction
- ✅ Dashboard UI integration (Analysis Tools page)
- ✅ Unit tests (47 tests covering ID mapping and enrichment)

**Key Files**:
- `amprenta_rag/analysis/id_mapping.py` - ID mapping services (449 lines)
- `amprenta_rag/analysis/pathway/mapping.py` - Pathway mapping (360 lines)
- `amprenta_rag/analysis/pathway/enrichment.py` - Enrichment analysis (352 lines)
- `scripts/dashboard/pages/analysis.py` - Dashboard UI (480 lines)
- `amprenta_rag/tests/analysis/test_id_mapping.py` - ID mapping tests (25 tests)
- `amprenta_rag/tests/analysis/test_pathway_enrichment.py` - Enrichment tests (22 tests)

**Usage**:
- Dashboard: Analysis Tools → Pathway Enrichment
- CLI: `python scripts/test_pathway_enrichment.py`

---

### 2. Chemistry & HTS Integration ✅ COMPLETE

**Status**: Fully implemented and tested (December 2025)

**Implemented Features**:
- ✅ SQLite chemistry database (compounds, HTS campaigns, biochemical results)
- ✅ SMILES normalization with RDKit fallback
- ✅ Screening ingestion pipeline
- ✅ Compound-signature linking (reverse matching)
- ✅ RAG integration for chemistry queries
- ✅ RAG Integration: Compound embeddings to Pinecone (`upsert_compound_to_pinecone`, `query_compounds_vector`, `batch_upsert_compounds`)
- ✅ Dashboard chemistry page with Signature Links tab
- ✅ Unit tests for compound linking

**Key Files**:
- `amprenta_rag/chemistry/compound_linking.py`
- `amprenta_rag/query/chemistry_query.py`
- `amprenta_rag/tests/chemistry/test_compound_linking.py`
- `scripts/dashboard/pages/chemistry.py`

### 3. Experimental Design Metadata System ✅ COMPLETE

**Status**: Fully implemented (December 2025)

**Implemented Features**:
- ✅ Schema fields added (`design_type`, `sample_groups`, `timepoints`, etc.)
- ✅ Alembic migration created
- ✅ `design_extraction.py` for GEO/MW pattern detection
- ✅ `design_aware_stats.py` for design-specific analysis

**Key Files**:
- `amprenta_rag/ingestion/design_extraction.py`
- `amprenta_rag/analysis/design_aware_stats.py`

---

### 4. UX Polish & Enhancements ✅ COMPLETE

**Status**: Fully implemented (December 2025)

**Implemented Features**:
- ✅ **UX Polish**: Error Messages (`errors.py`), Progress Indicators (tqdm), Batch Summaries, Configuration Validation (`config_check.py`)
- ✅ **PDF Export**: Report generation (`pdf_export.py`), Dashboard/CLI integration
- ✅ **Visualization**: Viz helpers (`viz_helpers.py`), Pathway Network Visualization (force-directed)
- ✅ **Additional Visualizations**: Program-Signature heatmap, Dataset Similarity dendrogram, Pathway Enrichment bars, Cross-Omics Convergence
- ✅ **Performance Utilities**: Parallel processing helpers (`parallel_map`, `chunked_parallel`)

### 5. Advanced Analytics ✅ COMPLETE

**Status**: Fully implemented (December 2025)

**Implemented Features**:
- ✅ Differential Expression Analysis (`differential_expression.py`)
- ✅ Time-series Analysis (`timeseries.py`)
- ✅ Cohort Comparison (`cohort_comparison.py`)

### 6. ArrayExpress Repository ✅ COMPLETE

**Status**: Fully implemented (December 2025)

**Implemented Features**:
- ✅ New repository integration for transcriptomics data

### 7. ELN Experiment Type Integration ✅ COMPLETE
- ✅ Edit Design tab for manual design_type configuration
- ✅ Auto-detect Design Types batch operation  
- ✅ Import from Repository (GEO, MW, PRIDE)
- ✅ Auto-creates experiments with design_type detection
- ✅ Verified working with real GEO study (GSE153873 -> case_control)

### 8. Multi-User Authentication ✅ COMPLETE
- ✅ User model with bcrypt password hashing
- ✅ Login page with session management
- ✅ Role-based access (admin, researcher, viewer)
- ✅ Admin-only user registration page
- ✅ Logout functionality in sidebar
- ✅ DISABLE_AUTH env flag for testing

### 9. Audit & Data Ownership ✅ COMPLETE
- ✅ AuditLog model tracking user actions
- ✅ Audit log viewer dashboard (admin-only, filters, CSV export)
- ✅ Auto-logging on login/logout
- ✅ created_by_id FK on Program, Experiment, Dataset, Signature

### 11. RAG Hybrid Search ✅ COMPLETE
- ✅ tsvector search column on RAGChunk with GIN index
- ✅ BM25 full-text search via PostgreSQL  
- ✅ Reciprocal Rank Fusion (RRF) for result merging
- ✅ use_hybrid and hybrid_alpha params in query_rag()

### 12. Automated Discovery Workflow ✅ COMPLETE
- ✅ DiscoveryJob/DiscoveredStudy models
- ✅ Repository scanning service (GEO, MetaboLights)
- ✅ Dashboard: Run scans, review/import studies, job history

### 13. Chemical Registration System ✅ COMPLETE
- ✅ Corporate ID generation (AMP-XXXXX)
- ✅ Duplicate checking via SMILES/InChIKey
- ✅ PostgreSQL-backed compound storage (migrated from SQLite)
- ✅ Dashboard registration UI
- ✅ Playwright E2E test
- ✅ Structure Search (substructure + similarity) with RDKit
- ✅ SAR Analysis with Lipinski Ro5 compliance
- ✅ SAR analysis, structure search, compound linking all use PostgreSQL
- ✅ SQLite chemistry/database.py deprecated

### 14. Lead Optimization Data Models ✅ COMPLETE
- ✅ BiochemicalAssay and ActivityResult models
- ✅ Activity cliff detection algorithm
- ✅ Activity Cliffs Network visualization (circular graph)
- ✅ ADMEResult model (permeability, stability, CYP inhibition)
- ✅ PKStudy model (AUC, Cmax, half-life, bioavailability)
- ✅ ToxicologyResult model (hERG, Ames, cytotoxicity)
- ✅ Lead Optimization dashboard tab with ADME/PK/Tox entry forms
- ✅ Compound profile view with liability alerts
- ✅ Playwright error detection utility added

### 15. Scientific Q&A Tracker ✅ COMPLETE
- ✅ SavedQuestion/SavedAnswer models with versioning
- ✅ Dashboard page with Ask, Browse, Re-run, Export tabs
- ✅ Persistent RAG answer storage

---

## 🎯 IMMEDIATE NEXT STEPS (Priority Order)

### 1. Integration Enhancements

**Public Repositories**:
- [ ] Add more repository types (others)
- [x] Automated discovery workflows
- [ ] Scheduled harvesting
- [ ] **Publication & Supplementary Data Extraction**:
  - [ ] Extract experiment details from PDF publications
  - [ ] Parse supplementary tables (Excel/CSV) for data
  - [ ] Link publications to repository datasets

**Security & Collaboration**:
- [ ] Team-based data access
- [ ] **User Feedback System**:
  - [ ] Error Reporting (Automatic capture + user description)
  - [ ] Feature Request Tracker (Vote/Prioritize)
  - [ ] Admin Dashboard for triage
- [ ] **In-App User Guidance**:
  - [ ] Interactive Feature Tours (e.g., "How to create a protocol")
  - [ ] Contextual Help/Tooltips (hover for definitions)
  - [ ] Embedded Documentation Widget (searchable guide)
  - [ ] **AI Help Assistant** (Chatbot for "How do I...?" questions)
- [ ] **Navigation & UI Organization**:
  - [ ] Functional Grouping (Discovery, ELN, Analysis, Admin)
  - [ ] Collapsible Sidebar Sections
  - [ ] Favorites/Recent items list

### 2. ELN & Workflow Features

**Sample Inventory**:
- [x] Physical location tracking (freezers/boxes)
- [x] Sample lineage and parent/child tracking
- [x] Barcode generation support

**Compliance & Integrity**:
- [ ] 21 CFR Part 11 readiness (audit trails)
- [ ] Electronic signatures (future)
- [ ] Data immutability logs

**Scientific Q&A & Insight Tracker**:
- [x] Question Registry (CRUD for scientific questions)
- [x] Persistent RAG Answers (Save outputs + evidence)
- [x] Versioning (Re-run questions as data updates)
- [x] Export capabilities (PDF/CSV of Q&A reports)

### 3. Chemistry Informatics Enhancements

**Chemical Registration**:
- [x] Corporate ID generation (e.g. AMP-001)
- [x] Salt and batch handling
- [x] Duplicate checking logic

**Cheminformatics Search**:
- [x] Substructure search (SMARTS)
- [x] Similarity search (Tanimoto)
- [ ] Pharmacophore search

**SAR Analysis**:
- [x] Calculated properties (LogP, MW, TPSA)
- [x] Activity cliffs visualization
- [ ] R-group decomposition

**Lead Optimization Data Models**:
- [x] ADME/DMPK Assays (Stability, Permeability, CYP inhibition)
- [x] In Vivo PK Studies (AUC, Cmax, Tmax, Bioavailability)
- [x] Safety & Toxicology (hERG, Ames, Cytotoxicity)

**Structural Biology & Virtual Screening**:
- [ ] Protein Structure Visualization (Mol* / PyMOL integration)
- [ ] Virtual Docking Pipeline (AutoDock Vina / Smina)
- [ ] Binding Site Analysis & Pocket Detection
- [ ] Ligand-Protein Interaction Profiler

### 4. RAG Maturity & Enhancements

**Retrieval Optimization**:
- [x] Hybrid Search (Sparse/BM25 + Dense Vectors)
- [ ] Reranking (Cross-encoder re-ranking)
- [ ] Query Expansion/Transformation (HyDE)

**Quality & Attribution**:
- [ ] Citation/Source Attribution
- [ ] Evaluation Framework (RAGAS/TruLens)
- [ ] Hallucination detection

**Advanced Capabilities**:
- [ ] Semantic Caching
- [ ] Agentic RAG (multi-step reasoning)

### 5. Future Innovations

**Experimental Design Assistant**:
- [ ] Power Analysis for sample size estimation
- [ ] Design recommendation engine
- [ ] Confounding variable detection

**Advanced Analytics**:
- [ ] Multi-dataset pattern detection

---

## 🚀 QUICK WINS (Can Do Now)

*All identified quick wins (PDF Export, Visualization Helpers, Progress Bars, Error Messages) have been completed.*

---

## 📋 RECOMMENDED WORKFLOW

### Phase 1: Integration Enhancements (Next)
1. Public Repository Expansion
2. Security & Multi-User Support (Team Access)

### Phase 2: ELN & Workflow Features
1. Scientific Q&A & Insight Tracker
2. Sample Inventory System
3. Audit Trails & Compliance

### Phase 3: Chemistry Informatics
1. Lead Optimization Data Models (ADME/PK)
2. Chemical Registration System
3. Advanced Structure Search
4. Structural Biology Tools (Docking/Viz)
5. SAR Analysis Tools

### Phase 4: RAG Maturity
1. Hybrid Search & Reranking
2. Citation & Attribution System
3. Evaluation Pipeline

### Phase 5: Future Innovations
1. Multi-dataset Pattern Detection
2. Experimental Design Assistant & Power Analysis

---

## 🎯 IMMEDIATE ACTION ITEMS

**This Week**:
1. Review Public Repository Expansion plan
2. Start ArrayExpress integration research

---

## 💡 SUGGESTIONS

**For Maximum Impact**:
- Focus on Advanced Analytics (enables deeper biological insights)
- Then Integration Enhancements (expands data sources and security)

**For Quick Wins**:
- Automated Discovery Workflows (Public repos)

---

## 📞 QUESTIONS TO CONSIDER

1. **Priority**: Which analytic capability is most urgent?
   - Differential Expression?
   - Time-series Analysis?
   - Cohort Comparison?

2. **Resources**: Do you have access to:
   - Longitudinal datasets for time-series testing?
   - Large cohorts for comparison testing?
   - Feedback on report formats?

3. **Timeline**: What's the target for:
   - Multi-user rollout?
   - Public release?

---

**Current System Status**: ✅ Production Ready (Full ELN System Complete)  
**Architecture**: Postgres-only (Notion removed December 2025)  
**Next Milestone**: Chemistry Informatics & RAG Enhancements  
**Estimated Time**: 3-4 weeks
