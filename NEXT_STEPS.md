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
- ✅ **Data Export Suite**: CSV/JSON/Excel export for all entity types
- ✅ **Data Import Suite**: CSV/JSON bulk import with validation
- ✅ **Dashboard Themes**: Light/Dark mode support (via Streamlit config)
- ✅ **Notifications System**: In-app alerts for job completions (toasts/status indicators)
- ✅ **Activity Dashboard**: Recent activity feed, stats, quick actions
- ✅ **Dashboard Widgets**: Metric cards on overview page (experiments, compounds, samples, discoveries)
- ✅ **System Health Dashboard**: DB stats, API status, system metrics
- ✅ **Keyboard Shortcuts**: Vim-style navigation and global hotkeys
- ✅ **Global Search**: Unified cross-entity search ("Spotlight" style)
- ✅ **Saved Filters**: Save and reuse filter presets for experiments/compounds
- ✅ **Bookmarks System**: Pin items for quick access (sidebar favorites)
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
- ✅ Session Timeout (auto-logout security)
- ✅ API Rate Limiting (abuse protection)
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
- ✅ Cross-encoder reranking (ms-marco-MiniLM model)
- ✅ Semantic Cache (embedding similarity, 1hr TTL, 0.92 threshold)
- ✅ HyDE (Hypothetical Document Embeddings) query expansion
- ✅ Hallucination Detection (Groundedness scoring, unsupported claims)
- ✅ RAGAS Evaluation Framework (Faithfulness, Relevance, Context Precision)
- ✅ use_hybrid, hybrid_alpha, use_rerank, use_cache, use_hyde, check_hallucination, evaluate params

### 12. Automated Discovery Workflow ✅ COMPLETE
- ✅ DiscoveryJob/DiscoveredStudy models
- ✅ HarvestSchedule model for automated scans (APScheduler)
- ✅ Repository scanning service (GEO, MetaboLights)
- ✅ Dashboard: Run scans, review/import studies, job history, schedules
- ✅ Background polling with manual review workflow

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

### 16. Collaboration & Teams ✅ COMPLETE
- ✅ Team, TeamMember, Project models
- ✅ Role-based permissions (owner/admin/member/viewer)
- ✅ Teams & Projects dashboard page

### 17. User Feedback System ✅ COMPLETE
- ✅ ErrorReporting (auto-stack trace) and UserFeedback models
- ✅ Feedback submission widget
- ✅ Admin triage dashboard with vote tracking

### 18. In-App User Guidance ✅ COMPLETE
- ✅ Contextual Help/Tooltips (hover for definitions)
- ✅ Embedded Documentation Widget (searchable guide)
- ✅ Interactive Feature Tours (via driver.js integration)

### 19. Workflow Automation Engine ✅ COMPLETE
- ✅ WorkflowRule and WorkflowExecution models
- ✅ Trigger types: experiment_created, compound_registered, discovery_imported, sample_transferred
- ✅ Action types: send_notification, add_note, run_validation
- ✅ Dashboard page for creating/managing rules and viewing execution history

---

## 🎯 IMMEDIATE NEXT STEPS (Priority Order)

### 1. Integration Enhancements

**Public Repositories**:
- [ ] Add more repository types (others)
- [ ] **Integrate Semantic Scholar / OpenAlex**:
  - [ ] Use Semantic Scholar API for citation graphs & AI relevance
  - [ ] Use OpenAlex API for author/institution metadata
  - [ ] Power "Automated Literature Analysis" with these sources
- [x] Automated discovery workflows
- [x] Scheduled harvesting
- [ ] **Publication & Supplementary Data Extraction**:
  - [ ] Extract experiment details from PDF publications
  - [ ] Parse supplementary tables (Excel/CSV) for data
  - [ ] Link publications to repository datasets

**Security & Collaboration**:
- [x] Team-based data access
- [ ] **Multi-Company Support (Multi-Tenancy)**:
  - [ ] Company Model & Super Admin Role
  - [ ] Data Segregation (Row-level Security / Separate Schemas)
  - [ ] Company-specific Settings & Branding
- [x] **Feature Visibility Controls**:
  - [x] Admin UI to toggle features per user/role
  - [x] Role-based UI customization (e.g., Hide "Chemistry" for Biologists)
  - [x] Granular permission sets for menu items
- [x] **User Feedback System**:
  - [x] Error Reporting (Automatic capture + user description)
  - [x] Feature Request Tracker (Vote/Prioritize)
  - [x] Admin Dashboard for triage
- [x] **In-App User Guidance**:
  - [x] Interactive Feature Tours (e.g., "How to create a protocol")
  - [x] Contextual Help/Tooltips (hover for definitions)
  - [x] Embedded Documentation Widget (searchable guide)
  - [x] **AI Help Assistant** (Chatbot for "How do I...?" questions)
- [ ] **Navigation & UI Organization**:
  - [ ] Functional Grouping (Discovery, ELN, Analysis, Admin)
  - [ ] Collapsible Sidebar Sections
  - [ ] Favorites/Recent items list
- [ ] **Advanced UI/UX Enhancements**:
  - [ ] "Scientist's Cockpit" Dashboard (Widgets for recent data, alerts, tasks)
  - [ ] Global "Spotlight" Search (Cmd+K for everything)
  - [ ] Advanced Data Grids (Ag-Grid: filtering, sorting, pinning)
  - [ ] Theme Density Control (Compact vs. Comfortable)

### 2. ELN & Workflow Features

**Sample Inventory**:
- [x] Physical location tracking (freezers/boxes)
- [x] Sample lineage and parent/child tracking
- [x] Barcode generation support

**Compliance & Integrity**:
- [ ] 21 CFR Part 11 readiness (audit trails)
- [ ] Electronic signatures (future)
- [ ] Data immutability logs
- [ ] **IP & Patent Tracking**:
  - [ ] Invention Disclosure Registry (Idea submission, inventors)
  - [ ] Patent Portfolio Manager (Filings, dates, status)
  - [ ] Experiment-to-IP Linking (Trace data to claims)

**General Experimental Data Support**:
- [ ] Generic Assay Result Model (for non-standard assays)
- [ ] Imaging Data Support (Microscopy/HCS metadata)
- [ ] Flow Cytometry / FACS Data Ingestion
- [ ] Biophysical Assay Support (SPR, MST, DSC)
- [ ] **Genomic & Variant Tracking**:
  - [ ] Mutation/Genotype Model (Gene, Variant, Zygosity)
  - [ ] Cell Line/Organism Genotype Mapping (e.g., "HeLa -> p53 null")
  - [ ] Variant-based Experiment Filtering

**Resource Management**:
- [ ] **Cost Tracking & Budgeting**:
  - [ ] Project Budget Allocation (CapEx/OpEx tracking)
  - [ ] Experiment Cost Calculator (Reagents + Labor + Outsourcing)
  - [ ] Purchase Request Approval Workflow
- [ ] **Experiment Scheduling**:
  - [ ] Shared Calendar / Equipment Booking
  - [ ] Gantt Chart view for Protocol timelines
  - [ ] Personal Task Scheduler (integrates with Cockpit)
- [ ] **Chemistry Procurement**:
  - [ ] Vendor Catalog Search (MolPort, Mcule, Enamine)
  - [ ] Availability & Pricing Check
  - [ ] Shopping Cart & PO Generation

**Scientific Q&A & Insight Tracker**:
- [x] Question Registry (CRUD for scientific questions)
- [x] Persistent RAG Answers (Save outputs + evidence)
- [x] Versioning (Re-run questions as data updates)
- [x] Export capabilities (PDF/CSV of Q&A reports)
- [ ] **Literature Critical Analysis**:
  - [x] Auto-generated critiques (Strengths, Weaknesses, Limitations)
  - [x] Unanswered Questions Extraction
  - [x] Contradiction Detection (Paper A vs Paper B)
- [ ] **Automated Study Critique**:
  - [x] Quality Assessment for internal/imported studies
  - [x] Design Flaw Detection (e.g., "Low N", "Missing Control")
  - [x] Data Gap Identification (Unanswered questions)

### 3. Chemistry Informatics Enhancements

**Chemical Registration**:
- [x] Corporate ID generation (e.g. AMP-001)
- [x] Salt and batch handling
- [x] Duplicate checking logic

**Cheminformatics Search**:
- [x] Substructure search (SMARTS)
- [x] Similarity search (Tanimoto)
- [x] Pharmacophore search

**SAR Analysis**:
- [x] Calculated properties (LogP, MW, TPSA)
- [x] Activity cliffs visualization
- [x] R-group decomposition

**Lead Optimization Data Models**:
- [x] ADME/DMPK Assays (Stability, Permeability, CYP inhibition)
- [x] In Vivo PK Studies (AUC, Cmax, Tmax, Bioavailability)
- [x] Safety & Toxicology (hERG, Ames, Cytotoxicity)
- [ ] **Candidate Selection Workflow**:
  - [x] Target Product Profile (TPP) Definition (Criteria & Thresholds)
  - [x] Development Candidate (DC) Nomination Checklist
  - [x] Traffic Light Scoring (Green/Yellow/Red vs TPP)

**Structural Biology & Virtual Screening**:
- [ ] Protein Structure Visualization (Mol* / PyMOL integration)
- [ ] Virtual Docking Pipeline (AutoDock Vina / Smina)
- [ ] Binding Site Analysis & Pocket Detection
- [ ] Ligand-Protein Interaction Profiler

### 4. RAG Maturity & Enhancements

**Retrieval Optimization**:
- [x] Hybrid Search (Sparse/BM25 + Dense Vectors)
- [x] Reranking (Cross-encoder re-ranking)
- [x] Query Expansion/Transformation (HyDE)

**Quality & Attribution**:
- [x] Citation/Source Attribution
- [x] Evaluation Framework (RAGAS/TruLens)
- [x] Hallucination detection

**Advanced Capabilities**:
- [x] Semantic Caching
- [x] Agentic RAG (multi-step reasoning)
- [x] **Hierarchical RAG Reasoning**:
  - [x] Data Trust Scoring (Validated Internal > External > General)
  - [x] Hybrid Prompting (Allow general knowledge as fallback)
  - [x] Source-weighted Generation
- [ ] **Multi-Model Intelligence**:
  - [x] Model Selection UI (GPT-4o, Claude 3.5, Gemini 1.5)
  - [x] Parallel Reasoning (Run multiple models -> Synthesize)
  - [x] Chain-of-Thought (CoT) Prompting options

### 6. Operational Maturity & Governance
**Data Governance**:
- [x] Data Lineage Visualization (Provenance graph)
- [ ] Retention & Archival Policies
- [ ] Internal Ontology Management
- [x] **Concurrent Editing Safety**:
  - [x] Optimistic Locking (Versioning/ETags for all models)
  - [x] Conflict Resolution UI ("Your changes conflict with...")
  - [x] Real-time Presence (Who is viewing this page?) - partial, conflict detection only

**Infrastructure**:
- [ ] Job Queue System (Celery/Redis)
- [ ] S3/Tiered Storage Integration
- [ ] **Automated Backup & Disaster Recovery**:
  - [ ] Daily Full Backups to S3
  - [ ] Point-in-Time Recovery (WAL Archiving)
  - [ ] User-initiated Data Export (Project zip download)
- [ ] **System Administration**:
  - [ ] Cache Management UI (Clear Redis/Semantic caches)
  - [ ] System Health Dashboard (CPU/Memory/Queue depth)

**Collaboration**:
- [ ] "Project" Workspace Abstraction
- [x] Contextual Commenting (on graphs/features)
- [x] **Email & Notifications**:
  - [x] Share Experiment Results via Email (PDF/Link)
  - [x] Automated Alerts (New data in subscribed project)
  - [x] Daily/Weekly Digest Summaries
- [x] **Automated Reporting & Slides**:
  - [x] Auto-generate PowerPoint (.pptx) from Experiment/Dataset
  - [x] "Smart Slides": Auto-layout plots, methods, and conclusions
  - [x] Meeting Deck Builder: Aggregating results across projects

### 7. Future Innovations

**Experimental Design Assistant**:
- [ ] Power Analysis for sample size estimation
- [x] Design recommendation engine
- [x] Confounding variable detection

**Advanced Analytics**:
- [x] Multi-dataset pattern detection
- [ ] **Bayesian Inference & Optimization**:
  - [ ] Probabilistic Dose-Response (PyMC/Stan)
  - [ ] Bayesian Experimental Design (Next-best experiment suggestion)
  - [ ] Uncertainty Quantification for Predictions
- [ ] **Machine Learning & AI**:
  - [ ] Predictive ADMET (Graph Neural Networks / Chemprop)
  - [ ] Generative Chemistry (De Novo Design / VAEs)
  - [ ] Image Analysis Pipeline (CellPose for HCS)

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
3. General Experimental Data Support
4. Audit Trails & Compliance

### Phase 3: Chemistry Informatics
1. Lead Optimization Data Models (ADME/PK)
2. Chemical Registration System
3. Advanced Structure Search
4. Structural Biology Tools (Docking/Viz)
5. SAR Analysis Tools

### Phase 4: RAG Maturity
1. Hybrid Search & Reranking
2. Hierarchical Reasoning & Trust Scoring
3. Citation & Attribution System
4. Evaluation Pipeline

### Phase 5: Future Innovations
1. Multi-dataset Pattern Detection
2. Experimental Design Assistant & Power Analysis
3. Genome Browser Integration (IGV.js)

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

**Current System Status**: ✅ Production Ready (Advanced Analytics & AI Features Live - 30+ Completed Features)
**Architecture**: Postgres-only (Notion removed December 2025)
**Next Milestone**: Operational Maturity & Scalability
**Estimated Time**: 2-3 weeks
