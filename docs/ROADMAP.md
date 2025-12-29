# ROADMAP (Single Source of Truth)

**Last Updated**: 2025-12-28 (Scientific Paper Ingestion added)

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
- ✅ **RAG operational** (pgvector-backed retrieval + LLM synthesis)

---

## ✅ DONE (Key Capabilities)

### Platform & Architecture
- ✅ Postgres-only architecture (Programs/Experiments/Datasets/Features/Signatures/Compounds/HTS)
- ✅ Notion removed (no runtime dependencies)
- ✅ Streamlit dashboard: modular page architecture + sessions/auth
- ✅ FastAPI service layer operational
- ✅ Environment Management (2025-12-24)
  - environment.yml with conda-forge for RDKit/OpenMM
  - requirements-pip.txt for pip-only packages
  - README updated with conda setup instructions

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
- ✅ mwTab API optimization (caching, Postgres storage, parallel fetching)
- ✅ Repository Interface (External Data Catalog, Subscriptions, Alerts, Health Monitoring) - 2024-12-18

### Jupyter Notebook & Analytics
- ✅ JupyterHub deployment with DockerSpawner and user isolation
- ✅ API client library (amprenta-client) with 6 resource clients
- ✅ Write endpoints for annotations
- ✅ Voila dashboards: HTS Plate Viewer, Compound Triage, SAR Delta Explorer
- ✅ Publishing automation (artifact registry, scheduled reports via Papermill)
- ✅ Context parameter passing and "Open in Jupyter" integration

### Advanced Analytics Features
- ✅ Signature Match Explainability (per-feature contributions, direction concordance)
- ✅ One-Click Narrative Reports (Papermill-parameterized templates)
- ✅ HTS QC & Triage Assistant (Z' factor, traffic-light scoring)
- ✅ Data Quality Watcher (monitoring and alerts)
- ✅ Protocol Version Diff & Deviations Audit
- ✅ Cross-Omics Pathway Analysis (enrichment, network proximity)
- ✅ MOA Inference from HTS + Multi-Omics (evidence fusion, probability scores)
- ✅ Multi-Objective Bayesian Optimization (2025-12-27)
  - qNoisyExpectedHypervolumeImprovement acquisition function
  - ModelListGP for multiple objectives (potency, hERG, logS, etc.)
  - Pareto front computation with is_non_dominated
  - 2D/3D Pareto visualizations in Experiment Optimizer dashboard
  - Direction controls (maximize/minimize per objective)
  - API endpoint: /api/v1/bayesian-optimization/recommend-multi-objective
  - 7 unit tests covering validation and multi-objective scenarios

### Collaboration & User Experience
- ✅ Activity Feed + Notifications (2025-12-27)
  - activity_events table for tracking actions
  - Unified notification system (merged with alerts)
  - API: 6 endpoints (feed, events, notifications, read, read-all, count)
  - Dashboard: Enhanced alerts bell with Activity tab
  - Activity Feed page with filters and timeline
  - Integration: model training, experiments, notebook reviews
  - Phase 2: compound_added, hit_confirmed, status_changed integration
- ✅ Comments + @Mentions (2025-12-27) - Inline Annotations deferred to Phase 2
  - Contextual commenting on entities (datasets, experiments, compounds, signatures)
  - @mention parsing and notification system
  - Edit and delete functionality (author-only)
  - Threaded replies with parent-child relationships
  - @mention highlighting in dashboard with styled HTML
  - Page integration (experiments, datasets, compounds)
  - API: 4 endpoints (create, list, update, delete)
  - 20 tests (8 service + 8 API + 4 E2E) - all passing

### Visualization & UI/UX
- ✅ Chemical Sketcher (2025-12-29)
  - Ketcher structure editor (CDN-based, no backend required)
  - Draw molecules to generate SMILES
  - Compound registration from drawings
  - Structure search from drawings
  - postMessage API integration for seamless data exchange
  - Dashboard page registered in PAGE_REGISTRY under Chemistry category
  - 17 tests (7 unit + 6 E2E + 4 integration) - 100% pass rate
  - P1 bug fixed: compounds.py import path issue
- ✅ High-Dimensional Projector (2025-12-29)
  - Interactive 3D scatter plots with UMAP, t-SNE, PCA algorithms
  - Real dataset integration with feature loading
  - Color-by support for categorical and continuous variables
  - 3 API endpoints (POST /compute, GET /datasets, POST /export)
  - Dashboard page with Plotly 3D visualization
  - Algorithm selector with parameter controls
  - 25 tests (10 engine + 9 API + 6 E2E) - 100% pass rate
  - P3 fixes applied: Real dataset loading, color-by visualization, dataset listing
  - Registered in PAGE_REGISTRY under Visualization category

### Knowledge Management & Literature
- ✅ Publication & Supplementary Data Extraction (2025-12-29)
  - PDF experiment extraction with LLM (PyMuPDF + GPT-4)
  - Supplementary file parsing with intelligent schema detection (Excel/CSV)
  - 4 API endpoints (POST /upload, /extract, /parse-supplementary, /link-dataset)
  - Dashboard integration: Publication upload page with extraction workflow
  - Database schema: PublicationExtraction model for tracking extraction jobs
  - 43 tests (12 PDF + 17 supplementary + 8 API + 6 E2E) - 100% pass rate
  - Automated extraction of: title, authors, methods, results, figures, tables
  - Schema auto-detection: column types, units, experiment metadata
- ✅ Semantic Scholar / OpenAlex Integration (2025-12-29)
  - Citation graph analysis and paper enrichment
  - Two repository implementations (SemanticScholarRepository, OpenAlexRepository)
  - Citation tracking database schema (PaperCitation model with direction field)
  - 3 new API endpoints (GET /citations, GET /references, POST /enrich)
  - 23 unit tests (12 repository + 6 API + 5 integration) - 100% pass rate
  - Bug fixed during review (citation direction field inverted)
  - Papers can be enriched with citation data, author networks, venue information
- ✅ Scientific Paper Ingestion (2025-12-28)
  - CRUD API for papers, authors, affiliations (7 endpoints)
  - Ingestion pipeline with reference parsing (DOI, PubMed ID, arXiv)
  - Author and affiliation automatic linking with deduplication
  - Full-text search and export functionality
  - 22 tests (7 API + 10 functional + 5 E2E) - all passing
  - Commits: b2dde2e, 8e5f81b, ede9605

### Collaboration & User Experience
- ✅ Collaboration Features MVP (2025-12-28)
  - Entity-scoped sharing (view/edit/admin permissions)
  - Review workflows (draft→submitted→in_review→approved/rejected state machine)
  - Workspaces Hub dashboard for teams and projects
  - API: 11 endpoints (5 sharing + 6 reviews)
  - 38 tests (6 service sharing + 8 service reviews + 10 API sharing + 14 API reviews)

### Testing & Development
- ✅ UI Development for API Routers (2025-12-28)
  - 4 new UI pages (Screening, Predictors, Scoring, Phenotypes)
  - 40 tests (24 E2E + 16 API)
  - 100% pass rate, zero skipped tests
  - All 4 routers that lacked UI now have full dashboard pages
  - Pages registered in PAGE_REGISTRY under Analysis Pages
  - E2E coverage: ~56% → ~58% (additional 24 tests)
  - API coverage: ~90% → ~92% (additional 16 tests)
- ✅ Test Coverage Remediation (2025-12-28)
  - 25 new test files (15 E2E + 10 API)
  - 153 tests with 100% pass rate
  - E2E coverage: 36% → ~56% (83 new tests covering Activity Feed, Digest Manager, Model Monitoring, MOA Inference, Biomarker Discovery, AI Extraction, Sync Monitor, Programs, Datasets, Experiments, Review Queue, Notebook Generator, Notebook Copilot, Nextflow Orchestrator, Pipeline Runner)
  - API coverage: 71% → ~90% (68 new tests covering digests, extraction, sync, automl, batch, pathways, protocols, quality, reports, scoring)
  - Zero skipped tests (No Bandaids policy enforced)
  - Technical debt cleaned: 29 networkidle → domcontentloaded replacements
  - Bug fixed: pathways.py endpoint List[UUID] query parameter
- ✅ Comprehensive Test Data Seeding Suite (all domains, size presets, deterministic)
- ✅ Phase 2 Code Quality (2025-12-19): ruff.toml config, 6,252 whitespace fixes, 213 F401 fixes, 22 E712 fixes, model re-exports restored
- ✅ Test Suite Baseline: 796 unit/integration passed, 32 skipped, 76 E2E tests (requires_server marker, in CI)
- ✅ Playwright E2E Separation (2025-12-19): chromium installed, server-dependent tests marked and deselected from default runs
- ✅ Mypy Type Checking Setup (2025-12-19): gradual typing with strict API layer, 54% error reduction (597→273), type stubs installed
- ✅ Code Quality Optimization (2025-12-21): 70% coverage achieved, type ignores reduced 69%, pre-commit hooks, Bandit security scanning, E2E tests in CI
- ✅ Dependency Audit (2025-12-24)
  - Security scan: 4 CVEs identified (no patches available yet)
  - 33 outdated packages documented
  - 115 unused imports removed
- ✅ Dead Code Detection (2025-12-24)
  - Vulture analysis: 99 high-confidence candidates triaged
  - 23 true dead code items removed
  - Maintenance scripts + reports created
- ✅ Comprehensive Test Data Seeding Suite (2025-12-24)
  - 13 domain seeders (added: genomics, single-cell, structures)
  - seed_utils.py: shared CLI parser, progress bars, schema detection
  - Enhanced seed_all.py: validation, summary, graceful skip
  - Idempotent design: safe to re-run without --reset
  - CI workflow: .github/workflows/seed_test_data.yml
- ✅ Technical Debt Cleanup (2025-12-27)
  - Removed pytest.skip from API/E2E tests (12 skips → assertions/fixtures)
  - Replaced 7 hardcoded auth IDs with get_current_user dependency
  - Fixed cache timing test (functional, not timing)
  - Comprehensive schema audit (48 tables added)
  - Migration reset (clean baseline e9386d67cc54 + b5d4975de50f)
  - Auto-migrate test fixture in conftest.py
  - Fixed Activity API session management
  - Fixed Datasets API schema validation
  - Fixed 4 E2E test selector failures (brittle input[aria-label] → robust get_by_text with .or_() fallbacks)
  - E2E tests: 66 total, 64 passing, 2 legitimate skips (optional ADMET UI, HPO data)
  - 250+ tests passing (API: 149, Integration: 10, Utils: 91)

### Testing & Development (Completed)
- ✅ Pytest Warnings Cleanup (2025-12-28) - 0 warnings achieved
  - Pydantic V2 migration: ConfigDict, Field validators (24+ fixes across 16 files)
  - Fixed matplotlib deprecation (cm.get_cmap)
  - Fixed duplicate router registration in main.py
  - Added filterwarnings for botorch/alembic
  - All external warnings suppressed via pytest.ini
  - Zero-warning policy enforced
- ✅ Pinecone Removal (2025-12-28)
  - Migrated to pgvector-only backend
  - Removed pinecone packages from requirements.txt
  - Created amprenta_rag/utils/metadata.py (extracted sanitize_metadata)
  - Stubbed pinecone functions with deprecation warnings
  - 1498 tests collecting successfully

### Testing & Development (Future)
- ❌ Integration Tests with Real Database
  - Current API tests use mocked DB (risk of mock drift)
  - Add integration tests that verify endpoints with real Postgres
  - CI pipeline with test database for subset of critical paths
  - Reduces "mocked tests pass, production fails" risk
- ❌ Split requirements.txt - move pdbfixer/openmm to requirements-structural.txt (pdbfixer requires conda-forge, not pip-installable)
- ❌ OOP Refactoring Review (code structure improvements)

### Future Backlog (Migrated from NEXT_STEPS.md 2025-12-28)

**Integration & Data Sources:**
- ❌ Imaging Data Support (Microscopy/HCS metadata)
- ❌ Flow Cytometry / FACS Data Ingestion
- ❌ Biophysical Assay Support (SPR, MST, DSC)

**Visualization & UI/UX:**
- ❌ Navigation & UI Organization
  - Functional grouping (Discovery, ELN, Analysis, Admin)
  - Collapsible sidebar sections
- ❌ Scientist's Cockpit Dashboard (widgets for recent data, alerts, tasks)

**Compliance & IP:**
- ❌ 21 CFR Part 11 Readiness
  - Electronic signatures
  - Data immutability logs
- ❌ IP & Patent Tracking
  - Invention Disclosure Registry
  - Patent Portfolio Manager
  - Experiment-to-IP Linking

**Infrastructure & Operations:**
- ❌ Job Queue System (Celery/Redis)
- ❌ Automated Backup & Disaster Recovery
  - Daily full backups to S3
  - Point-in-Time Recovery (WAL archiving)
  - User-initiated data export (project zip)
- ❌ Enhanced System Administration Tools
  - Cache Management UI
  - Extended health monitoring (CPU/Memory/Queue)

**Future ML/AI Innovations:**
- ❌ Generative Chemistry (De Novo Design / VAEs)
- ❌ Image Analysis Pipeline (CellPose for HCS)
- ❌ Power Analysis for Sample Size Estimation

---

## ⏳ NEXT UP (Prioritized)

### 1) Jupyter Notebook Integration

#### Phase 1: API Client Library ✓
- amprenta-client package with 6 resource clients

#### Phase 2: Write Endpoints ✓
- Annotation endpoints for datasets, experiments, signatures, compounds

#### Phase 3: JupyterHub Deployment ✓
- DockerSpawner with user isolation
- RDKit + data science stack
- Local development working

#### Phase 4: SSO Integration ✓
- ✓ JWT-based authentication from Streamlit
- ✓ TokenAuthenticator with /hub/token-login endpoint

#### Phase 5: Templates + Launch ✓
- "Open in Jupyter" button in dashboard
- Starter templates: getting_started, molecule_analysis, signature_explorer
- Context parameter passing (entity IDs)

#### Phase 6: Write-Back + Context ✓
- Pre-load notebooks with current entity context
- Save cell outputs as annotations
- Notebook attachments to entities
- "Save to RAG" magic command

#### Phase 7: Advanced Templates ✓
- Pathway enrichment analysis
- Dose-response curve fitting
- Compound similarity clustering
- Signature heatmaps with RDKit

#### Phase 8: Publishing + Automation ✅ DONE

Core Features (Completed):
- ✅ Executable Evidence Reports (Papermill → PDF/MD)
- ✅ Parameterized Voila Dashboards
- ✅ Scheduled Notebook Jobs (event + cron triggered)
- ✅ Results Write-Back + Provenance Artifacts
- ✅ Shared Template Library with Parameters

AI/LLM Integrations (Next):
- ✅ Notebook Co-Pilot (cell synthesizer/fixer) (2025-12-22)
- ✅ Query→Notebook Generator (2025-12-22)
- ✅ Auto-Explain + Summarize Notebook (2025-12-22)

Drug Discovery-Specific (Completed):
- ✅ Weekly Compound Triage Digest
- ✅ HTS QC Autopublisher
- ✅ SAR Delta Explorer (Voila)

##### Phase 8B: Voila foundation ✅ DONE

### Voila Dashboards
1. ✅ HTS Plate Viewer - plate heatmaps, Z' factor, QC at a glance
2. ✅ Compound Triage Dashboard - traffic-light scoring, Pareto trade-offs
3. ✅ SAR Delta Explorer - matched molecular pairs, R-group grids, py3Dmol
4. ✅ Signature Validation Console - per-feature contributions, approve/reject workflow
5. ✅ Pathway Impact Explorer - Cytoscape graph overlays, cross-omics toggles

### Papermill + Voila Integration ✅ DONE
- ✅ Parameterized notebooks with JSON schema
- ✅ Event triggers (new dataset/signature) + cron scheduling
- ✅ Artifact registry (HTML/JSON/PNGs) with provenance
- ✅ Cache by (notebook_path, param_hash, code_hash)

### Quick Wins (Completed 2025-12-25)
- ✅ "Open as Dashboard" button (Voila /voila/render endpoint)
- ✅ "Pin to Program" - lock dashboard to Program with config
- ✅ Weekly Executive Digests (Papermill + email/Slack notifications)
- ✅ Approval gates with signed review cards (HMAC-SHA256 signatures)
- ✅ AutoML Notebook Templates (classification, regression, clustering)

### Advanced
- ✅ Jupyter Advanced Features (2025-12-24)
  - Dashboard Discovery: Notebook Gallery with 13 templates, search/filter, JupyterHub deep links
  - SAR What-If Designer: SMILES input, real-time ADMET prediction, 2D structure rendering
  - Scaffold Hopping: benzene→pyridine, cyclohexane→piperidine transformations
  - Risk assessment with color-coded GREEN/YELLOW/RED heuristics
- ✅ AutoML Notebook Templates (2025-12-25)
  - XGBoost classification/regression
  - K-means clustering with elbow plot
  - Papermill parameters + MLModelRegistry integration

### Cross-Platform Context Pattern
URL: `?ctx=<base64url(JSON)>&sig=<HMAC>`

JSON: `{entityType, entityId, campaignId?, plateId?, version, ts}`
- DashboardRun provenance table (params_hash, code_hash, artifacts)
- Round-trip: publish writes run_id + artifacts back to source entity

### SAR What-If Designer ✅ (covered by Jupyter Advanced Features 2025-12-24)
- ✅ Sketch R-group substitutions live
- ✅ Predicted potency + ADMET traffic lights with uncertainty
- ✅ 3D conformer overlay (py3Dmol)
- ✅ One-click "save candidate" to Compounds

### Dashboard Discovery ✅ (covered by Jupyter Advanced Features 2025-12-24)
- ✅ Streamlit entity pages: "Open as Dashboard" button
- ✅ JupyterHub launcher: pinned icons
- ✅ Global Catalog page: searchable with tags
- ✅ Program/Team Library: curated per workspace

## External Repository Integration (Future)

### Tier 1 - High Impact
- ✅ External Sync Orchestrator MVP (2025-12-24)
  - Schema: SyncJob, SyncRecord, SyncConflict, ChEMBLActivity, PubChemBioassay
  - Adapters: ChEMBL (bioactivity), PubChem (bioassay)
  - Incremental sync with checksum-based change detection
  - API: 5 endpoints (run, status, list, conflicts, resolve)
  - Dashboard: Sync Monitor (jobs, conflicts, stats)
  - APScheduler integration for scheduled syncs
  - Phase 2: UniProt, KEGG, GEO incremental, auto-conflict resolution
- ❌ Provenance Ledger + Checksums (audit trail, diffs)
- ❌ GEO Incremental Harvester + Metadata Normalization
- ✅ ChEMBL/PubChem Bioactivity Sync (2025-12-24) - covered by External Sync Orchestrator MVP

### Tier 2 - Strategic
- ❌ UniProt/KEGG Mapping Refresh
- ✅ Repository Subscriptions & Alerts (saved queries) - 2024-12-18
- ✅ External Data Catalog Page (health, freshness dashboard) - 2024-12-18

### Tier 3 - AI-Assisted (Completed 2025-12-28)
- ✅ AI Dataset Finder (2025-12-28) - NL → cross-repo query
- ✅ Metadata Enrichment from Abstracts (2025-12-28) - LLM extraction
- ✅ AI Relevance & Novelty Scoring (2025-12-28)

## Machine Learning (Future)

### Tier 1 - High Impact
- ✅ ADMET SHAP Explainability (2025-12-25)
  - Feature name mapping: amprenta_rag/ml/admet/features.py (2054 names)
  - Ensemble SHAP explainer: amprenta_rag/ml/admet/explainer.py
  - API: POST /api/admet/explain
  - Dashboard: ADMET Predictor tabs (Explain + Global Importance)
  - Docs: docs/ADMET_SHAP_GUIDE.md
- ✅ Structural Alerts MVP (2025-12-25)
  - Filters: PAINS (RDKit FilterCatalog), Brenk (RDKit FilterCatalog), Lilly (SMARTS subset)
  - API: POST /api/alerts/check, POST /api/alerts/batch, GET /api/alerts/filters
  - Dashboard: Structural Alerts page + compound detail integration
- ✅ Target-Activity QSAR (2025-12-25)
  - ChEMBL integration for training data
  - TargetDatasetLoader with unit normalization + median IC50 aggregation
  - BootstrapEnsemble per target + isotonic calibration + applicability checks
  - API: /api/qsar/targets, /api/qsar/predict
  - Dashboard: Target QSAR page
  - Training CLI: scripts/train_qsar_models.py
- ✅ Signature→MOA Inference (2025-12-17) - see MOA Inference from HTS + Multi-Omics (Advanced Analytics Features)
- ✅ Biomarker Discovery (2025-12-25)
  - Statistical tests (t-test, Mann-Whitney, ANOVA) with FDR correction
  - Stability selection (bootstrap + LassoCV)
  - Cross-validated importance (RandomForest)
  - Consensus ranking across methods
  - API: POST /api/biomarker/discover, GET /api/biomarker/methods
  - Dashboard: Biomarker Discovery page (setup, methods, results tabs)

### Tier 2 - Strategic
- ✅ Potency & Liability Multi-objective Ranking (2025-12-27)
  - 5 objectives: potency, hERG, alerts, logS, logP
  - Multi-type potency normalization (IC50/EC50/Ki/Kd)
  - Weighted alert scoring (PAINS > Brenk/Lilly)
  - 4 presets: balanced, potency-first, safety-first, CNS-optimized
  - 2D Pareto visualization (liability vs potency)
  - API: 3 endpoints (score, pareto, presets)
- ✅ Assay Outcome Predictors (2025-12-28) - program-specific
- ✅ Active Learning for Screening (2025-12-28)

### Tier 3 - Platform
- ✅ Lightweight Model Registry + Reproducibility (2025-12-23)
  - MLModelRegistry in amprenta_rag/ml/registry.py
  - Postgres metadata + joblib artifacts
- ✅ AutoML Notebook Templates (2025-12-25) - see Quick Wins section
- ✅ Drift & Calibration Monitoring (2025-12-27)
  - PSI drift detection for molecular descriptors
  - FP aggregate drift (tanimoto, hamming)
  - ECE calibration for classification models
  - API: 5 endpoints (drift, calibration, log, feedback, health)
  - Dashboard: Model Health Overview, Drift Analysis, Calibration tabs
- ✅ Universal Unstructured AI Extraction (AI ETL) MVP (2025-12-24)
  - Parsers: DOCX, PPTX, Excel, CSV (PDF existing)
  - Structured LLM extraction with Pydantic schemas
  - Batch upload API + background processing
  - Dashboard: AI Extraction page (3 tabs)
  - Phase 2: OCR, web scraping, entity normalization APIs

## Bayesian Framework

### Tier 1 - High Impact (Completed 2025-12-22)
- ✅ Batch Bayesian Optimization for Screening (BoTorch qEI, single-objective)
- ✅ Hierarchical Bayesian Dose-Response (PyMC, credible intervals, ArviZ diagnostics)
- ✅ Bayesian ADMET with Uncertainty Quantification (2025-12-24)
  - Bootstrap ensemble (5 XGBoost models)
  - Isotonic/Platt calibration
  - Applicability domain checker (Tanimoto)
  - API endpoint POST /api/admet/predict
  - ADMET Predictor dashboard (3 tabs)
  - Phase 2: Train models with ChEMBL data (deferred)
- ✅ Signature→MOA Bayesian Evidence Fusion (Beta regression)

### Tier 2 - Strategic (Completed 2025-12-28)
- ✅ Prior Builder (2025-12-28) - expert/literature priors
- ✅ Posterior Inclusion Probabilities for Biomarkers (2025-12-28)
- ✅ Multi-Objective BO (qNEHVI, Pareto front visualization) (2025-12-27)
  - ModelListGP with qNoisyExpectedHypervolumeImprovement
  - Objective direction handling (maximize/minimize)
  - Auto reference point computation
  - 2D/3D Pareto visualization in dashboard
  - API endpoint: POST /api/v1/bayesian-optimization/recommend-multi-objective
- ✅ Dashboard toggle for Bayesian dose-response in HTS QC page (2025-12-28)
- ✅ E2E tests for Bayesian endpoints (2025-12-28)

### Libraries
- BoTorch + Ax + GPyTorch (optimization) ✅ Installed
- PyMC + ArviZ (probabilistic modeling) ✅ Installed

## Advanced Visualization (Future)

### Tier 1 - High Impact
- ✅ 3D Molecule Viewer + Conformer Overlay (2025-12-25)
  - Conformer generation with RDKit ETKDG + MMFF optimization
  - py3Dmol 3D visualization with style options (stick, sphere, cartoon)
  - Multi-conformer overlay and molecule comparison
  - API: POST /api/viz3d/conformers, POST /api/viz3d/overlay, GET /api/viz3d/protein/{id}
  - Dashboard: Molecule Viewer page + Compound/Protein 3D View expanders
  - Fallback to 2D image when py3Dmol unavailable
- ✅ Compound-Target Network (2025-12-25)
  - Cytoscape.js visualization with drug discovery styling
  - Node shapes by entity type (ellipse=compound, diamond=target)
  - Edge width by pIC50, color by activity type
  - IC50 range filtering and activity type filters
  - API: POST /api/network/compound-target, expand endpoints
  - Dashboard: Compound-Target Network page (3 tabs)
  - Integration: "View in Network" from compound detail
- ✅ Pathway Maps with Data Overlays (2025-12-25)
  - KEGG KGML pathway structure parsing with coordinate normalization
  - Expression data overlay with diverging colormaps
  - File-based caching with 7-day TTL
  - Node shapes: gene (rectangle), compound (ellipse), sub-pathway (hexagon)
  - Edge styles: activation (green), inhibition (red dashed)
  - API: 5 endpoints (structure, overlay, search, enriched, expression)
  - Dashboard: Pathway Map Viewer page with 4 tabs
- ✅ Dose-Response Explorer (2025-12-25)
  - 3PL/4PL/Bayesian curve fitting with scipy and PyMC
  - Bootstrap confidence intervals for parameter estimation
  - Time series trend analysis, smoothing, changepoint detection
  - Plotly visualizations with log-scale and CI bands
  - API: 4 endpoints for fitting, comparison, analysis
  - Dashboard: Dose-Response Explorer page with 4 tabs
- ✅ Multi-Omics Integration View (2025-12-27)
  - Alluvial (Sankey) diagram for feature flow visualization
  - UpSet plot for set intersections across datasets
  - Integrated into existing MOFA page as Tabs 5-6
  - API: 3 endpoints (datasets, alluvial, upset)

### Tier 2 - Strategic
- ✅ SAR Grid & R-Group Explorer (2025-12-27)
  - SAR matrix builder with pivot table
  - Plotly heatmap with pIC50 coloring
  - Murcko scaffold extraction
  - Integrated into existing SAR Analysis tab
  - API: 2 endpoints (scaffolds, grid)
- ✅ Interactive Plate Maps (HTS QC) (2025-12-27)
  - Plotly heatmap with click/box selection
  - Color modes: Normalized, Z-Score, Hit Status, Raw
  - Well detail panel with compound info
  - CSV export for selected wells
  - Auto-detect 96/384/1536-well formats
- ❌ Program/Experiment Voila Dashboards

### Libraries
- Networks: Cytoscape.js
- 3D Molecules: py3Dmol, mol*
- Charts: Plotly, Altair

## Collaboration (Future)

### Tier 1 - High Impact (Completed as MVP 2025-12-28)
- ✅ Entity-Scoped Sharing & Roles (per-entity ACLs)
- ✅ Review & Approval Workflows (datasets, signatures, protocols)
- ✅ Shared Notebooks & Project Library
- ✅ Project/Team Workspaces Hub

### Tier 2 - Strategic
- ❌ Collaborative Notebook Editing (RTC)
- ❌ Notebook Review Threads + Diffs
- ❌ Voila Share Links with Expiring Tokens
- ❌ Scheduled Review Cycles & SLAs

### Approved Features (2025-12-13) - ✅ Completed
- ✅ HTS QC & Triage Assistant (3-5 days)
- ✅ Signature Match Explainability (2-4 days)
- ✅ MOA Inference from HTS + Multi-Omics (5-8 days)
- ✅ Cross-Omics Pathway Analysis (5-8 days)
- ✅ One-Click Narrative Reports (2-3 days)
- ✅ Data Quality Watcher (2-3 days)
- ✅ Protocol Version Diff & Deviations Audit (2-4 days)

#### Feature Specs (Approved)

##### MOA Inference (5-8 days MVP)
Data Inputs:
- HTS: HTSResult, BiochemicalResult (IC50/EC50/Ki/Kd)
- Omics: Dataset features, Signature library
- Mappings: Feature↔Target↔Pathway, gene↔protein
- Literature: RAG chunks with target/signature co-mentions

Approach:
- Candidate generation: targets, pathways, signatures
- Evidence features: omics concordance, bioactivity support, pathway enrichment, network proximity, literature score, chemotype similarity
- Fusion: Calibrated logistic regression (MVP) → P(MOA)
- v1: Bayesian logistic (PyMC) with posteriors and CIs

Output JSON:
- candidate_id, type, probability, ci_low/high, rank
- contributions: {feature_name, value, weight}
- evidence_links: plots, tables, citations

Reuses: `signature_scoring.py`, `pathway/enrichment.py`, `query/rag/match_processing.py`

##### Signature Match Explainability (2-4 days)
- Per-feature contribution: sign × weight × presence × direction_match
- Visuals: lollipop chart (top ± contributors), direction concordance bar
- Extend `signature_scoring.py` to emit contributions
- UI expander on signature/dataset pages

##### One-Click Narrative Reports (2-3 days)
Sections:
- Executive summary (go/no-go bullets)
- Data provenance (tables, checksums)
- Key findings (signatures, pathways, MOA hypotheses)
- Visuals (auto-insert plots)
- Methods (dataset, stats)
- Next steps

Implementation:
- Papermill-parameterized template
- nbconvert → HTML/PDF
- Store artifact + run metadata
- "Publish to Program" button

### 2) AWS Deployment / Infrastructure Hardening
- ✅ IaC (Terraform) — Lightsail + RDS configuration (Phase 1)
- ✅ ECS/RDS/ALB deployment architecture (Phase 2, 2025-12-21)
  - VPC with public/private subnets
  - ECS Fargate cluster (API + Dashboard services)
  - Application Load Balancer with path-based routing
  - RDS PostgreSQL in private subnets
  - AWS Secrets Manager integration
  - CloudWatch Logs and monitoring
- ✅ CI/CD pipeline (GitHub Actions → deploy)
  - **Note**: CD workflow requires GitHub secrets (LIGHTSAIL_HOST, LIGHTSAIL_SSH_KEY) - pending infrastructure setup

### 3) Advanced Visualization Suite
- ✅ Cytoscape.js networks (relationships/lineage)
- ✅ IGV.js genome/variant browser
- ✅ Ag-Grid for large-table browsing/editing

### 4) Testing Closure (Implemented-but-Untested)
- ✅ Playwright E2E: **Concurrent Editing / optimistic locking** flow (2025-12-22)

### 5) Multi-Company Support (Multi-Tenancy)
- ✅ Multi-Company Support (Multi-Tenancy) (2025-12-23)
  - Company model with subdomain routing, branding, limits
  - Row-Level Security (RLS) policies on 14 tenant tables
  - Company context dependency for transparent data isolation
  - 7 API endpoints for company + user management
  - Company Settings dashboard (2 tabs)

---

## ✅ DONE (Completed from Backlog)

### Infrastructure (Completed 2025-12-22)
- ✅ Bioinformatics Pipeline Runner Dashboard - Salmon/Kallisto quantification UI with job tracking
- ✅ Bioinformatics pipeline orchestrator (Nextflow/Snakemake)
- ✅ Data Version Control (DVC) - MVP local versioning; S3 remote Phase 2
- ✅ Pinecone → pgvector migration - Self-host vector search in Postgres
- ✅ Bayesian inference & optimization workflows
  - Hierarchical Bayesian Dose-Response (PyMC, ArviZ)
  - Bayesian Optimization for Screening (BoTorch qEI)
  - MOA Bayesian Evidence Fusion (PyMC Beta regression)
- ✅ ML/AI model registry + predictive ADMET (2025-12-23)
  - Model Registry: Postgres metadata + joblib artifacts
  - ADMET Predictor: Morgan FP + RDKit descriptors

### Bioinformatics & Integration (Completed 2025-12-23)

#### High Priority

#### Evidence Graph + GraphRAG (Track A, 2-3 weeks)
- ✅ A1: Graph Schema + ETL (2025-12-23)
  - GraphEdge table with composite indexes
  - 4 edge types: Compound↔Target, Target↔Pathway, Signature↔Pathway, Compound↔Signature
  - EdgeBuilder service with upsert + neighbor queries
  - API: GET /api/graph/neighbors
- ✅ A2: Graph API + Visualization (2025-12-23)
  - k-hop traversal (BFS, max 500 nodes)
  - Shortest path (bidirectional BFS, 5s timeout)
  - Cytoscape.js JSON with entity-type colors
  - Graph Explorer page (streamlit-agraph)
  - Deep-link support + "View in Graph" integration
- ✅ A3: GraphRAG Reranker (2025-12-23)
  - Path-aware boosting for RAG queries
  - Hybrid entity extraction (LLM + regex fallback)
  - Alpha blending: score * (1 + 0.3 * graph_boost)
  - LRU cache (20k entries) for performance
- ✅ A4: Graph Analytics (2025-12-23)
  - Degree centrality + Louvain community detection
  - POST /api/graph/analytics endpoint
  - Graph Explorer: "Size by Degree" + "Color by Community" toggles
  - Modularity legend display

- ✅ Connectivity Mapping (LINCS/CMap) (2025-12-23)
  - Data: LINCS L1000 Level 5 (978 landmark genes)
  - Schema: LINCSGene, LINCSSignature, ConnectivityScore tables
  - Algorithm: Weighted cosine similarity (CMap standard)
  - API: 6 endpoints (ingest, list, compute, reversals, mimics, mechanisms)
  - Dashboard: Connectivity Map (3 tabs)
  - Tests: 4 passed (scorer, parser, E2E)
- ✅ Single-Cell Omics Integration (2025-12-23)
  - Scanpy pipeline: QC, normalize, HVG, PCA, UMAP, Leiden clustering
  - 8 API endpoints: ingest, list, clusters, recluster, markers, umap, expression
  - Dashboard: 4-tab viewer (Dataset, UMAP overlay, Clusters, Markers)
  - Gene-to-Feature mapping for marker discovery
- ✅ Batch Effect Correction (2025-12-23) - ComBat harmonization pipeline
- ✅ Sphingolipid Pathway Imbalance Scoring (2025-12-23) - Enzyme coverage + metabolite ratio scoring

#### Medium Priority
- ✅ Variant Interpretation (VEP/ClinVar) (2025-12-23)
  - VEP TSV parser with SIFT/PolyPhen extraction
  - ClinVar integration (download + match)
  - Gene burden scoring (pathogenic×10 + VUS)
  - 6 API endpoints + 4-tab dashboard with pathway enrichment
- ✅ Phenotype→Gene Mapping (2025-12-23) - HPO gene associations, query expansion API
- ✅ CRISPR Screen Analysis (MAGeCK) (2025-12-23)
  - MAGeCK test (RRA) with Docker-first runner
  - Gene-level results with FDR-based hit calling
  - 6 API endpoints + 3-tab dashboard (screens, results, volcano)
- ✅ Multi-Omics Latent Factors (MOFA2) (2025-12-23)
  - MOFA2 Bayesian factor analysis with sample alignment
  - Variance explained per omics type
  - 7 API endpoints + 4-tab dashboard (setup, variance, scatter, loadings)
- ✅ Lipidomics Spectral Matching (2025-12-23)
  - MGF parser for LipidBlast/GNPS spectral libraries
  - Greedy cosine similarity with precursor filtering (10 ppm)
  - Confidence scoring (cosine > 0.7, m/z error < 10 ppm)
  - 5 API endpoints + 2-tab dashboard with Accept/Reject review

#### Structural Biology & Virtual Screening (Track B, 2 weeks)
- ✅ B1: Protein Structure Store (2025-12-23)
  - Fetch: RCSB PDB + AlphaFold DB (v4 models)
  - Preparation: pdbfixer (add H, missing residues, chain selection)
  - Metadata: sequence (SEQRES), resolution, method extraction
  - API: 6 endpoints (list, get, fetch, upload, prepare, download)
  - Dashboard: Fetch UI, prep button, download links
  - Tests: 7 passed (fetchers, storage, prep, parser, E2E)
- ✅ B2: Binding Site Detection (2025-12-23)
  - Tool: fpocket (geometry-based, Docker-first with local fallback)
  - Output: Top 5 pockets sorted by druggability score
  - Storage: BindingSite table with centroid coords + pocket PDB
  - API: detect-pockets, list, download endpoints
  - Dashboard: Structure selector, detect button, table view
  - Tests: 5 passed (unit + E2E)
- ✅ B3: Structure-Based Virtual Screening (2025-12-23)
  - Tool: AutoDock Vina (Docker-first, CPU-based)
  - Ligand prep: RDKit ETKDG + Meeko/OpenBabel PDBQT
  - Batch: ThreadPoolExecutor (4 workers), atomic progress
  - API: 5 endpoints (create run, list, poses, download)
  - Dashboard: Runs page + Triage page (hit filtering)
  - Tests: 9 passed (ligand prep, vina, API, E2E)
- ✅ B4: Pose QC & Interactions (2025-12-23)
  - Tool: PLIP (Docker-first, 7 interaction types)
  - QC: Ligand efficiency, clash detection
  - Storage: PoseQuality + PoseInteraction tables
  - API: 3 endpoints (analyze, quality, interactions)
  - Tests: 3 passed (PLIP runner/parser, E2E)

---

## Notes

- The full, detailed product checklist remains in `NEXT_STEPS.md`.
- Historical roadmap documents are archived under `md_archive/roadmaps/`.

