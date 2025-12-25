# ROADMAP (Single Source of Truth)

**Last Updated**: 2025-12-25

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

### Testing & Development
- ✅ Comprehensive Test Data Seeding Suite (all domains, size presets, deterministic)
- ✅ Phase 2 Code Quality (2025-12-19): ruff.toml config, 6,252 whitespace fixes, 213 F401 fixes, 22 E712 fixes, model re-exports restored
- ✅ Test Suite: 796 unit/integration passed, 32 skipped, 76 E2E tests (requires_server marker, in CI)
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

### Testing & Development (Future)
- ❌ OOP Refactoring Review (code structure improvements)

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

### Tier 3 - AI-Assisted
- ❌ AI Dataset Finder (NL → cross-repo query)
- ❌ Metadata Enrichment from Abstracts (LLM extraction)
- ❌ AI Relevance & Novelty Scoring

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
- ❌ Potency & Liability Multi-objective Ranking
- ❌ Assay Outcome Predictors (program-specific)
- ❌ Active Learning for Screening

### Tier 3 - Platform
- ✅ Lightweight Model Registry + Reproducibility (2025-12-23)
  - MLModelRegistry in amprenta_rag/ml/registry.py
  - Postgres metadata + joblib artifacts
- ✅ AutoML Notebook Templates (2025-12-25) - see Quick Wins section
- ❌ Drift & Calibration Monitoring
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

### Tier 2 - Strategic (Deferred)
- ❌ Prior Builder (expert/literature priors)
- ❌ Posterior Inclusion Probabilities for Biomarkers
- ❌ Multi-Objective BO (qNEHVI, Pareto front visualization)
- ❌ Dashboard toggle for Bayesian dose-response in HTS QC page
- ❌ E2E tests for Bayesian endpoints

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
- ❌ Pathway Maps with Data Overlays
- ❌ Dose-Response & Time-Series Explorer (Plotly)
- ❌ Multi-Omics Integration View (Alluvial/UpSet)

### Tier 2 - Strategic
- ❌ SAR Grid & R-Group Explorer
- ❌ Interactive Plate Maps (HTS QC)
- ❌ Program/Experiment Voila Dashboards

### Libraries
- Networks: Cytoscape.js
- 3D Molecules: py3Dmol, mol*
- Charts: Plotly, Altair

## Collaboration (Future)

### Tier 1 - High Impact
- ❌ Entity-Scoped Sharing & Roles (per-entity ACLs)
- ❌ Comments, Mentions, Inline Annotations
- ❌ Review & Approval Workflows (datasets, signatures, protocols)
- ❌ Activity Feed + Notifications
- ❌ Shared Notebooks & Project Library
- ❌ Project/Team Workspaces Hub

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

