# Unified Strategic Roadmap - Amprenta Multi-Omics Platform

**Status**: Comprehensive roadmap integrating strategic analysis + ChatGPT roadmap

**Purpose**: Single source of truth for all future enhancements, prioritized by value and aligned with platform goals.

---

## 🎯 Roadmap Organization

This roadmap integrates:
- ✅ Strategic feature analysis (10 high-value enhancements)
- ✅ ChatGPT's detailed roadmap (7 major categories)
- ✅ Implementation priorities
- ✅ Notion agent instructions where needed
- ✅ Dependencies and sequencing

---

## 🔥 TIER 1: IMMEDIATE HIGH VALUE (Next 2-4 Weeks)

### 1.1 Enhanced Dataset Feature Extraction & Caching ⚡ **PERFORMANCE CRITICAL**

**Priority**: 🔥🔥🔥 **HIGHEST**

**Strategic Value**:
- **10-100x performance improvement** for signature scoring
- Reduces Notion API load significantly
- Enables real-time batch operations

**Current Gap**:
- Every signature scoring operation queries Notion APIs on-demand
- No caching of extracted feature sets
- Limited direction/fold-change extraction

**Implementation**:
```python
# New module: amprenta_rag/ingestion/dataset_feature_cache.py
- In-memory cache with TTL
- Cache dataset feature sets by omics type
- Batch pre-loading support
- File-based feature extraction option
```

**Integration Points**:
- `amprenta_rag/ingestion/multi_omics_scoring.py`
- `amprenta_rag/ingestion/signature_matching.py`

**Notion Agent**: None needed ✅

**Dependencies**: None

**Estimated Effort**: Medium (3-5 days)

---

### 1.2 Batch Ingestion Framework ⚡ **OPERATIONAL EFFICIENCY**

**Priority**: 🔥🔥🔥 **HIGH**

**Strategic Value**:
- Faster bulk ingestion
- Better operational efficiency
- Progress visibility

**ChatGPT Roadmap Alignment**: Section 2.2

**Implementation**:
```python
# New script: scripts/batch_ingest_omics.py
- Auto-detect omics type from files
- Parallelize chunk creation + upsert
- Directory scanning
- Progress bars/logging
- Error aggregation
- Per-file summaries
```

**Features**:
- Ingest multiple files from directory
- Auto-detect omics type (lipidomics, metabolomics, proteomics, transcriptomics)
- Parallel processing where safe
- Comprehensive logging

**Notion Agent**: None needed ✅

**Dependencies**: None

**Estimated Effort**: Low-Medium (2-3 days)

---

### 1.3 Enhanced Cross-Omics RAG Reasoning ⚡ **ALREADY IMPLEMENTED, ENHANCE**

**Priority**: 🔥🔥 **MEDIUM-HIGH**

**Strategic Value**:
- Better scientific summaries
- Disease/model context awareness
- More actionable insights

**ChatGPT Roadmap Alignment**: Section 1.3

**Current Status**: Functions exist, need enhancement

**Enhancements Needed**:
```python
# Enhance amprenta_rag/query/cross_omics_reasoning.py
- Include disease context in summaries
- Model system awareness (in vitro, in vivo, patient)
- Matrix context (CSF, plasma, tissue)
- Comparative analysis (disease vs control)
- Multi-modal retrieval builder utilities
- Prompt templates for each use case
- Unit tests around retrieval/summary logic
```

**Notion Agent**: None needed (uses existing properties) ✅

**Dependencies**: Existing cross-omics functions

**Estimated Effort**: Low (2-3 days for enhancements)

---

## ⭐ TIER 2: STRATEGIC CAPABILITIES (Next 4-8 Weeks)

### 2.1 Automated Signature Discovery from Datasets 🔬 **DISCOVERY ACCELERATION**

**Priority**: 🔥🔥🔥 **HIGH**

**Strategic Value**:
- Accelerate signature library growth 10x
- Discover patterns humans might miss
- Enable data-driven signature generation

**ChatGPT Roadmap Alignment**: Section 3.3

**Current Gap**:
- Signatures must be manually defined
- No automatic pattern detection across datasets

**Implementation**:
```python
# New module: amprenta_rag/signatures/signature_discovery.py
# (Also: amprenta_rag/signatures/auto_discovery.py per ChatGPT roadmap)
- Analyze multiple datasets for recurring patterns
- Detect co-occurring features across omics
- Cluster features across datasets
- Identify recurring cross-omics patterns
- Direction consistency checking
- Generate candidate signature files
- Write to Notion as "Suggested Signatures" (optional)
```

**Features**:
- Statistical pattern detection
- Feature co-occurrence analysis
- Direction consistency scoring
- Export candidate signatures as TSV
- Optional Notion integration for suggestions

**Notion Agent**: Optional - Add "Suggested Signatures" status to Signatures DB

**Dependencies**: Statistical libraries (scipy, numpy)

**Estimated Effort**: High (5-7 days)

---

### 2.2 Evidence Report Engine 📊 **REPORTING CAPABILITY**

**Priority**: 🔥🔥 **MEDIUM-HIGH**

**Strategic Value**:
- Automated evidence summaries
- Program/Experiment/Dataset reports
- PDF export capability (future)

**ChatGPT Roadmap Alignment**: Section 1.2

**Implementation**:
```python
# New module: amprenta_rag/reporting/evidence_report.py
- Generate cross-omics summaries for:
  - Program
  - Experiment
  - Dataset
  - Signature
  - Feature (gene/protein/metabolite/lipid)
- Output to:
  - Notion pages (rich text)
  - PDF (optional future)
- Integrate with RAG Engine for context retrieval
```

**Features**:
- Reuse cross-omics reasoning functions
- Structured report templates
- Notion page creation/updates
- Export capabilities

**Notion Agent**: Optional - Create "Evidence Reports" database

**Dependencies**: Cross-omics reasoning (Tier 1.3)

**Estimated Effort**: Medium (4-5 days)

---

### 2.3 Program-Level Multi-Omics Signature Maps 🗺️ **DASHBOARD INTELLIGENCE**

**Priority**: 🔥🔥 **MEDIUM-HIGH**

**Strategic Value**:
- Visual signature-program relationships
- Omics coverage visibility
- Cross-omics convergence indicators

**ChatGPT Roadmap Alignment**: Section 1.1

**Implementation**:
```python
# New module: amprenta_rag/analysis/program_signature_maps.py
- Compute Program × Signature matrices for all modalities
- Add scoring visualizations and summaries
- Populate Notion dashboard views:
  - Strongest matching signatures
  - Omics coverage of each program
  - Cross-omics convergence indicators
- Automatically recompute when datasets are updated
```

**Features**:
- Program-signature scoring matrices
- Coverage analysis
- Convergence detection
- Notion dashboard updates

**Notion Agent**: Optional - Add summary properties to Programs DB

**Dependencies**: Signature scoring, cross-omics reasoning

**Estimated Effort**: Medium (4-5 days)

---

### 2.4 Dataset Comparison & Clustering 📊 **ANALYTICS CAPABILITY**

**Priority**: 🔥🔥 **MEDIUM**

**Strategic Value**:
- Identify similar datasets across studies
- Find complementary datasets
- Better experimental planning

**Implementation**:
```python
# New module: amprenta_rag/analysis/dataset_comparison.py
- Compare 2+ datasets across all omics types
- Compute dataset similarity scores
- Find shared/differential features
- Cluster datasets by similarity
- Generate comparison reports
- Feature set comparison algorithms
- Jaccard similarity for feature overlap
- Direction-aware comparison
```

**Notion Agent**: 
- Add "Similar Datasets" relation to Experimental Data Assets
- Add "Dataset Comparison Score" property (optional)

**Dependencies**: Feature extraction

**Estimated Effort**: Medium (3-4 days)

---

### 2.5 Cross-Omics Pathway Analysis 🧬 **BIOLOGICAL INTERPRETABILITY**

**Priority**: 🔥🔥🔥 **HIGH STRATEGIC VALUE**

**Strategic Value**:
- Biological context for multi-omics data
- Pathway-level insights (more interpretable than features)
- Better scientific reasoning

**ChatGPT Roadmap Alignment**: Section 3.2 (Enrichment Analysis)

**Implementation**:
```python
# New module: amprenta_rag/analysis/pathway_analysis.py
# New module: amprenta_rag/analysis/enrichment.py (per ChatGPT roadmap)
- Map features to known pathways (KEGG, Reactome, etc.)
- Detect pathway-level disruptions across omics
- Cross-omics pathway enrichment
- Pathway-based dataset/signature scoring
- Overrepresentation tests
- Link enriched features back to Notion (optional)
```

**Features**:
- KEGG/Reactome API integration
- Feature-to-pathway mapping
- Pathway enrichment analysis
- Pathway-aware cross-omics summaries
- Enrichment statistics

**Notion Agent**: **YES - REQUIRED** ✅
- Create "Pathways" database
- Add pathway relations to features
- Add pathway properties to datasets/signatures

**Dependencies**: Pathway database APIs (KEGG, Reactome)

**Estimated Effort**: High (6-8 days)

---

### 2.6 Experimental Design Metadata System 🔬 **STRATEGIC CAPABILITY**

**Priority**: 🔥🔥🔥 **HIGH**

**Strategic Value**:
- Enable proper case/control comparisons
- Support time course analysis
- Track intervention groups (drug A vs B vs control)
- Understand experimental context for accurate interpretation
- Enable design-aware statistical analysis

**Current Gap**:
- No structured capture of experimental design
- Case/control relationships not explicit
- Time course data not tracked
- Intervention groups not modeled

**Implementation**:

#### Phase 1: Database Schema (1-2 days)

```python
# Extend experiment/dataset models:
- design_type (enum: case_control, time_course, intervention, dose_response, multi_factorial)
- sample_groups (JSON array of groups with labels and sample IDs)
- timepoints (JSON array for time course studies)
- interventions (JSON array for treatment groups)
- control_group (reference to control group ID)
```

**Schema Extensions**:
- `experiments` table:
  - `design_type` (text)
  - `design_metadata` (JSONB)
- `datasets` table:
  - `sample_group` (text, e.g., "case", "control", "drug_a", "timepoint_0h")
  - `timepoint` (text, for time course)
  - `intervention` (text, for intervention studies)

#### Phase 2: Repository Metadata Extraction (2-3 days)

```python
# New module: amprenta_rag/ingestion/design_extraction.py
- Parse GEO sample groups from SOFT files
- Parse MW study design from metadata
- LLM-based design extraction as fallback
- Automatic case/control detection
- Timepoint extraction from sample names
```

**Features**:
- GEO GSM group parsing
- MW analysis name parsing
- Pattern matching for common designs
- GPT-4 extraction for complex designs

#### Phase 3: Design-Aware Statistical Analysis (2 days)

```python
# New module: amprenta_rag/analysis/design_aware_stats.py
- Case vs control comparisons (respects design)
- Time course analysis (repeated measures, trend tests)
- Intervention group comparisons (multi-group ANOVA)
- Pairwise comparisons with design context
```

**Design Types Supported**:

1. **Case vs Control**
   - Two-group comparisons (t-test, Mann-Whitney)
   - Explicit case/control labels
   - Example: ALS patients (n=30) vs healthy controls (n=30)

2. **Time Course (Longitudinal)**
   - Repeated measures ANOVA
   - Trend analysis over time
   - Example: Drug response at 0h, 4h, 8h, 24h

3. **Intervention Groups**
   - Multi-group comparisons
   - Drug A vs Drug B vs control
   - Example: Vehicle, Compound 1 (10µM), Compound 1 (100µM)

4. **Dose Response**
   - Dose-dependent trend tests
   - EC50/IC50 estimation
   - Example: 0, 1, 10, 100, 1000 nM doses

5. **Multi-factorial**
   - Two-way ANOVA
   - Interaction effects
   - Example: Disease (ALS vs control) × Treatment (treated vs untreated)

**Notion Agent**: **YES - OPTIONAL** ✅
- Add "Design Type" property to Experiments/Datasets
- Add "Sample Groups" property (rich text/JSON)
- Add "Timepoints" property (for time course)
- Add "Control Group" relation

**Dependencies**: 
- Repository metadata parsers (GEO, MW)
- Statistical libraries (scipy, statsmodels)
- Optional: GPT-4 for complex design extraction

**Estimated Effort**: Medium (5-6 days total)

---

### 2.7 Chemistry & HTS Integration 🧪 **NEW ARCHITECTURE**

**Priority**: 🔥🔥🔥 **HIGH STRATEGIC VALUE**

**Strategic Value**:
- Complete drug discovery pipeline integration
- HTS data management (up to 1M molecules)
- Compound-to-biology linking
- Campaign-to-program integration

**Architecture Principles**:
1. **SQLite is the system of record** for chemistry and screening data
2. **Notion is the knowledge lens** - only "promoted" compounds appear in Notion
3. **Pinecone receives summaries only** - not full HTS data matrices
4. **HTS data NOT fully mirrored** - only summaries and selected hits

**Implementation**:

#### Phase 1: SQLite Chemistry & Screening Layer

```python
# New modules:
- amprenta_rag/chemistry/compound_registry.py
- amprenta_rag/chemistry/screening_db.py
- amprenta_rag/chemistry/normalization.py (RDKit-based, optional)
```

**Schema (First Pass)**:
- `compounds` table (compound_id, smiles, inchi_key, salt_stripped_smiles, source)
- `libraries` table (library_id, name, vendor, description)
- `hts_campaigns` table (campaign_id, program_id, assay_name, target, library_id, screen_date, num_molecules, hit_cutoff, num_hits, qc_metrics_json)
- `hts_results` table (campaign_id, compound_id, plate_id, well_id, raw_signal, normalized_value, zscore, is_hit)
- `biochemical_results` table (assay_id, campaign_id, compound_id, ic50, ec50, curve_fit_quality, replicate_count, qc_flags, stage)
- `compound_program` table (compound_id, program_id, role)

**Features**:
- RDKit-based SMILES normalization (if available)
- Internal compound ID generation (AMPR-000001 pattern)
- Indexes on compound_id, campaign_id, program_id

**Notion Agent**: **YES - REQUIRED** ✅
- Create "🧪 Compound Features" database
- Create "🧪 HTS Campaigns" database
- Create "🧪 Biochemical Hits" database (optional)
- See detailed schema in `CHEMISTRY_HTS_ARCHITECTURE.md`

**Dependencies**: SQLite, RDKit (optional), Notion schema updates

**Estimated Effort**: High (8-10 days)

#### Phase 2: Screening Ingestion Pipeline

```python
# New modules:
- amprenta_rag/ingestion/screening_ingestion.py
- scripts/ingest_screening.py
```

**Features**:
- Ingest HTS campaign summary files (metadata, QC)
- Ingest hit lists and dose-response data into SQLite
- Create/refresh HTS Campaigns DB entries in Notion (summary only)
- Maintain vendor ID → internal compound_id mapping

**CLI Use Cases**:
```bash
# Register/Update Campaign
python scripts/ingest_screening.py --campaign-metadata-file <path>

# Ingest Hit List
python scripts/ingest_screening.py --hit-list-file <path> --campaign-id <id>

# Ingest Dose-Response
python scripts/ingest_screening.py --dose-response-file <path> --campaign-id <id>

# Promote Compounds
python scripts/ingest_screening.py --promote-compounds --program-id <id> --stage <stage>
```

**Notion Agent**: Uses existing Compound Features, HTS Campaigns, Biochemical Hits DBs

**Dependencies**: Phase 1 (SQLite layer)

**Estimated Effort**: Medium (4-5 days)

#### Phase 3: Compound-Biology Linking

**Features**:
- Link compounds to Programs/Experiments/Datasets
- Link compounds to Signatures (if compound affects signature features)
- Cross-reference compound features with omics features
- RAG embedding of compound summaries (campaign summaries, biological assay summaries)

**Integration Points**:
- `amprenta_rag/ingestion/compound_linking.py`
- Extend existing feature linking to support compounds
- Extend RAG embedding to include compound context

**Notion Agent**: Uses existing relation properties

**Dependencies**: Phase 1, Phase 2

**Estimated Effort**: Medium (3-4 days)

#### Phase 4: RAG Integration for Chemistry

**Features**:
- Embed HTS campaign summaries into Pinecone
- Embed biochemical hit summaries
- Query support: "Which compounds match this signature?", "What HTS campaigns tested this target?"
- Cross-omics + chemistry reasoning

**Integration Points**:
- `amprenta_rag/query/chemistry_query.py`
- Extend `cross_omics_reasoning.py` to include compound context

**Dependencies**: Phase 1, Phase 2, Phase 3

**Estimated Effort**: Medium (3-4 days)

**Total Estimated Effort**: 18-23 days (Weeks 13-16)

**See**: `context/CHEMISTRY_HTS_ARCHITECTURE.md` for complete details

---

### 2.8 Public Repository Ingestion for Multi-Omics 🌐 **DATA EXPANSION**

**Priority**: 🔥🔥 **MEDIUM-HIGH**

**Strategic Value**:
- Expand dataset coverage across all omics types
- Leverage public data repositories
- Automate discovery and ingestion
- Build comprehensive multi-omics knowledge base

**Current State**:
- ✅ Metabolomics Workbench (MW) - lipidomics and metabolomics
- ✅ `scripts/harvest_mw_studies.py` - working MW ingestion
- ⏳ GEO (transcriptomics) - not yet implemented
- ⏳ PRIDE (proteomics) - not yet implemented
- ⏳ MetaboLights (metabolomics) - not yet implemented

**Implementation**:

#### Phase 1: Repository Abstraction Layer

```python
# New modules:
- amprenta_rag/ingestion/repositories/base.py
- amprenta_rag/ingestion/repositories/geo.py
- amprenta_rag/ingestion/repositories/pride.py
- amprenta_rag/ingestion/repositories/metabolights.py
- amprenta_rag/ingestion/repositories/mw.py (refactored)
- amprenta_rag/ingestion/repositories/discovery.py
```

**Features**:
- Unified repository interface
- StudyMetadata model
- Unified discovery script
- Repository-specific harvest scripts

**Estimated Effort**: Medium (3-4 days)

#### Phase 2: GEO (Transcriptomics) Integration

**Repository**: Gene Expression Omnibus (GEO)
- Entrez E-utilities API
- GSE/GSM identifiers
- Microarray and RNA-seq data

**Implementation**:
- GEO API client
- Study discovery and metadata fetching
- Data file download (processed DGE tables)
- Integration with transcriptomics ingestion

**Estimated Effort**: Medium (4-5 days)

#### Phase 3: PRIDE (Proteomics) Integration

**Repository**: PRIDE Archive
- PRIDE REST API
- Mass spectrometry proteomics
- Protein identification and quantification

**Implementation**:
- PRIDE API client
- Study discovery and metadata fetching
- Data file download (protein tables)
- Integration with proteomics ingestion

**Estimated Effort**: Medium (4-5 days)

#### Phase 4: MetaboLights Integration

**Repository**: MetaboLights
- MetaboLights REST API
- Metabolomics data (NMR and MS)

**Implementation**:
- MetaboLights API client
- Study discovery and metadata fetching
- Data file download (metabolite tables)
- Integration with metabolomics ingestion

**Estimated Effort**: Medium (3-4 days)

#### Phase 5: Extended MW Support

**Enhancement**:
- Extend MW discovery beyond lipidomics
- Support full metabolomics studies
- Better keyword filtering

**Estimated Effort**: Low (1-2 days)

**Total Estimated Effort**: 15-20 days (Weeks 17-20)

**Notion Agent**: None needed (uses existing Experimental Data Assets DB) ✅

**Dependencies**: 
- Existing ingestion pipelines
- Feature linking system
- Signature scoring system

**See**: `context/PUBLIC_REPOSITORY_INGESTION.md` for complete details

---

## 🏗️ TIER 3: ARCHITECTURE EVOLUTION (Future)

### 3.1 Postgres + FastAPI + Frontend Architecture Evolution 🏛️ **MAJOR ARCHITECTURAL CHANGE**

**Priority**: 🔥🔥 **MEDIUM-HIGH** (Future Evolution)

**Strategic Value**:
- More robust, scalable architecture
- Better multi-user support
- Cleaner separation of concerns
- Foundation for frontend development
- Production-ready infrastructure

**Current Architecture**:
- Notion: Primary knowledge graph and ELN
- SQLite: Chemistry and screening data
- Pinecone: Semantic index for RAG
- Python scripts: Ingestion pipelines

**Target Architecture**:
- **Postgres**: System of record for all structured data
- **FastAPI**: Service/API layer
- **Frontend**: Next.js/React (future)
- **Notion**: Knowledge/documentation layer (not primary DB)
- **SQLite**: Optional cache or prototype DB
- **Pinecone**: Semantic index (unchanged)

**Implementation Phases**:

#### Phase 1: Domain Model Extraction
- Extract and formalize domain models (Pydantic/dataclasses)
- Create stable abstraction layer
- Update ingestion modules to use domain models

**Estimated Effort**: Medium (3-4 days)

#### Phase 2: Postgres Schema Design
- Design normalized Postgres schema
- Create SQLAlchemy models
- Set up Alembic migrations
- Maintain external IDs for migration support

**Estimated Effort**: High (5-7 days)

#### Phase 3: FastAPI Service Layer
- Implement core API endpoints (Programs, Experiments, Datasets, Features, Signatures, Compounds, Screening)
- Set up FastAPI application structure
- Add authentication (basic, later OIDC)

**Estimated Effort**: High (6-8 days)

#### Phase 4: Migration Utilities
- Export Notion + SQLite data
- Bootstrap Postgres with initial data
- Implement dual-write capability (transition phase)

**Estimated Effort**: Medium (6-8 days)

#### Phase 5: RAG Integration with Postgres
- Update RAG builders for Postgres + Notion hybrid
- DB for structured details, Notion for narratives
- Update metadata to use Postgres IDs

**Estimated Effort**: Medium (4-5 days)

#### Phase 6: Transition to Postgres SoT
- Pivot ingestion to Postgres as primary
- Keep Notion as documentation layer
- Final testing and validation

**Estimated Effort**: Low (1-2 days)

**Total Estimated Effort**: 25-34 days (Weeks 17-23)

**Key Principles**:
1. Postgres becomes system of record for structured data
2. SQLite remains internal detail or stepping stone
3. Notion remains human-friendly overlay (ELN, SOPs, dashboards)
4. Pinecone remains semantic index
5. Keep ID semantics stable across systems

**Notion Agent**: Schema changes may be needed for read-only mirroring

**Dependencies**: 
- Core features stable
- Chemistry & HTS integration complete
- Domain model understanding

**See**: `context/ARCHITECTURE_EVOLUTION_ROADMAP.md` for complete details

---

## 💡 TIER 4: QUALITY & OPERATIONS (Ongoing)

### 3.1 Quality Control Extraction 📋 **DATA QUALITY**

**Priority**: 🔥 **MEDIUM**

**Strategic Value**:
- QC metrics visibility
- Data quality tracking
- Better dataset curation

**ChatGPT Roadmap Alignment**: Section 2.1

**Implementation**:
```python
# Enhance existing ingestion modules
- Extract QC data from files/metadata:
  - RNA-seq: RIN, mapping rate, library complexity
  - Proteomics: identification rates, LFQ CVs
  - Lipidomics: TIC normalization, missingness rates
  - Metabolomics: RT drift, baseline intensity
- Store QC values in Notion properties
```

**Notion Agent**: **YES**
- Add "QC Metrics" property (JSON or rich text)
- Add "QC Status" property (select: Pass, Warning, Fail)

**Dependencies**: Omics-specific QC parsers

**Estimated Effort**: Medium (3-4 days per omics type)

---

### 3.2 Signature Validation & Quality Metrics ✅ **SIGNATURE QUALITY**

**Priority**: 🔥 **MEDIUM**

**Strategic Value**:
- Higher quality signatures
- Confidence scores for matches
- Better scientific rigor

**ChatGPT Roadmap Alignment**: Section 7.4 (Model-Based Validation)

**Implementation**:
```python
# New module: amprenta_rag/signatures/signature_validation.py
- Compute signature quality metrics:
  - Coverage (how many datasets match)
  - Specificity (false positive rate)
  - Reproducibility (consistency across datasets)
- Validate signatures against known ground truth
- Cross-validation scoring
- Bootstrap-based confidence scores
- Empirical p-values for signature scores
- Suggest signature improvements
```

**Notion Agent**: **YES**
- Add "Quality Score" property to Signatures DB
- Add "Validation Status" property (Validated, Candidate, Suggested)

**Dependencies**: Statistical validation libraries

**Estimated Effort**: Medium (4-5 days)

---

### 3.3 Cross-Feature Mapping 🔗 **FEATURE GRAPH**

**Priority**: 🔥 **MEDIUM**

**Strategic Value**:
- Richer feature relationships
- Pathway-level connections
- Biological context

**ChatGPT Roadmap Alignment**: Section 3.1

**Implementation**:
```python
# Enhance feature linking
- Gene ↔ Protein (HGNC/UniProt mapping)
- Gene ↔ Metabolite (pathway-level)
- Metabolite ↔ Lipid (shared synthesis pathways)
- Optional cross-feature relations in Notion
```

**Notion Agent**: Optional
- Add cross-feature relation properties to feature databases
- Use existing relations where possible

**Dependencies**: External mapping databases (HGNC, UniProt)

**Estimated Effort**: Medium (3-4 days)

---

### 3.4 Caching Feature Lookups 🚀 **PERFORMANCE**

**Priority**: 🔥🔥 **MEDIUM-HIGH**

**Strategic Value**:
- Faster feature lookups
- Reduced Notion API calls
- Better performance

**ChatGPT Roadmap Alignment**: Section 6.2

**Implementation**:
```python
# Enhance feature extraction modules
- In-memory mapping of name → page_id
- Avoid duplicate Notion searches
- Cache feature page lookups
- TTL-based expiration
```

**Integration Points**:
- `amprenta_rag/ingestion/feature_extraction.py`
- All omics ingestion modules

**Notion Agent**: None needed ✅

**Dependencies**: None

**Estimated Effort**: Low (1-2 days)

---

### 3.5 Retry Logic for Notion/Pinecone Errors 🔄 **RELIABILITY**

**Priority**: 🔥 **MEDIUM**

**Strategic Value**:
- Better error handling
- Resilience to rate limits
- Improved reliability

**ChatGPT Roadmap Alignment**: Section 6.3

**Implementation**:
```python
# Add retry utilities
- Exponential backoff for:
  - 429 rate limits
  - Intermittent network failures
- Configurable retry counts
- Clear error logging
```

**Integration Points**:
- Notion API client
- Pinecone client
- All ingestion modules

**Notion Agent**: None needed ✅

**Dependencies**: None

**Estimated Effort**: Low (2 days)

---

### 3.6 Performance Logging 📈 **OBSERVABILITY**

**Priority**: 🔥 **LOW-MEDIUM**

**Strategic Value**:
- Performance visibility
- Bottleneck identification
- Optimization guidance

**ChatGPT Roadmap Alignment**: Section 6.4

**Implementation**:
```python
# Add timestamps around:
- Parsing
- Chunking
- Embedding
- Upsert time
- Notion update delays
- Aggregate performance metrics
```

**Notion Agent**: None needed ✅

**Dependencies**: None

**Estimated Effort**: Low (1-2 days)

---

### 3.7 Sync Back-Pressure Indicators 🔄 **OPERATIONS**

**Priority**: 🔥 **LOW-MEDIUM**

**Strategic Value**:
- Ingestion status visibility
- Better operational awareness

**ChatGPT Roadmap Alignment**: Section 4.1

**Notion Agent**: **YES**
- Add "Ingestion Status" property (select: Pending, In Progress, Complete, Failed)
- Add "Last Ingested" property (date) - may already exist as "Last Embedded"

**Dependencies**: None

**Estimated Effort**: Low (1 day)

---

### 3.8 Auto-Link Experiments ↔ Programs 🤖 **AUTOMATION**

**Priority**: 🔥 **LOW-MEDIUM**

**Strategic Value**:
- Automatic relationship creation
- Reduced manual linking
- Better data consistency

**ChatGPT Roadmap Alignment**: Section 4.2

**Implementation**:
```python
# Enhance ingestion modules
- Try to infer Program/Experiment from metadata
- Auto-link when unambiguous
- Log when auto-linking occurs
- Manual override always available
```

**Notion Agent**: None needed (uses existing relations) ✅

**Dependencies**: Metadata extraction

**Estimated Effort**: Medium (2-3 days)

---

## 📊 TIER 4: DASHBOARD & VISUALIZATION (Optional)

### 4.1 Multi-Omics Coverage Maps 📊

**Priority**: 💡 **LOW**

**ChatGPT Roadmap Alignment**: Section 5.1

**Implementation**:
```python
# Generate Notion dashboards showing:
- Omics completeness per Program
- Dataset availability per modality
- Top-scoring signatures per disease
```

**Notion Agent**: Manual dashboard creation, or automated summary properties

**Estimated Effort**: Low-Medium (2-3 days)

---

### 4.2 Feature Recurrence Visualization 📈

**Priority**: 💡 **LOW**

**ChatGPT Roadmap Alignment**: Section 5.2

**Implementation**:
```python
# Create summary tables/JSON artifacts:
- Most recurrent genes/proteins/metabolites/lipids
- Cross-omics recurrence matrices
- Feature ↔ Omics heatmaps
- Generate during ingestion
```

**Notion Agent**: Optional summary pages

**Estimated Effort**: Medium (3-4 days)

---

## 🔮 TIER 5: STRETCH GOALS (Future)

### 5.1 Multi-Omics Embedding Fusion 🧠

**ChatGPT Roadmap Alignment**: Section 7.1

**Implementation**:
- Multi-modal embeddings for datasets
- Cross-omics similarity search
- Signature projection in embedding spaces

**Estimated Effort**: High (research + implementation)

---

### 5.2 Pathway-Level Models 🧬

**ChatGPT Roadmap Alignment**: Section 7.2

**Implementation**:
- Pathway signatures
- Mechanistic hypotheses
- Graph-guided RAG prompts

**Estimated Effort**: High (research + implementation)

---

### 5.3 Longitudinal Omics Integration ⏱️

**ChatGPT Roadmap Alignment**: Section 7.3

**Implementation**:
- Time-series transcriptomics
- Proteomics time courses
- Dynamic modeling

**Estimated Effort**: High (new capabilities)

---

### 5.4 Public Repository Ingestion 🌐

**ChatGPT Roadmap Alignment**: Section 2.3

**Implementation**:
```python
# Ingestion modules for:
- GEO / ArrayExpress (transcriptomics)
- PRIDE (proteomics)
- Metabolomics Workbench expansion (full API support)
- SRA / ENA (optional future)
```

**Features**:
- `fetch_public_dataset(<accession>)`
- Download + normalize into CSV/TSV
- Run dataset-specific ingestion pipeline

**Estimated Effort**: High (per repository)

---

## 📋 Implementation Roadmap Timeline

### Phase 1: Performance & Operations (Weeks 1-2)
1. ✅ Enhanced Dataset Feature Extraction & Caching
2. ✅ Batch Ingestion Framework
3. ✅ Caching Feature Lookups
4. ✅ Retry Logic for API Errors

### Phase 2: Discovery & Analytics (Weeks 3-4)
1. ✅ Automated Signature Discovery
2. ✅ Enhanced Cross-Omics Reasoning
3. ✅ Dataset Comparison & Clustering

### Phase 3: Reporting & Intelligence (Weeks 5-6)
1. ✅ Evidence Report Engine
2. ✅ Program-Level Signature Maps
3. ✅ Enhanced Cross-Omics Reasoning (context-aware)

### Phase 4: Biological Context (Weeks 7-8)
1. ✅ Cross-Omics Pathway Analysis
2. ✅ Enrichment Analysis Module
3. ✅ Cross-Feature Mapping

### Phase 5: Quality & Operations (Weeks 9-10)
1. ✅ Quality Control Extraction
2. ✅ Signature Validation & Quality Metrics
3. ✅ Sync Back-Pressure Indicators
4. ✅ Performance Logging

### Phase 6: Automation & Optimization (Weeks 11-12)
1. ✅ Auto-Link Experiments ↔ Programs
2. ✅ Code refactoring for shared utilities
3. ✅ Dashboard enhancements

### Phase 7: Chemistry & HTS Integration (Weeks 13-16)
1. ⏳ SQLite Chemistry & Screening Layer (Week 13-14)
2. ⏳ Screening Ingestion Pipeline (Week 14-15)
3. ⏳ Compound-Biology Linking (Week 15)
4. ⏳ RAG Integration for Chemistry (Week 16)

### Phase 8: Public Repository Ingestion (Weeks 17-20) 🌐 **NEW**
1. ⏳ Repository Abstraction Layer (Week 17)
2. ⏳ GEO (Transcriptomics) Integration (Week 18)
3. ⏳ PRIDE (Proteomics) Integration (Week 19)
4. ⏳ MetaboLights Integration + Extended MW (Week 20)

### Phase 9: Architecture Evolution (Weeks 21-27) ⭐ **FUTURE EVOLUTION**
1. ⏳ Domain Model Extraction (Week 21)
2. ⏳ Postgres Schema Design (Week 22)
3. ⏳ FastAPI Service Layer (Week 23-24)
4. ⏳ Migration Utilities (Week 25)
5. ⏳ RAG Integration with Postgres (Week 26)
6. ⏳ Transition to Postgres SoT (Week 27)

---

## 🎯 Priority Matrix

| Enhancement | Value | Effort | Dependencies | Priority |
|------------|-------|--------|--------------|----------|
| Feature Caching | ⭐⭐⭐⭐⭐ | Medium | None | 🔥🔥🔥 |
| Batch Ingestion | ⭐⭐⭐⭐ | Low-Medium | None | 🔥🔥🔥 |
| Chemistry & HTS Integration | ⭐⭐⭐⭐⭐ | High | SQLite, RDKit, Notion | 🔥🔥🔥 |
| Public Repository Ingestion | ⭐⭐⭐⭐ | Medium | GEO, PRIDE APIs | 🔥🔥 |
| Architecture Evolution | ⭐⭐⭐⭐⭐ | Very High | Postgres, FastAPI | 🔥🔥 |
| Signature Discovery | ⭐⭐⭐⭐⭐ | High | Stats libs | 🔥🔥 |
| Pathway Analysis | ⭐⭐⭐⭐⭐ | High | APIs + Notion | 🔥🔥 |
| Evidence Reports | ⭐⭐⭐⭐ | Medium | Cross-omics | 🔥🔥 |
| Program Maps | ⭐⭐⭐⭐ | Medium | Scoring | 🔥🔥 |
| QC Extraction | ⭐⭐⭐ | Medium | Per omics | 🔥 |
| Feature Lookup Cache | ⭐⭐⭐ | Low | None | 🔥 |

---

## 📋 Notion Agent Instructions Summary

### **Required Schema Changes**:

1. **Chemistry & HTS Databases** (for Chemistry Integration) ⭐ **NEW**
   - **🧪 Compound Features Database**
     - Name (title) - internal ID (AMPR-000123)
     - SMILES (text)
     - InChIKey (text)
     - Series (select)
     - Stage (select: HTS Hit, Biochemical Hit, Cell-Based, Lead, Candidate, Tool)
     - Programs (relation → Pipeline/Programs)
     - Experiments (relation → Experiments)
     - Datasets (relation → Experimental Data Assets)
     - Signatures (relation → Signatures)
     - Notes (text)
   - **🧪 HTS Campaigns Database**
     - Campaign Name (title)
     - Program (relation → Pipeline/Programs)
     - Assay Name (text)
     - Target (text)
     - Library (text or relation)
     - Date (date)
     - Molecules Screened (number)
     - Hit Cutoff (text/number)
     - Hits Found (number)
     - QC Summary (text)
     - Data Files (files or links)
     - Notes (text)
   - **🧪 Biochemical Hits Database** (Optional)
     - Compound (relation → Compound Features)
     - Program (relation → Pipeline/Programs)
     - Assay Name (text)
     - IC50 / EC50 (number)
     - Curve Quality (number/text)
     - Decision (select: Promote, Drop, Investigate)
     - Notes (text)

2. **Pathways Database** (for Pathway Analysis)
   - Create new database
   - Properties: Name, Pathway ID, Type, Source, Relations to features
   - See: `FEATURE_ENHANCEMENT_TOP_3_DETAILED.md` Section 3.3

3. **QC Properties** (for Quality Control)
   - Add "QC Metrics" (JSON/rich text)
   - Add "QC Status" (select)

4. **Signature Quality** (for Validation)
   - Add "Quality Score" (number)
   - Add "Validation Status" (select)

5. **Ingestion Status** (for Operations)
   - Add "Ingestion Status" (select)
   - Verify "Last Ingested" exists

6. **Dataset Relations** (for Comparison)
   - Add "Similar Datasets" (relation)
   - Optional: "Dataset Comparison Score" (number)

---

## ✅ Code Quality & Refactoring (Ongoing)

### Shared Utilities Consolidation

**ChatGPT Roadmap Alignment**: Section 6.1

**Refactoring Opportunities**:
- Normalization helpers → `amprenta_rag/utils/normalization.py`
- Pinecone upsert code → `amprenta_rag/ingestion/embedding_utils.py`
- Notion CRUD patterns → Already exists, enhance
- RAG chunking → Already exists, enhance

**Estimated Effort**: Ongoing, incremental

---

## 📝 Notes for Implementation

1. **Maintain Patterns**: Follow existing ingestion module patterns
2. **Idempotency**: All operations must be idempotent
3. **Error Handling**: Robust, non-blocking
4. **Logging**: Consistent prefixes (`[INGEST]`, `[FEATURE]`, `[RAG]`)
5. **Batching**: Use for Pinecone upserts
6. **Environment Config**: Respect DB IDs from config
7. **Notion Schema**: Only modify when explicitly instructed

---

## 🚀 Next Steps

1. **Review this unified roadmap**
2. **Prioritize based on immediate needs**
3. **Start with Tier 1 items (Performance & Operations)**
4. **Coordinate Notion agent for schema changes**
5. **Implement incrementally with testing**

---

**This roadmap is a living document. Update as priorities evolve and features are completed.**

