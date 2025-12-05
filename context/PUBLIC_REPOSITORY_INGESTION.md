# Public Repository Ingestion for Multi-Omics Data

**Status**: 📋 Roadmap Item - Tier 2.7 (Promoted from Tier 5.4)

**Priority**: 🔥🔥 **MEDIUM-HIGH**

**Strategic Value**:
- Expand dataset coverage across all omics types
- Leverage public data repositories
- Automate discovery and ingestion
- Build comprehensive multi-omics knowledge base

---

## 0. CURRENT STATE

### Existing Implementation

**Metabolomics Workbench (MW)**:
- ✅ `scripts/harvest_mw_studies.py` - Harvests MW studies
- ✅ `scripts/discover_lipidomics_studies.py` - Discovers lipidomics studies
- ✅ Supports lipidomics and metabolomics data
- ✅ Creates Notion pages with mwTab data
- ✅ Triggers dataset ingestion pipeline

**Pattern to Extend**:
1. Discover studies by keyword/metadata
2. Fetch study metadata and data
3. Create/update Notion Dataset pages
4. Trigger appropriate ingestion pipeline (lipidomics, metabolomics, etc.)

---

## 1. TARGET REPOSITORIES

### 1.1 Transcriptomics: Gene Expression Omnibus (GEO)

**Repository**: https://www.ncbi.nlm.nih.gov/geo/

**API**: 
- Entrez E-utilities API
- GEO Datasets API
- Series (GSE) and Samples (GSM) identifiers

**Data Types**:
- Microarray data
- RNA-seq data
- Single-cell RNA-seq
- DGE tables

**Key Endpoints**:
- `https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi` - Search studies
- `https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi` - Fetch metadata
- `https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc={GSE_ID}&targ=self&form=text&view=full` - Full metadata

**Implementation**:
- Search by keywords (disease, organism, platform)
- Fetch GSE metadata
- Download processed data (if available)
- Create transcriptomics dataset pages
- Trigger `ingest_transcriptomics_file`

**Estimated Effort**: Medium (4-5 days)

---

### 1.2 Proteomics: PRIDE Archive

**Repository**: https://www.ebi.ac.uk/pride/archive

**API**: 
- PRIDE REST API
- PRIDE Archive API v2

**Data Types**:
- Mass spectrometry proteomics
- Protein identification
- Quantification data (LFQ, iBAQ, etc.)

**Key Endpoints**:
- `https://www.ebi.ac.uk/pride/archive/projects` - List projects
- `https://www.ebi.ac.uk/pride/archive/projects/{project_id}` - Project details
- `https://www.ebi.ac.uk/pride/archive/projects/{project_id}/files` - File listing

**Implementation**:
- Search by keywords (disease, organism, instrument)
- Fetch project metadata
- Download processed protein tables (if available)
- Create proteomics dataset pages
- Trigger `ingest_proteomics_file`

**Estimated Effort**: Medium (4-5 days)

---

### 1.3 Metabolomics: MetaboLights

**Repository**: https://www.ebi.ac.uk/metabolights/

**API**: 
- MetaboLights REST API

**Data Types**:
- Metabolomics data
- NMR and MS data
- Metabolite identification

**Key Endpoints**:
- `https://www.ebi.ac.uk/metabolights/webservice/study/list` - List studies
- `https://www.ebi.ac.uk/metabolights/webservice/study/{study_id}` - Study details

**Implementation**:
- Search by keywords
- Fetch study metadata
- Download metabolite tables
- Create metabolomics dataset pages
- Trigger `ingest_metabolomics_file`

**Estimated Effort**: Medium (3-4 days)

---

### 1.4 Extended Metabolomics Workbench

**Current**: Only lipidomics-focused discovery

**Enhancement**:
- Extend to full metabolomics (not just lipidomics)
- Better keyword filtering
- Support for polar metabolomics studies

**Estimated Effort**: Low (1-2 days)

---

## 2. UNIFIED ARCHITECTURE

### 2.1 Repository Abstraction Layer

Create a unified interface for all repositories:

**Location**: `amprenta_rag/ingestion/repositories/`

**Structure**:
```
amprenta_rag/ingestion/repositories/
├── __init__.py
├── base.py              # Base repository interface
├── geo.py              # GEO (transcriptomics)
├── pride.py            # PRIDE (proteomics)
├── metabolights.py     # MetaboLights (metabolomics)
├── mw.py               # Metabolomics Workbench (refactored)
└── discovery.py        # Unified discovery across repositories
```

**Base Interface**:
```python
class RepositoryInterface:
    def search_studies(self, keywords: List[str], filters: Dict) -> List[StudyMetadata]
    def fetch_study_metadata(self, study_id: str) -> StudyMetadata
    def fetch_study_data(self, study_id: str) -> Optional[DataFile]
    def get_omics_type(self) -> str  # "transcriptomics", "proteomics", etc.
```

### 2.2 Study Metadata Model

**Location**: `amprenta_rag/models/repository.py`

```python
@dataclass
class StudyMetadata:
    study_id: str
    repository: str  # "GEO", "PRIDE", "MetaboLights", "MW"
    title: str
    summary: str
    doi: Optional[str]
    disease: List[str]
    organism: List[str]
    sample_type: List[str]
    omics_type: str
    data_files: List[DataFile]
    publication_date: Optional[datetime]
```

### 2.3 Unified Discovery Script

**Location**: `scripts/discover_omics_studies.py`

**Features**:
- Search across all repositories
- Filter by omics type, disease, organism
- Unified output format
- Option to harvest discovered studies

**Usage**:
```bash
# Discover transcriptomics studies
python scripts/discover_omics_studies.py \
  --omics-type transcriptomics \
  --keywords "ALS" "Alzheimer" \
  --repository GEO

# Discover proteomics studies
python scripts/discover_omics_studies.py \
  --omics-type proteomics \
  --keywords "cancer" \
  --repository PRIDE

# Discover across all repositories
python scripts/discover_omics_studies.py \
  --keywords "disease" \
  --all-repositories
```

---

## 3. IMPLEMENTATION PHASES

### Phase 1: Repository Abstraction (Week 1)

**Tasks**:
- Create base repository interface
- Refactor MW ingestion to use new interface
- Create StudyMetadata model
- Create unified discovery script

**Estimated Effort**: Medium (3-4 days)

---

### Phase 2: GEO (Transcriptomics) Integration (Week 2)

**Tasks**:
- Implement GEO repository interface
- GEO API client
- Study discovery and metadata fetching
- Data file download (processed DGE tables)
- Integration with transcriptomics ingestion

**Estimated Effort**: Medium (4-5 days)

---

### Phase 3: PRIDE (Proteomics) Integration (Week 3)

**Tasks**:
- Implement PRIDE repository interface
- PRIDE API client
- Study discovery and metadata fetching
- Data file download (protein tables)
- Integration with proteomics ingestion

**Estimated Effort**: Medium (4-5 days)

---

### Phase 4: MetaboLights Integration (Week 4)

**Tasks**:
- Implement MetaboLights repository interface
- MetaboLights API client
- Study discovery and metadata fetching
- Data file download (metabolite tables)
- Integration with metabolomics ingestion

**Estimated Effort**: Medium (3-4 days)

---

### Phase 5: Extended MW Support (Week 4)

**Tasks**:
- Extend MW discovery beyond lipidomics
- Support full metabolomics studies
- Better keyword filtering

**Estimated Effort**: Low (1-2 days)

---

**Total Estimated Effort**: 15-20 days (Weeks 1-4)

---

## 4. INTEGRATION WITH EXISTING PIPELINES

### 4.1 Workflow

1. **Discovery**: `discover_omics_studies.py` finds relevant studies
2. **Harvest**: Repository-specific harvest script creates Notion pages
3. **Ingestion**: Existing ingestion pipeline processes data files
4. **Feature Linking**: Automatic feature linking (already implemented)
5. **Signature Scoring**: Automatic signature matching (already implemented)
6. **RAG Embedding**: Automatic RAG embedding (already implemented)

### 4.2 Repository-Specific Harvest Scripts

**Pattern** (similar to `harvest_mw_studies.py`):
- `scripts/harvest_geo_studies.py` - Harvest GEO studies
- `scripts/harvest_pride_studies.py` - Harvest PRIDE studies
- `scripts/harvest_metabolights_studies.py` - Harvest MetaboLights studies

**Common Features**:
- Fetch study metadata
- Create/update Notion Dataset pages
- Download data files (if available)
- Trigger appropriate ingestion pipeline
- Idempotent operations

---

## 5. DATA QUALITY & STANDARDIZATION

### 5.1 Metadata Normalization

- Normalize disease names across repositories
- Normalize organism names
- Standardize sample type classifications
- Map repository-specific fields to Notion schema

### 5.2 Data Format Handling

- Handle different file formats (CSV, TSV, Excel, etc.)
- Parse repository-specific data structures
- Normalize to ingestion pipeline expectations
- Quality checks before ingestion

---

## 6. CONFIGURATION

### 6.1 Repository Configuration

Add to `config.py`:
```python
class RepositoryConfig:
    GEO_API_KEY: Optional[str] = None  # Optional, for higher rate limits
    PRIDE_API_BASE_URL: str = "https://www.ebi.ac.uk/pride/archive"
    METABOLIGHTS_API_BASE_URL: str = "https://www.ebi.ac.uk/metabolights/webservice"
    MW_API_BASE_URL: str = "https://www.metabolomicsworkbench.org/rest"
    
    # Rate limiting
    GEO_RATE_LIMIT: int = 3  # requests per second
    PRIDE_RATE_LIMIT: int = 10
    METABOLIGHTS_RATE_LIMIT: int = 5
```

### 6.2 Discovery Configuration

- Default keywords per omics type
- Default filters (disease, organism, etc.)
- Maximum results per repository
- Cache settings for API responses

---

## 7. TESTING STRATEGY

### 7.1 Unit Tests

- Test repository API clients
- Test metadata extraction
- Test data file parsing
- Test Notion page creation

### 7.2 Integration Tests

- End-to-end discovery → harvest → ingestion
- Test with real repository data
- Verify feature linking works
- Verify signature scoring works

### 7.3 Sample Studies

- Identify test studies from each repository
- Use for validation and testing
- Document expected results

---

## 8. USAGE EXAMPLES

### 8.1 Discover Transcriptomics Studies

```bash
# Discover ALS-related transcriptomics studies in GEO
python scripts/discover_omics_studies.py \
  --omics-type transcriptomics \
  --repository GEO \
  --keywords "ALS" "amyotrophic lateral sclerosis" \
  --disease "ALS" \
  --organism "Homo sapiens" \
  --max-results 50
```

### 8.2 Harvest Specific Study

```bash
# Harvest a GEO study
python scripts/harvest_geo_studies.py \
  --study-id GSE12345 \
  --create-notion \
  --ingest

# Harvest a PRIDE project
python scripts/harvest_pride_studies.py \
  --project-id PXD012345 \
  --create-notion \
  --ingest
```

### 8.3 Batch Harvest

```bash
# Harvest multiple studies from discovery results
python scripts/discover_omics_studies.py \
  --omics-type transcriptomics \
  --keywords "Alzheimer" \
  --output studies.json

python scripts/harvest_geo_studies.py \
  --study-list studies.json \
  --create-notion \
  --ingest
```

---

## 9. PRIORITIZATION

### Recommended Order

1. **GEO (Transcriptomics)** - High value, well-documented API
2. **PRIDE (Proteomics)** - High value, good API support
3. **MetaboLights (Metabolomics)** - Complements MW
4. **Extended MW** - Low effort, high impact

---

## 10. DEPENDENCIES

- Existing ingestion pipelines (already implemented)
- Feature linking system (already implemented)
- Signature scoring system (already implemented)
- RAG embedding system (already implemented)
- Notion API client (already implemented)

**New Dependencies**:
- `biopython` (for GEO Entrez API)
- `requests` (already used)
- Optional: `pandas` (already used)

---

## 11. CHALLENGES & CONSIDERATIONS

### 11.1 API Rate Limits

- Implement rate limiting per repository
- Cache API responses
- Batch requests where possible
- Respect repository terms of service

### 11.2 Data Format Variability

- Different repositories use different formats
- Some studies have processed data, others raw
- Need robust parsing and normalization

### 11.3 Metadata Quality

- Repository metadata quality varies
- Some fields may be missing
- Need fallback strategies

### 11.4 Large Datasets

- Some studies have very large files
- May need streaming or chunked downloads
- Consider storage requirements

---

**Last Updated**: 2025-12-04
**Status**: Ready for Implementation

