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

## 💡 TIER 3: QUALITY & OPERATIONS (Ongoing)

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

---

## 🎯 Priority Matrix

| Enhancement | Value | Effort | Dependencies | Priority |
|------------|-------|--------|--------------|----------|
| Feature Caching | ⭐⭐⭐⭐⭐ | Medium | None | 🔥🔥🔥 |
| Batch Ingestion | ⭐⭐⭐⭐ | Low-Medium | None | 🔥🔥🔥 |
| Signature Discovery | ⭐⭐⭐⭐⭐ | High | Stats libs | 🔥🔥 |
| Pathway Analysis | ⭐⭐⭐⭐⭐ | High | APIs + Notion | 🔥🔥 |
| Evidence Reports | ⭐⭐⭐⭐ | Medium | Cross-omics | 🔥🔥 |
| Program Maps | ⭐⭐⭐⭐ | Medium | Scoring | 🔥🔥 |
| QC Extraction | ⭐⭐⭐ | Medium | Per omics | 🔥 |
| Feature Lookup Cache | ⭐⭐⭐ | Low | None | 🔥 |

---

## 📋 Notion Agent Instructions Summary

### **Required Schema Changes**:

1. **Pathways Database** (for Pathway Analysis)
   - Create new database
   - Properties: Name, Pathway ID, Type, Source, Relations to features
   - See: `FEATURE_ENHANCEMENT_TOP_3_DETAILED.md` Section 3.3

2. **QC Properties** (for Quality Control)
   - Add "QC Metrics" (JSON/rich text)
   - Add "QC Status" (select)

3. **Signature Quality** (for Validation)
   - Add "Quality Score" (number)
   - Add "Validation Status" (select)

4. **Ingestion Status** (for Operations)
   - Add "Ingestion Status" (select)
   - Verify "Last Ingested" exists

5. **Dataset Relations** (for Comparison)
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

