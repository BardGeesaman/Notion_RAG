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

---

## 🎯 IMMEDIATE NEXT STEPS (Priority Order)

### 1. Integration Enhancements

**Public Repositories**:
- [x] ArrayExpress integration
- [ ] Add more repository types (others)
- [ ] Automated discovery workflows
- [ ] Scheduled harvesting

**Security & Collaboration**:
- [ ] Multi-User Support (RBAC) - Added to roadmap Tier 3.2
- [ ] Team-based data access
- [ ] Audit logging

### 2. ELN & Workflow Features

**Protocol Management**:
- [ ] Structured protocol templates
- [ ] Version control for protocols
- [ ] Link protocols to experiments

**Sample Inventory**:
- [ ] Physical location tracking (freezers/boxes)
- [ ] Sample lineage and parent/child tracking
- [ ] Barcode generation support

**Compliance & Integrity**:
- [ ] 21 CFR Part 11 readiness (audit trails)
- [ ] Electronic signatures (future)
- [ ] Data immutability logs

### 3. Chemistry Informatics Enhancements

**Chemical Registration**:
- [ ] Corporate ID generation (e.g. AMP-001)
- [ ] Salt and batch handling
- [ ] Duplicate checking logic

**Cheminformatics Search**:
- [ ] Substructure search (SMARTS)
- [ ] Similarity search (Tanimoto)
- [ ] Pharmacophore search

**SAR Analysis**:
- [ ] Calculated properties (LogP, MW, TPSA)
- [ ] Activity cliffs visualization
- [ ] R-group decomposition

### 4. RAG Maturity & Enhancements

**Retrieval Optimization**:
- [ ] Hybrid Search (Sparse/BM25 + Dense Vectors)
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
2. Security & Multi-User Support

### Phase 2: ELN & Workflow Features
1. Protocol Management & Templates
2. Sample Inventory System
3. Audit Trails & Compliance

### Phase 3: Chemistry Informatics
1. Chemical Registration System
2. Advanced Structure Search
3. SAR Analysis Tools

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

**Current System Status**: ✅ Production Ready (ELN Integration Complete)  
**Architecture**: Postgres-only (Notion removed December 2025)  
**Next Milestone**: Integration Enhancements  
**Estimated Time**: 3-4 weeks

