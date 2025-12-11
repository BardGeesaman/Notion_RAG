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

### 5. Advanced Analytics ✅ COMPLETE

**Status**: Fully implemented (December 2025)

**Implemented Features**:
- ✅ Differential Expression Analysis (`differential_expression.py`)
- ✅ Time-series Analysis (`timeseries.py`)
- ✅ Cohort Comparison (`cohort_comparison.py`)

---

## 🎯 IMMEDIATE NEXT STEPS (Priority Order)

### 1. Enhancements & Polish (Remaining) 🔥

**Performance**:
- [ ] Further optimize batch operations
- [ ] Parallel processing for large datasets

**RAG Integration**:
- [ ] Embed promoted compounds into Pinecone (Enhancement)
- [ ] Add compound summaries to RAG index
- [ ] Enable compound queries in RAG

---

## 📊 MEDIUM PRIORITY (Future Enhancements)

### 3. Advanced Analytics

**Signature Discovery Enhancements**:
- [ ] Multi-dataset pattern detection
- [ ] Temporal signature evolution
- [ ] Signature validation metrics

**Dataset Analysis**:
- [ ] Differential expression analysis
- [ ] Time-series analysis
- [ ] Cohort comparison

### 4. User Experience

**CLI Improvements**:
- [ ] Interactive mode for complex operations
- [ ] Progress bars for long operations
- [ ] Better error messages and recovery

**Documentation**:
- [ ] API documentation (Sphinx/ReadTheDocs)
- [ ] Tutorial notebooks
- [ ] Video walkthroughs

### 5. Integration Enhancements

**Public Repositories**:
- [ ] Add more repository types (ArrayExpress, etc.)
- [ ] Automated discovery workflows
- [ ] Scheduled harvesting

**Security & Collaboration**:
- [ ] Multi-User Support (RBAC) - Added to roadmap Tier 3.2
- [ ] Team-based data access
- [ ] Audit logging

---

## 🚀 QUICK WINS (Can Do Now)

*All identified quick wins (PDF Export, Visualization Helpers, Progress Bars, Error Messages) have been completed.*

---

## 📋 RECOMMENDED WORKFLOW

### Phase 1: Advanced Analytics (Next)
1. Differential Expression Analysis
2. Time-series Analysis
3. Cohort Comparison

### Phase 2: Integration Enhancements
1. Public Repository Expansion
2. Security & Multi-User Support

---

## 🎯 IMMEDIATE ACTION ITEMS

**This Week**:
1. Review Advanced Analytics requirements
2. Plan Differential Expression implementation

---

## 💡 SUGGESTIONS

**For Maximum Impact**:
- Focus on Advanced Analytics (enables deeper biological insights)
- Then Integration Enhancements (expands data sources and security)
- Then RAG Integration for compounds (links chemistry to biology)

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

**Current System Status**: ✅ Production Ready (Experimental Design & UX Polish Complete)  
**Architecture**: Postgres-only (Notion removed December 2025)  
**Next Milestone**: Advanced Analytics  
**Estimated Time**: 2-3 weeks

