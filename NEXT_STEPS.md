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

## 🎯 IMMEDIATE NEXT STEPS (Priority Order)

### 1. Experimental Design Metadata System (5-6 days) 🔥🔥🔥

**Priority**: HIGH - Enables proper comparative analysis

**What's Needed**:

**Phase 1: Database Schema (1-2 days)**
- Extend experiment/dataset models with design fields
- Add `design_type` (case_control, time_course, intervention, dose_response, multi_factorial)
- Add `design_metadata` (JSONB for sample groups, timepoints, interventions)
- Add `sample_group`, `timepoint`, `intervention` to datasets table

**Phase 2: Repository Metadata Extraction (2-3 days)**
- Parse GEO sample groups from SOFT files
- Parse Metabolomics Workbench study design
- Implement pattern matching for common designs (case/control detection)
- Add LLM-based extraction as fallback for complex designs
- Auto-detect timepoint labels in sample names

**Phase 3: Design-Aware Statistical Analysis (2 days)**
- Implement design-aware t-tests and ANOVA
- Add repeated measures analysis for time courses
- Add multi-group comparisons for intervention studies
- Add dose-response trend tests
- Integrate with existing statistical analysis modules

**Files to Create/Update**:
- `amprenta_rag/ingestion/design_extraction.py` (new)
- `amprenta_rag/analysis/design_aware_stats.py` (new)
- `amprenta_rag/models/experiment.py` (extend)
- `amprenta_rag/models/dataset.py` (extend)
- `amprenta_rag/ingestion/repositories/geo.py` (update)
- `amprenta_rag/ingestion/repositories/mw.py` (update)

**Impact**: Enables proper case/control comparisons, time course analysis, intervention group tracking

---

### 2. Chemistry & HTS Integration (3-4 days) 🧪

**Goal**: Ingest High-Throughput Screening data (up to 1M molecules)

**What's Needed**:
- Implement SQLite layer for chemistry data
- Build screening ingestion pipeline
- Link compounds to biological signatures they reverse

**Impact**: Expands platform from "biology-only" to "drug discovery"

---

### 3. Enhancements & Polish (Various) 🔥

**Evidence Reports**:
- [ ] PDF export functionality
- [ ] HTML report generation
- [ ] Email distribution

**Visualizations**:
- [ ] Program-signature map heatmaps
- [ ] Dataset similarity dendrograms
- [ ] Pathway enrichment bar charts
- [ ] Cross-omics convergence plots

**Performance**:
- [ ] Further optimize batch operations
- [ ] Parallel processing for large datasets

**RAG Integration**:
- [ ] Embed promoted compounds into Pinecone
- [ ] Add compound summaries to RAG index
- [ ] Enable compound queries in RAG

---

## 📊 MEDIUM PRIORITY (Future Enhancements)

### 4. Advanced Analytics

**Signature Discovery Enhancements**:
- [ ] Multi-dataset pattern detection
- [ ] Temporal signature evolution
- [ ] Signature validation metrics

**Dataset Analysis**:
- [ ] Differential expression analysis
- [ ] Time-series analysis
- [ ] Cohort comparison

### 5. User Experience

**CLI Improvements**:
- [ ] Interactive mode for complex operations
- [ ] Progress bars for long operations
- [ ] Better error messages and recovery

**Documentation**:
- [ ] API documentation (Sphinx/ReadTheDocs)
- [ ] Tutorial notebooks
- [ ] Video walkthroughs

### 6. Integration Enhancements

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

1. **Add PDF Export to Evidence Reports** (2-3 hours)
   - Use `reportlab` or `weasyprint`
   - Simple implementation

2. **Create Visualization Helpers** (1 day)
   - Matplotlib/Plotly wrappers
   - Standard chart templates

3. **Add Progress Bars** (2-3 hours)
   - Use `tqdm` for long operations
   - Better user feedback

4. **Improve Error Messages** (2-3 hours)
   - More descriptive errors
   - Recovery suggestions

---

## 📋 RECOMMENDED WORKFLOW

### Week 1: Experimental Design System
1. Extend database schema (design_type, sample_groups, timepoints)
2. Implement repository metadata extraction (GEO, MW)
3. Add design-aware statistical tests

### Week 2: Chemistry Integration
1. SQLite layer for HTS data
2. Screening ingestion pipeline
3. Compound-signature linking

### Week 3: Testing & Polish
1. Comprehensive testing
2. Fix any issues
3. Add quick wins (PDF export, progress bars)

### Week 4: Enhancements
1. Visualizations
2. RAG integration for compounds
3. Advanced analytics

---

## 🎯 IMMEDIATE ACTION ITEMS

**Today/Tomorrow**:
1. Start Experimental Design Metadata schema design
2. Review GEO/MW metadata extraction patterns

**This Week**:
1. Complete Experimental Design Phase 1 (schema)
2. Begin Phase 2 (metadata extraction)

**Next Week**:
1. Complete Experimental Design system
2. Begin Chemistry/HTS integration

---

## 💡 SUGGESTIONS

**For Maximum Impact**:
- Focus on Experimental Design (enables proper statistical comparisons)
- Then Chemistry integration (expands to drug discovery)
- Then enhancements (improves usability)

**For Quick Wins**:
- PDF export (highly requested)
- Progress bars (better UX)
- Better error messages (reduces support burden)

---

## 📞 QUESTIONS TO CONSIDER

1. **Priority**: What's most important right now?
   - Experimental Design metadata?
   - Chemistry/HTS integration?
   - Quick wins and polish?

2. **Resources**: Do you have access to:
   - Example GEO/MW datasets with complex designs?
   - HTS screening data for testing?
   - User feedback for prioritization?

3. **Timeline**: What's the target for:
   - Production deployment?
   - Feature completion?
   - User rollout?

---

**Current System Status**: ✅ Production Ready (Pathway Analysis Complete)  
**Architecture**: Postgres-only (Notion removed December 2025)  
**Next Milestone**: Experimental Design Metadata System  
**Estimated Time**: 5-6 days for experimental design

