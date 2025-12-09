# Next Steps - Implementation Roadmap

**Last Updated**: 2025-12-04  
**Status**: Core features complete, ready for enhancements

---

## 🎯 IMMEDIATE NEXT STEPS (Priority Order)

### 1. Complete Pathway Analysis ID Mapping (2-3 days) 🔥🔥🔥

**Status**: Framework complete, needs ID mapping services

**What's Needed**:
- Implement UniProt mapping for proteins
- Implement gene symbol → KEGG ID mapping
- Implement metabolite name → KEGG compound ID mapping
- Complete KEGG/Reactome API integration

**Files to Update**:
- `amprenta_rag/analysis/pathway_analysis.py` (ID mapping functions)
- Consider adding `amprenta_rag/analysis/id_mapping.py` for centralized mapping

**Resources**:
- UniProt REST API: https://www.uniprot.org/help/api
- KEGG REST API: https://www.kegg.jp/kegg/rest/keggapi.html
- Reactome API: https://reactome.org/ContentService/
- MyGene.info for gene mapping: https://docs.mygene.info/

**Impact**: Enables full pathway enrichment analysis for all omics types

---

### 2. Test All New Features with Real Data (1-2 days) 🧪

**What to Test**:
- [ ] Program signature maps with real programs
- [ ] Dataset comparison with multiple datasets
- [ ] Evidence reports for all entity types
- [ ] Chemistry compound promotion workflow
- [ ] Pathway enrichment (once ID mapping is complete)

**Test Scripts**:
- Create integration tests for each major feature
- Test with real Notion data
- Verify performance improvements

---

### 3. Pathway Analysis - Notion Integration (1-2 days) 🔥🔥

**Status**: Pathways database created, needs integration

**What's Needed**:
- Add pathway relations to feature pages
- Add pathway properties to datasets/signatures
- Create pathway summary pages
- Link enriched pathways back to Notion

**Files to Create/Update**:
- `amprenta_rag/analysis/pathway_notion_integration.py`
- Update `amprenta_rag/analysis/enrichment.py`

---

### 4. Experimental Design Metadata System (5-6 days) 🔥🔥🔥

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

### 5. Enhancements & Polish (Various) 🔥

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
- [ ] Add caching for pathway queries
- [ ] Parallel processing for large datasets

**RAG Integration**:
- [ ] Embed promoted compounds into Pinecone
- [ ] Add compound summaries to RAG index
- [ ] Enable compound queries in RAG

---

## 📊 MEDIUM PRIORITY (Future Enhancements)

### 6. Advanced Analytics

**Signature Discovery Enhancements**:
- [ ] Multi-dataset pattern detection
- [ ] Temporal signature evolution
- [ ] Signature validation metrics

**Dataset Analysis**:
- [ ] Differential expression analysis
- [ ] Time-series analysis
- [ ] Cohort comparison

### 7. User Experience

**CLI Improvements**:
- [ ] Interactive mode for complex operations
- [ ] Progress bars for long operations
- [ ] Better error messages and recovery

**Documentation**:
- [ ] API documentation (Sphinx/ReadTheDocs)
- [ ] Tutorial notebooks
- [ ] Video walkthroughs

### 8. Integration Enhancements

**Public Repositories**:
- [ ] Add more repository types (ArrayExpress, etc.)
- [ ] Automated discovery workflows
- [ ] Scheduled harvesting

**Notion Enhancements**:
- [ ] Dashboard views
- [ ] Automated report generation
- [ ] Notification system

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

### Week 1: Complete Pathway Analysis
1. Implement ID mapping services
2. Test with real data
3. Add Notion integration

### Week 2: Experimental Design System
1. Extend database schema (design_type, sample_groups, timepoints)
2. Implement repository metadata extraction (GEO, MW)
3. Add design-aware statistical tests

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
1. ✅ Complete pathway ID mapping (highest priority)
2. ✅ Test pathway enrichment with real data
3. ✅ Add pathway Notion integration

**This Week**:
1. Comprehensive feature testing
2. Fix any bugs found
3. Add PDF export to reports

**Next Week**:
1. Visualizations
2. RAG compound integration
3. Performance optimizations

---

## 💡 SUGGESTIONS

**For Maximum Impact**:
- Focus on pathway analysis completion first (unlocks full multi-omics analysis)
- Then testing (ensures reliability)
- Then enhancements (improves usability)

**For Quick Wins**:
- PDF export (highly requested)
- Progress bars (better UX)
- Better error messages (reduces support burden)

---

## 📞 QUESTIONS TO CONSIDER

1. **Priority**: What's most important right now?
   - Completing pathway analysis?
   - Testing existing features?
   - Adding new capabilities?

2. **Resources**: Do you have access to:
   - UniProt/KEGG API keys (if needed)?
   - Test datasets for validation?
   - User feedback for prioritization?

3. **Timeline**: What's the target for:
   - Production deployment?
   - Feature completion?
   - User rollout?

---

**Current System Status**: ✅ Production Ready (core features)  
**Next Milestone**: Pathway Analysis Complete  
**Estimated Time**: 2-3 days for pathway completion

