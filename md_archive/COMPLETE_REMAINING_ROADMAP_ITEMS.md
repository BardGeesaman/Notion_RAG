# Complete Remaining Roadmap Items

**Last Updated**: 2025-01-XX  
**Status**: Comprehensive list of all remaining work

## Overview

This document provides a complete list of all remaining roadmap items across all tiers, including items that may have been missed in initial assessments.

---

## ⏳ TIER 4: QUALITY & OPERATIONS (5 Items Remaining)

### 3.1 Quality Control Extraction ❌ NOT IMPLEMENTED

**Priority**: 🔥 MEDIUM  
**Effort**: 3-4 days per omics type

**What's Needed**:
- Extract QC metrics from files/metadata:
  - RNA-seq: RIN, mapping rate, library complexity
  - Proteomics: identification rates, LFQ CVs
  - Lipidomics: TIC normalization, missingness rates
  - Metabolomics: RT drift, baseline intensity
- Store QC values in Postgres (Dataset model)
- Add QC properties to Postgres schema

**Files to Create/Update**:
- `amprenta_rag/ingestion/qc_extraction.py` (new)
- Enhance each omics ingestion module
- Add QC fields to `amprenta_rag/database/models.py`

**Postgres Schema Changes**:
- Add `qc_metrics` (JSONB) to Dataset model
- Add `qc_status` (String) to Dataset model

---

### 3.2 Signature Validation & Quality Metrics ❌ NOT IMPLEMENTED

**Priority**: 🔥 MEDIUM  
**Effort**: 4-5 days

**What's Needed**:
- Compute signature quality metrics:
  - Coverage (how many datasets match)
  - Specificity (false positive rate)
  - Reproducibility (consistency across datasets)
- Validate signatures against known ground truth
- Cross-validation scoring
- Bootstrap-based confidence scores
- Empirical p-values for signature scores
- Suggest signature improvements

**Files to Create**:
- `amprenta_rag/signatures/signature_validation.py` (new)

**Postgres Schema Changes**:
- Add `quality_score` (Float) to Signature model
- Add `validation_status` (String) to Signature model

---

### 3.3 Cross-Feature Mapping ❌ NOT IMPLEMENTED

**Priority**: 🔥 MEDIUM  
**Effort**: 3-4 days

**What's Needed**:
- Gene ↔ Protein mapping (HGNC/UniProt)
- Gene ↔ Metabolite mapping (pathway-level)
- Metabolite ↔ Lipid mapping (shared synthesis pathways)
- Store cross-feature relationships in Postgres

**Files to Create**:
- `amprenta_rag/analysis/cross_feature_mapping.py` (new)

**Postgres Schema Changes**:
- Add `related_features` (JSONB or junction table) to Feature model
- Or create `feature_relationships` table

---

### 3.7 Sync Back-Pressure Indicators ❌ NOT IMPLEMENTED

**Priority**: 🔥 LOW-MEDIUM  
**Effort**: 1 day

**What's Needed**:
- Add "Ingestion Status" property to Postgres Dataset model
- Track ingestion progress (Pending, In Progress, Complete, Failed)
- Update status during ingestion

**Files to Update**:
- `amprenta_rag/database/models.py` - Add `ingestion_status` field
- Update ingestion modules to set status

**Postgres Schema Changes**:
- Add `ingestion_status` (String) to Dataset model
- Verify `last_ingested` (DateTime) exists

---

### 3.8 Auto-Link Experiments ↔ Programs ❌ NOT IMPLEMENTED

**Priority**: 🔥 LOW-MEDIUM  
**Effort**: 2-3 days

**What's Needed**:
- Infer Program/Experiment from metadata during ingestion
- Auto-link when unambiguous
- Log when auto-linking occurs
- Manual override always available

**Files to Create/Update**:
- `amprenta_rag/ingestion/auto_linking.py` (new)
- Enhance ingestion modules to use auto-linking

---

## 📊 TIER 4: DASHBOARD & VISUALIZATION (2 Items Remaining)

### 4.1 Multi-Omics Coverage Maps ❌ NOT IMPLEMENTED

**Priority**: 💡 LOW  
**Effort**: 2-3 days

**What's Needed**:
- Visual representation of omics data across programs/datasets
- Show which omics types are available for each program
- Cross-omics integration visualization
- Generate visualizations for Streamlit dashboard

**Files to Create**:
- `amprenta_rag/analysis/coverage_maps.py` (new)
- Add visualization page to Streamlit dashboard

---

### 4.2 Feature Recurrence Visualization ❌ NOT IMPLEMENTED

**Priority**: 💡 LOW  
**Effort**: 3-4 days

**What's Needed**:
- Most recurrent genes/proteins/metabolites/lipids
- Cross-omics recurrence matrices
- Feature ↔ Omics heatmaps
- Generate during ingestion or on-demand

**Files to Create**:
- `amprenta_rag/analysis/feature_recurrence.py` (new)
- Add visualization page to Streamlit dashboard

---

## 🚨 CRITICAL PRIORITY: AI Agent Tools UI Integration

**Status**: ⚠️ **CRITICAL - MUST HAVE**  
**Priority**: 🔥🔥🔥 HIGHEST  
**Effort**: 1-2 weeks

### What's Needed

**Phase 1: Streamlit Dashboard Integration (IMMEDIATE)**

1. **RAG Query Interface Page** ❌ NOT IMPLEMENTED
   - Natural language search with AI-synthesized answers
   - Filter sidebar (source type, disease, target, lipid, signature)
   - Results display with relevance scores and citations
   - **Effort**: 2-3 days

2. **Cross-Omics Reasoning Widgets** ❌ NOT IMPLEMENTED
   - "Ask AI" buttons on Program, Signature, Feature, Dataset pages
   - AI-powered summaries using existing cross-omics functions
   - **Effort**: 2-3 days

3. **Signature Matching Interface** ❌ NOT IMPLEMENTED
   - "Find Matching Signatures" on Dataset pages
   - Visual similarity scoring
   - **Effort**: 1-2 days

**Files to Create/Update**:
- `scripts/run_dashboard.py` - Add AI tools pages
- `amprenta_rag/api/routers/ai.py` - FastAPI endpoints for AI tools (new)
- `amprenta_rag/api/services/ai.py` - AI service layer (new)

**FastAPI Endpoints Needed**:
- `POST /api/v1/ai/query` - RAG query endpoint
- `POST /api/v1/ai/cross-omics/program/{id}` - Program summary
- `POST /api/v1/ai/cross-omics/signature/{id}` - Signature summary
- `POST /api/v1/ai/cross-omics/feature` - Feature summary
- `POST /api/v1/ai/cross-omics/dataset/{id}` - Dataset summary
- `GET /api/v1/ai/signatures/match/{dataset_id}` - Signature matching

**Documentation**: `docs/AI_AGENT_TOOLS_UI_PLAN.md`

---

## 🔧 Code Quality & Refactoring (Ongoing)

### Shared Utilities Consolidation ⏳ ONGOING

**Status**: Partially complete, needs verification

**Refactoring Opportunities**:
- ✅ Normalization helpers → Already consolidated in `amprenta_rag/ingestion/features/normalization.py`
- ⚠️ Pinecone upsert code → May need consolidation to `amprenta_rag/ingestion/embedding_utils.py`
- ✅ Notion CRUD patterns → Already exists
- ⚠️ RAG chunking → May need enhancement/consolidation

**Files to Verify**:
- Check if Pinecone upsert code is duplicated across modules
- Check if RAG chunking utilities need consolidation

**Estimated Effort**: Ongoing, incremental (1-2 days to verify and complete)

---

## 🌐 TIER 2: Public Repository Ingestion (Verification Needed)

### Status: ⚠️ PARTIALLY COMPLETE

**What Exists**:
- ✅ MW (Metabolomics Workbench) - Complete
- ⚠️ GEO (transcriptomics) - Implementation exists, needs verification
- ⚠️ PRIDE (proteomics) - Implementation exists, needs verification
- ⚠️ MetaboLights (metabolomics) - Implementation exists, needs verification

**What's Needed**:
- Verify implementations work correctly
- Test with real data
- Fix any bugs or gaps
- Complete documentation

**Files to Verify**:
- `amprenta_rag/ingestion/repositories/geo.py`
- `amprenta_rag/ingestion/repositories/pride.py`
- `amprenta_rag/ingestion/repositories/metabolights.py`

**Effort**: 2-3 days for verification and fixes

---

## 🏗️ TIER 3: Architecture Evolution (Future)

### Phase 9: Architecture Evolution ⏳ FUTURE

**Status**: Marked as "Future Evolution" in roadmap

**Note**: The timeline shows Phase 9 as future, but TIER 3 is actually COMPLETE. This appears to be a planning artifact. Phase 9 items are already done:
- ✅ Domain Model Extraction
- ✅ Postgres Schema Design
- ✅ FastAPI Service Layer
- ✅ Migration Utilities
- ✅ RAG Integration with Postgres
- ✅ Transition to Postgres SoT

**No work needed here** - this is already complete.

---

## 🔮 TIER 5: STRETCH GOALS (Future)

### 5.1 Multi-Omics Embedding Fusion ❌ NOT IMPLEMENTED

**Priority**: 💡 FUTURE  
**Effort**: High (research + implementation)

**What's Needed**:
- Multi-modal embeddings for datasets
- Cross-omics similarity search
- Signature projection in embedding spaces

---

### 5.2 Pathway-Level Models ❌ NOT IMPLEMENTED

**Priority**: 💡 FUTURE  
**Effort**: High (research + implementation)

**What's Needed**:
- Pathway signatures
- Mechanistic hypotheses
- Graph-guided RAG prompts

---

### 5.3 Longitudinal Omics Integration ❌ NOT IMPLEMENTED

**Priority**: 💡 FUTURE  
**Effort**: High (new capabilities)

**What's Needed**:
- Time-series transcriptomics
- Proteomics time courses
- Dynamic modeling

---

## 📋 Summary by Priority

### 🔥🔥🔥 CRITICAL (Do First)

1. **AI Agent Tools UI Integration** (1-2 weeks)
   - RAG Query Interface
   - Cross-Omics Reasoning Widgets
   - Signature Matching Interface
   - FastAPI AI endpoints

### 🔥🔥 HIGH PRIORITY (TIER 4 Quality)

2. **Quality Control Extraction** (3-4 days per omics type)
3. **Signature Validation & Quality Metrics** (4-5 days)
4. **Cross-Feature Mapping** (3-4 days)
5. **Auto-Link Experiments ↔ Programs** (2-3 days)
6. **Sync Back-Pressure Indicators** (1 day)

### 🔥 MEDIUM PRIORITY

7. **Public Repository Ingestion Verification** (2-3 days)
8. **Multi-Omics Coverage Maps** (2-3 days)
9. **Feature Recurrence Visualization** (3-4 days)
10. **Code Refactoring Verification** (1-2 days)

### 💡 LOW PRIORITY / FUTURE

11. **TIER 5 Stretch Goals** (Future research projects)

---

## 📊 Total Remaining Work

| Category | Items | Total Effort |
|----------|-------|--------------|
| **CRITICAL** | 1 (AI Tools UI) | 1-2 weeks |
| **TIER 4 Quality** | 5 | ~13-20 days |
| **TIER 4 Dashboard** | 2 | ~5-7 days |
| **Verification** | 1 (Repos) | 2-3 days |
| **Refactoring** | 1 (Verify) | 1-2 days |
| **TIER 5** | 3 | Future |
| **TOTAL** | **13 items** | **~4-6 weeks** |

---

## 🎯 Recommended Order

1. **Week 1-2**: AI Agent Tools UI Integration (CRITICAL)
2. **Week 3**: Sync Back-Pressure Indicators + Auto-Link (quick wins)
3. **Week 4-5**: Quality Control Extraction (start with one omics type)
4. **Week 6**: Signature Validation & Quality Metrics
5. **Week 7**: Cross-Feature Mapping
6. **Week 8**: Dashboard Visualizations (Coverage Maps, Feature Recurrence)
7. **Week 9**: Public Repository Verification + Code Refactoring

---

## Notes

- The roadmap timeline shows some items as complete (✅) but status documents say they're not implemented. Trust the status documents (`COMPREHENSIVE_ROADMAP_STATUS.md`, `ACTUAL_STATUS.md`) as they reflect actual implementation.
- AI Agent Tools UI is marked as CRITICAL in `AI_AGENT_TOOLS_UI_PLAN.md` - this should be highest priority.
- Code refactoring is marked as "ongoing" - needs verification to see what's actually left.

