# Implementation Plan - Signature Scoring & Dataset Matching

**Date**: December 2, 2025  
**Status**: In Progress

---

## 📋 **TASK BREAKDOWN**

This is a **large implementation** requiring multiple new capabilities. Breaking it into phases:

### Phase 1: Core Infrastructure ✅ (In Progress)
- [x] Configuration support added
- [x] Signature matching module structure created
- [ ] Implement signature loading from Notion
- [ ] Complete lipid mapping utilities

### Phase 2: Dataset Integration ⏳
- [ ] Integrate matching into dataset ingestion
- [ ] Extract dataset species from mwTab/features
- [ ] Run signature scoring during ingestion

### Phase 3: Notion Writebacks ⏳
- [ ] Update dataset pages with signature matches
- [ ] Add relation properties
- [ ] Add summary text properties
- [ ] Add score properties

### Phase 4: RAG Integration ⏳
- [ ] Extend RAG queries for signature scoring
- [ ] Add CLI arguments
- [ ] Implement ranking by score

---

## ⚠️ **COMPLEXITY NOTES**

### Signature Loading from Notion
**Challenge**: Signatures in Notion have components in a separate database

**Approach**:
1. Fetch signature pages from Lipid Signatures DB
2. For each signature, query Signature Components DB by "Signature" relation
3. Build Signature objects from components

**Alternative (Practical)**:
- Load signatures from SIGNATURES_DIR files (we know this works)
- Can enhance to Notion loading later

### Dataset Species Extraction
**Challenge**: Need to extract species from datasets consistently

**Approach**:
1. Use existing `extract_features_from_mwtab()` for mwTab datasets
2. Map metabolite features to lipid species using normalization
3. Handle raw lipid names with mapping utilities

---

## 🎯 **RECOMMENDATION**

Given the scope, implement in this order:

1. **Quick Win**: Use file-based signature loading (SIGNATURES_DIR)
   - Load signatures from TSV files
   - Score against datasets
   - Write back results

2. **Enhancement**: Add Notion-based signature loading
   - Query Signature Components DB
   - Build Signature objects dynamically

3. **Complete**: Full Notion integration
   - Writebacks
   - RAG queries
   - Summary generation

---

## 📝 **CURRENT STATUS**

**Files Created**:
- ✅ `signature_matching.py` - Module structure (needs completion)

**Files Modified**:
- ✅ `config.py` - Added configuration options

**Next Steps**:
1. Complete signature loading from Notion (or files)
2. Integrate into dataset ingestion
3. Add Notion writebacks
4. Extend RAG queries

---

**Estimated Completion**: This is a substantial feature. Should proceed incrementally.

