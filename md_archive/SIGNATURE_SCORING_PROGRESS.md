# Signature Scoring Implementation Progress

**Date**: December 2, 2025  
**Status**: ⏳ **Framework Complete, Integration In Progress**

---

## ✅ **COMPLETED WORK**

### 1. Configuration Support ✅
**File**: `amprenta_rag/config.py`

Added configuration options:
- `SIGNATURE_OVERLAP_THRESHOLD` (default: 0.3)
- `ENABLE_SIGNATURE_SCORING` (default: true)
- `ENABLE_LIPID_MAPPING` (default: true)
- Integrated into `PipelineConfig` dataclass

### 2. Core Module Created ✅
**File**: `amprenta_rag/ingestion/signature_matching.py` (~480 lines)

**Completed Functions**:
- ✅ `map_raw_lipid_to_canonical_species()` - Maps vendor formats to canonical names
- ✅ `fetch_all_signatures_from_notion()` - Fetches all signature pages
- ✅ `load_signature_from_notion_page()` - **COMPLETED** - Loads signature with components
- ✅ `score_signature_against_dataset()` - Wrapper around existing scoring engine
- ✅ `find_matching_signatures_for_dataset()` - Finds matches above threshold
- ✅ `update_dataset_with_signature_matches()` - Notion writeback framework

**Signature Loading Implementation**:
- Queries Signature Components DB for components linked to signature
- Parses component properties (Name, Direction, Weight)
- Builds Signature objects from components
- Handles Notion property value mapping (Up→↑, Down→↓, etc.)

---

## ⏳ **REMAINING WORK**

### 3. Dataset Integration ⏳
**File**: `amprenta_rag/ingestion/dataset_ingestion.py`

**Tasks**:
- [ ] Extract dataset species from mwTab/features
- [ ] Call `find_matching_signatures_for_dataset()` after Pinecone upsert
- [ ] Call `update_dataset_with_signature_matches()` with results
- [ ] Non-blocking error handling

**Estimated**: ~50 lines of integration code

### 4. Notion Writebacks ⏳
**File**: `amprenta_rag/ingestion/signature_matching.py`

**Tasks**:
- [ ] Verify exact Notion property names match schema
- [ ] Handle missing properties gracefully
- [ ] Format summary text nicely
- [ ] Ensure idempotency

**Status**: Framework exists, needs testing/refinement

### 5. RAG Query Extension ⏳
**Files**: 
- `amprenta_rag/query/rag_engine.py`
- `scripts/rag_query.py`

**Tasks**:
- [ ] Add `signature_similarity_query()` function
- [ ] Add `--signature-score` CLI argument
- [ ] Format output for signature rankings

**Estimated**: ~100 lines of code

### 6. Lipid Mapping Enhancement ⏳
**File**: `amprenta_rag/ingestion/signature_matching.py`

**Tasks**:
- [ ] Enhance vendor format handling
- [ ] Add Lipid Species DB lookup option
- [ ] Improve normalization rules

**Status**: Basic implementation exists, can be enhanced

---

## 🎯 **NEXT STEPS**

### Immediate Priority:
1. **Dataset Integration** (Critical)
   - This enables automatic signature matching during ingestion
   - ~50 lines of code to add
   - High impact

2. **Test End-to-End** (Verification)
   - Run dataset ingestion with signature matching
   - Verify Notion updates
   - Check logs

### Follow-up:
3. **RAG Query Extension** (User-facing feature)
   - Enables signature similarity queries
   - ~100 lines of code
   - Nice-to-have

4. **Polish & Enhance** (Quality)
   - Better error messages
   - More robust mapping
   - Summary formatting

---

## 📊 **COMPLETION ESTIMATE**

| Component | Status | Lines | Complexity |
|-----------|--------|-------|------------|
| Configuration | ✅ Complete | ~10 | Low |
| Core Module | ✅ Complete | ~480 | High |
| Dataset Integration | ⏳ Next | ~50 | Medium |
| Notion Writebacks | ⏳ Framework | ~100 | Medium |
| RAG Queries | ⏳ Pending | ~100 | Medium |
| **Total** | **~60%** | **~740** | **Medium-High** |

---

## 🔧 **TECHNICAL NOTES**

### Signature Loading Implementation
The `load_signature_from_notion_page()` function:
1. Extracts signature name and description from page properties
2. Queries Signature Components DB with filter: `{"property": "Signature", "relation": {"contains": signature_page_id}}`
3. Paginates through all components
4. Parses each component:
   - "Component Name" (title) → species
   - "Direction" (select) → direction (mapped: Up→↑, Down→↓, etc.)
   - "Weight" (number) → weight
5. Builds SignatureComponent objects
6. Returns Signature object

This implementation should work with the existing Notion schema.

### Next Integration Point
In `dataset_ingestion.py`, after successful Pinecone upsert and feature extraction:

```python
# Extract dataset species
dataset_species = set(extract_features_from_mwtab(mwtab_data))

# Find matching signatures
if cfg.pipeline.enable_signature_scoring:
    matches = find_matching_signatures_for_dataset(
        dataset_species=dataset_species,
        overlap_threshold=cfg.pipeline.signature_overlap_threshold,
    )
    
    if matches:
        update_dataset_with_signature_matches(
            dataset_page_id=page_id,
            matches=matches,
        )
```

---

## ⚠️ **KNOWN ISSUES / CONSIDERATIONS**

1. **Property Name Verification**: Need to verify exact Notion property names match the schema
   - "Related Lipid Signatures" (relation)
   - "Signature Overlap Summary" (rich_text)
   - "Signature Match Score" (number)

2. **Component Loading Performance**: Loading all signatures and their components on each dataset ingestion could be slow. Consider:
   - Caching signatures
   - Loading only active signatures
   - Background processing

3. **Error Handling**: Framework has try/except blocks, but needs testing to ensure graceful degradation

---

## 📝 **FILES TO MODIFY**

**Next**:
- `amprenta_rag/ingestion/dataset_ingestion.py` - Add matching call

**Later**:
- `amprenta_rag/query/rag_engine.py` - Add signature queries
- `scripts/rag_query.py` - Add CLI argument

---

**Recommendation**: Proceed with dataset integration next, then test end-to-end before adding RAG queries.

