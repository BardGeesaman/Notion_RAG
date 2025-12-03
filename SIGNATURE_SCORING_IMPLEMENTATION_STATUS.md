# Signature Scoring Implementation Status

**Date**: December 2, 2025  
**Status**: ⏳ **IN PROGRESS - Framework Created**

---

## ✅ **COMPLETED**

1. **Module Structure Created** (`signature_matching.py`)
   - ✅ Basic module structure
   - ✅ Dataclass definitions (SignatureMatchResult)
   - ✅ Helper function stubs
   - ⚠️ Signature loading from Notion needs completion

2. **Configuration Added** (`config.py`)
   - ✅ `SIGNATURE_OVERLAP_THRESHOLD` (default: 0.3)
   - ✅ `ENABLE_SIGNATURE_SCORING` (default: true)
   - ✅ `ENABLE_LIPID_MAPPING` (default: true)
   - ✅ Added to `PipelineConfig` dataclass

3. **Core Functions Started**
   - ✅ `map_raw_lipid_to_canonical_species()` - Basic implementation
   - ✅ `fetch_all_signatures_from_notion()` - Implemented
   - ⚠️ `load_signature_from_notion_page()` - **Needs completion**
   - ✅ `score_signature_against_dataset()` - Wrapper implemented
   - ✅ `find_matching_signatures_for_dataset()` - Framework ready
   - ✅ `update_dataset_with_signature_matches()` - Framework ready

---

## ⏳ **IN PROGRESS**

### 1. Signature Loading from Notion
**Status**: Needs implementation

**Challenge**: Signatures in Notion have components in a separate database (Signature Components DB). Need to:
1. Query Signature Components DB for components where "Signature" relation contains signature_page_id
2. Build SignatureComponent objects from component pages
3. Build Signature object from components

**Current**: Placeholder function exists but returns None

**Next**: Implement component fetching and Signature object construction

---

## 📋 **REMAINING TASKS**

### Phase 1: Complete Core Functions
1. ✅ Implement `load_signature_from_notion_page()`
   - Fetch components from Signature Components DB
   - Parse component properties (Name, Direction, Weight)
   - Build Signature object

2. ✅ Complete lipid mapping
   - Enhance `map_raw_lipid_to_canonical_species()` 
   - Add vendor format handling
   - Add Lipid Species DB lookup

### Phase 2: Dataset Integration
3. ✅ Integrate into `dataset_ingestion.py`
   - Extract dataset species from mwTab/features
   - Call `find_matching_signatures_for_dataset()`
   - Non-blocking error handling

4. ✅ Extract dataset species
   - Use `extract_features_from_mwtab()`
   - Map metabolite features to lipid species
   - Handle raw lipid names

### Phase 3: Notion Writebacks
5. ✅ Complete `update_dataset_with_signature_matches()`
   - Fix property names to match schema
   - Handle missing properties gracefully
   - Idempotent updates

6. ✅ Add summary generation
   - Format match summary text
   - Include scores and components
   - Append to existing summaries

### Phase 4: RAG Integration
7. ✅ Extend `rag_engine.py`
   - Add `signature_similarity_query()` function
   - Rank signatures by score
   - Return structured results

8. ✅ Update `rag_query.py` CLI
   - Add `--signature-score` argument
   - Format output nicely

---

## 🎯 **RECOMMENDATION**

Given the substantial scope, recommend implementing in this order:

1. **Complete signature loading from Notion** (critical blocker)
   - This is needed for matching to work
   - Can start with file-based loading as fallback

2. **Complete dataset integration** (core functionality)
   - Extract species from datasets
   - Run matching during ingestion
   - Basic Notion writebacks

3. **Enhance and polish** (nice-to-have)
   - Better lipid mapping
   - RAG queries
   - Summary generation

---

## 📝 **FILES CREATED/MODIFIED**

**Created**:
- ✅ `amprenta_rag/ingestion/signature_matching.py` (396 lines, needs completion)
- ✅ `IMPLEMENTATION_PLAN_SIGNATURE_SCORING.md`
- ✅ `SIGNATURE_SCORING_IMPLEMENTATION_STATUS.md`

**Modified**:
- ✅ `amprenta_rag/config.py` (added config options)

**To Modify**:
- ⏳ `amprenta_rag/ingestion/dataset_ingestion.py` (add matching call)
- ⏳ `amprenta_rag/query/rag_engine.py` (add signature queries)
- ⏳ `scripts/rag_query.py` (add CLI argument)

---

## 🔧 **TECHNICAL NOTES**

### Signature Loading Approach
**Option A**: Load from Notion (requires querying components)
- More dynamic
- Always up-to-date
- Requires complex queries

**Option B**: Load from files in SIGNATURES_DIR
- Simpler to implement
- Requires file-based signatures
- Already works (we know ingestion works)

**Recommendation**: Implement Option A (Notion), with Option B as fallback

### Component Query Pattern
```python
# Query Signature Components DB
filter = {
    "property": "Signature",
    "relation": {"contains": signature_page_id}
}
```

Then parse each component page:
- "Component Name" (title) → species
- "Direction" (select) → direction (Up/Down/NoChange)
- "Weight" (number) → weight

---

## ⚠️ **BLOCKERS**

1. **Signature Loading**: Need to implement Notion component fetching
2. **Property Names**: Need to verify exact Notion property names match
3. **Dataset Species Extraction**: Need consistent extraction from mwTab/features

---

**Next Step**: Complete signature loading from Notion, then proceed with dataset integration.

