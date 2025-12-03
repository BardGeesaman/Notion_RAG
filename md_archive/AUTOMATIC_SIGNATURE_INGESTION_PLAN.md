# Automatic Signature Ingestion Implementation Plan

## 🎯 Goal

Make lipid signature ingestion fully automatic, disease-agnostic, and integrated across ALL ingestion pipelines.

---

## 📋 Implementation Requirements

### 1. ✅ Signature Detection Module (COMPLETED)
**File**: `amprenta_rag/ingestion/signature_detection.py`

**Status**: ✅ Created with:
- Keyword detection
- Embedded table extraction
- File attachment detection
- URL detection
- Text table parsing

### 2. ⏳ Enhanced Signature Ingestion (IN PROGRESS)

**Add to `signature_ingestion.py`**:

A. **Reverse Linking Functions** (for all source types):
- `link_signature_to_source()` - Generalized reverse linking
  - Maps source_type → Notion relation property name
  - Links signature page to source page

B. **Metabolite Features Integration**:
- `link_component_to_metabolite_feature()` - Links signature components to Metabolite Features
- Ensures: Signature → Component → Lipid Species → Metabolite Feature chain

C. **Signature Embedding**:
- `embed_signature()` - Creates text representation, chunks, embeds, upserts to Pinecone
- Stores embeddings in RAG Engine DB

### 3. ⏳ Pipeline Integration (PENDING)

**Add signature detection to each ingestion module**:

**A. Literature (`zotero_ingest.py`)**:
- After Pinecone upsert, scan chunks for signature keywords
- Check for attached signature TSV/CSV files
- Extract embedded signature tables
- Call `ingest_signature_from_file()` if signature found
- Link signature back to literature page

**B. Dataset (`dataset_ingestion.py`)**:
- Check mwTab for signature-like structures
- Look for attached signature files
- Extract signature from dataset metadata if present
- Link signature back to dataset page

**C. Email (`email_ingestion.py`)**:
- Scan email body for signature keywords
- Extract embedded tables
- Link signature back to email page

**D. Experiments (`experiments_ingestion.py`)**:
- Scan experiment description for signature keywords
- Extract embedded tables
- Link signature back to experiment page

### 4. ⏳ Bulk Integration (PENDING)

**Integrate bulk ingestion into harvest scripts**:
- Add optional flag to `harvest_mw_studies.py`
- Or create unified `scripts/run_all_ingestion.py`

---

## 🔧 Technical Approach

### Detection Strategy

1. **Text-based Detection**:
   - Scan for signature keywords
   - Look for table-like structures
   - Parse embedded TSV/CSV content

2. **File-based Detection**:
   - Check attachments for signature files
   - Download from URLs if needed

3. **Metadata-based Detection**:
   - Extract from source page metadata
   - Use existing relation properties

### Integration Pattern

Each ingestion module follows this pattern:

```python
# After successful Pinecone upsert

# 1. Detect signatures
signature_files = detect_signature_files_in_content(...)
embedded_tables = extract_embedded_signature_tables(...)

# 2. Ingest each signature found
for sig_file in signature_files:
    result = ingest_signature_from_file(...)
    link_signature_to_source(result['signature_page_id'], source_page_id, source_type)

# 3. Link components to Metabolite Features
link_components_to_metabolite_features(...)
```

### Disease-Agnostic Approach

- ✅ No hardcoded disease keywords
- ✅ Extract disease from source metadata
- ✅ Infer from filename patterns
- ✅ Support any disease context

---

## 📝 Files to Create/Modify

### New Files
- ✅ `amprenta_rag/ingestion/signature_detection.py` (COMPLETED)

### Files to Enhance
- ⏳ `amprenta_rag/ingestion/signature_ingestion.py` - Add reverse linking, embedding, Metabolite Features
- ⏳ `amprenta_rag/ingestion/zotero_ingest.py` - Add signature detection
- ⏳ `amprenta_rag/ingestion/dataset_ingestion.py` - Add signature detection
- ⏳ `amprenta_rag/ingestion/email_ingestion.py` - Add signature detection
- ⏳ `amprenta_rag/ingestion/experiments_ingestion.py` - Add signature detection
- ⏳ `scripts/harvest_mw_studies.py` - Optional bulk integration

---

## ⚠️ Complexity Assessment

This is a **major enhancement** that touches:
- 5+ ingestion modules
- Multiple Notion databases
- Signature embedding pipeline
- Cross-linking between multiple systems

**Recommendation**: Implement in phases:
1. ✅ Core detection module (DONE)
2. ⏳ Reverse linking functions (IN PROGRESS)
3. ⏳ Integrate into one pipeline as proof-of-concept
4. ⏳ Roll out to all pipelines
5. ⏳ Add embedding capability
6. ⏳ Add Metabolite Features cross-linking

---

## 🎯 Next Steps

Given the complexity, should we:
- **Option A**: Implement fully across all pipelines (comprehensive)
- **Option B**: Implement core infrastructure + one pipeline as proof-of-concept
- **Option C**: Break into smaller phases with separate implementation cycles

The current state:
- ✅ Signature detection module: Complete
- ✅ Signature ingestion from files: Complete
- ✅ Bulk ingestion: Complete
- ⏳ Automatic detection in pipelines: Needs implementation
- ⏳ Reverse linking (all source types): Needs implementation
- ⏳ Signature embedding: Needs implementation
- ⏳ Metabolite Features cross-linking: Needs implementation

**What would you like me to prioritize?**

