# Status Update for ChatGPT: Dataset Ingestion Pipeline

## ✅ What Cursor Has Already Implemented

**File: `amprenta_rag/ingestion/dataset_ingestion.py`**

### 1. Metadata Extraction (≈80% Complete)
- ✅ `_extract_mwtab_from_page_content()` - Extracts mwTab JSON from Notion page content (code blocks under "mwTab Data" heading)
- ✅ `_extract_metadata_from_mwtab()` - Parses mwTab JSON and extracts:
  - ✅ `model_systems` (from SUBJECT.SUBJECT_SPECIES)
  - ✅ `disease_terms` (pattern matching from STUDY_TITLE/STUDY_SUMMARY)
  - ✅ `matrix_terms` (from SUBJECT.SUBJECT_TYPE)
  - ✅ `methods` (from TREATMENT, SAMPLEPREP, CHROMATOGRAPHY sections)
  - ✅ `summary` (from STUDY.STUDY_SUMMARY)
  - ✅ `results` (basic metabolite count info)
  - ❌ `conclusions` - NOT YET IMPLEMENTED
  - ❌ `data_origin` - NOT YET IMPLEMENTED
  - ❌ `dataset_source_type` - NOT YET IMPLEMENTED
  - ❌ `source_url` - NOT YET IMPLEMENTED

### 2. Notion Update Helper (≈70% Complete)
- ✅ `_update_experimental_data_asset_metadata()` - Updates Notion page with metadata
- ✅ Updates: Model Systems, Disease, Matrix (all as multi_select)
- ✅ Updates: Methods, Summary, Results (all as rich_text - matches schema)
- ❌ `Conclusions` - NOT YET IMPLEMENTED
- ❌ `Data Origin` (select) - NOT YET IMPLEMENTED
- ❌ `Dataset Source Type` (select) - NOT YET IMPLEMENTED
- ❌ `Source URL / DOI` (url) - NOT YET IMPLEMENTED
- ❌ `Full Text (hidden)` - NOT YET IMPLEMENTED (optional)
- ❌ `Chunks (hidden)` - NOT YET IMPLEMENTED (optional)

### 3. Integration (✅ Complete)
- ✅ Wired into `ingest_dataset()` flow
- ✅ Runs after successful Pinecone upsert
- ✅ Runs after Embedding IDs + Last Embedded update
- ✅ Error handling: Metadata failures logged as warnings, don't fail ingestion

### 4. Backfill Script (✅ Complete)
- ✅ `scripts/backfill_dataset_metadata.py` exists and works
- ✅ Finds pages missing metadata
- ✅ Processes pages with mwTab data

---

## 📋 What's Still Needed (Per ChatGPT's Instructions)

### Missing Fields in Metadata Extraction:
1. `conclusions` - Needs LLM generation or extraction from mwTab
2. `data_origin` - Should infer from source (e.g., "External – Open Dataset" for Metabolomics Workbench)
3. `dataset_source_type` - Should infer from mwTab structure (e.g., "Processed table" for ST004396)
4. `source_url` - Extract from mwTab or pass as parameter

### Missing Fields in Notion Update:
1. `Conclusions` (rich_text)
2. `Data Origin` (select) - Values: "Internal – Amprenta", "External – Published", "External – CRO/Vendor", "External – Open Dataset", "Collaboration"
3. `Dataset Source Type` (select) - Values: "Raw file", "Processed table", "Summary (stats only)", "DOI / link only", "Supplemental table"
4. `Source URL / DOI` (url)
5. `Full Text (hidden)` (rich_text) - Optional: store full mwTab text
6. `Chunks (hidden)` (rich_text) - Optional: store chunk text

---

## ✅ Confirmed Schema Details

Verified actual Notion schema:
- `Methods` = **rich_text** ✅ (we're using this correctly)
- `Summary` = **rich_text** ✅ (we're using this correctly)
- `Results` = **rich_text** ✅ (we're using this correctly)
- `Conclusions` = **rich_text** ✅ (need to add this)

---

## 🎯 Next Steps for Cursor

To fully match ChatGPT's specification, Cursor needs to:

1. **Enhance `_extract_metadata_from_mwtab()`** to extract/derive:
   - `conclusions` (can use LLM or leave empty for now)
   - `data_origin` (infer: "External – Open Dataset" for MW studies)
   - `dataset_source_type` (infer: "Processed table" if MS_METABOLITE_DATA present)
   - `source_url` (from mwTab METABOLOMICS WORKBENCH section or pass as param)

2. **Enhance `_update_experimental_data_asset_metadata()`** to set:
   - `Conclusions` (rich_text)
   - `Data Origin` (select)
   - `Dataset Source Type` (select)
   - `Source URL / DOI` (url)
   - Optionally: `Full Text (hidden)` and `Chunks (hidden)`

3. **Optional LLM Enhancement**: Consider using LLM to generate better `conclusions` and enhance `results` summaries.

---

## 📝 Current Implementation Status

**Summary**: We have ~75% of the required functionality implemented. The core pipeline works, but we're missing:
- Conclusions field
- Data Origin / Dataset Source Type / Source URL
- Optional Full Text / Chunks fields

The foundation is solid and follows the exact patterns you requested. We just need to add the remaining fields.

---

## ✅ Verification

To verify the current implementation works:
```bash
python scripts/ingest_dataset.py --dataset-page-id "2bdadf61-42ab-811c-b2b2-cbd014210210" --force
```

This will:
- ✅ Embed to Pinecone
- ✅ Set Embedding IDs + Last Embedded
- ✅ Set Model Systems, Disease, Matrix, Methods, Summary, Results
- ❌ Missing: Conclusions, Data Origin, Dataset Source Type, Source URL

