# Feature Extraction Implementation Status

## ✅ Completed

### 1. Core Module Created
- **File**: `amprenta_rag/ingestion/feature_extraction.py`
- **Functions Implemented**:
  - ✅ `normalize_metabolite_name()` - Normalizes metabolite names, handles prefixes, synonyms
  - ✅ `extract_features_from_mwtab()` - Extracts metabolite names from mwTab JSON structure
  - ✅ `extract_features_from_text()` - Lightweight text scanner for amino acids, nucleotides, ceramides, common metabolites
  - ✅ `link_features_to_notion_items()` - Creates/updates metabolite pages and adds relations
  - ✅ `_find_or_create_metabolite_page()` - Finds existing or creates new Metabolite Features pages
  - ✅ `_add_relation_to_metabolite_page()` - Adds relations from metabolite pages to target items

### 2. Configuration Updated
- **File**: `amprenta_rag/config.py`
- ✅ Added `NOTION_METABOLITE_FEATURES_DB_ID` environment variable support
- ✅ Added `metabolite_features_db_id` field to `NotionConfig` dataclass

### 3. Database Creation Script
- **File**: `scripts/create_metabolite_features_db.py`
- ✅ Script to create the Metabolite Features database via Notion API
- ✅ Creates database with required schema:
  - Name (title)
  - Class (select) - Amino Acid, Organic Acid, Nucleotide, Coenzyme, Sugar, Lipid, Other
  - Synonyms (rich_text)
  - Notes (rich_text)
  - Datasets (relation → Experimental Data Assets)
  - Literature Mentions (relation → Literature DB)
  - Emails / Notes (relation → Email DB)
  - Pathways (multi_select)
- ✅ Note: Experiments, Lipid Species, and Lipid Signature Components relations can be added manually once those DB IDs are available

### 4. Dataset Ingestion Integration
- **File**: `amprenta_rag/ingestion/dataset_ingestion.py`
- ✅ Imports feature extraction functions
- ✅ Extracts mwTab data (reuses existing extraction logic)
- ✅ Calls `extract_features_from_mwtab()` on mwTab JSON
- ✅ Calls `link_features_to_notion_items()` with `item_type="dataset"`
- ✅ Integrated after Pinecone upsert, non-blocking (warnings only)

## 🚧 Partially Complete

### 5. Literature Ingestion Integration
- **File**: `amprenta_rag/ingestion/zotero_ingest.py`
- ✅ Imports added
- ✅ Feature extraction code added at end of function
- ⚠️ **Issue**: Text collection logic needs refinement - currently re-extracts text which is inefficient
- **Status**: Functional but should be optimized to collect text during processing

### 6. Email/Notes Ingestion Integration
- **File**: `amprenta_rag/ingestion/email_ingestion.py`
- ⏳ **Not yet started** - ready to implement following same pattern

### 7. Experiments Ingestion Integration
- **File**: `amprenta_rag/ingestion/experiments_ingestion.py`
- ⏳ **Not yet started** - ready to implement following same pattern

## 📋 Remaining Work

### Immediate Tasks

1. **Optimize Literature Ingestion** (`zotero_ingest.py`)
   - Collect text during processing loops instead of re-extracting
   - Store `full_text` or chunks in a list as we process
   - Extract features from collected text at the end

2. **Email Ingestion Integration** (`email_ingestion.py`)
   - Add imports for feature extraction functions
   - Extract features from email body text
   - Link features after successful Pinecone upsert
   - Pattern: `extract_features_from_text(email_body)` → `link_features_to_notion_items(..., item_type="email")`

3. **Experiments Ingestion Integration** (`experiments_ingestion.py`)
   - Add imports for feature extraction functions
   - Extract features from experiment notes/descriptions
   - Link features after successful Pinecone upsert
   - Pattern: `extract_features_from_text(experiment_text)` → `link_features_to_notion_items(..., item_type="experiment")`

### Database Setup

4. **Create Metabolite Features Database**
   - Run: `python scripts/create_metabolite_features_db.py --parent-page-id <parent_page_id>`
   - Add returned database ID to `.env`: `NOTION_METABOLITE_FEATURES_DB_ID=<id>`
   - Optionally add Experiments, Lipid Species, Lipid Signature Components relations manually in Notion UI

## 🎯 Implementation Pattern

For each ingestion type, the integration follows this pattern:

```python
# 1. Add imports
from amprenta_rag.ingestion.feature_extraction import (
    extract_features_from_text,  # or extract_features_from_mwtab for datasets
    link_features_to_notion_items,
)

# 2. After successful Pinecone upsert, extract features
try:
    feature_names = extract_features_from_text(content_text)  # or mwTab for datasets
    if feature_names:
        link_features_to_notion_items(
            feature_names=feature_names,
            item_page_id=page_id,  # with dashes
            item_type="dataset",  # or "literature", "email", "experiment"
        )
except Exception as e:
    logger.warning("[INGEST][MODULE] Error extracting/linking features: %r", e)
    # Don't raise - feature extraction is non-critical
```

## ✅ Verification Checklist

After completing all integrations:

1. **Create Database**
   - [ ] Run database creation script
   - [ ] Add DB ID to `.env`
   - [ ] Verify database appears in Notion

2. **Test Dataset Ingestion** (ST004396)
   - [ ] Run: `python scripts/ingest_dataset.py --dataset-page-id 2bdadf61-42ab-811c-b2b2-cbd014210210 --force`
   - [ ] Check logs for metabolite extraction
   - [ ] Verify Metabolite Features pages created in Notion
   - [ ] Verify "Datasets" relation populated on metabolite pages
   - [ ] Verify dataset page shows reverse relation (if dual property)

3. **Test Literature Ingestion**
   - [ ] Run literature ingestion on a test item
   - [ ] Check logs for feature extraction
   - [ ] Verify metabolite pages linked to Literature page

4. **Test Email/Notes Ingestion**
   - [ ] Run email ingestion
   - [ ] Check logs for feature extraction
   - [ ] Verify metabolite pages linked

5. **Test Experiments Ingestion**
   - [ ] Run experiment ingestion
   - [ ] Check logs for feature extraction
   - [ ] Verify metabolite pages linked

## 📝 Notes

- Feature extraction is **non-critical** - failures are logged as warnings but don't block ingestion
- All feature extraction happens **after** successful Pinecone upsert
- Relations are **idempotent** - adding the same relation multiple times is safe
- Metabolite pages are created **on-demand** when first referenced
- Text extraction uses lightweight pattern matching - can be enhanced with ML/NER later

## 🔄 Next Steps

1. Complete remaining ingestion integrations (email, experiments)
2. Optimize literature ingestion text collection
3. Create and configure Metabolite Features database
4. Test end-to-end with ST004396
5. Verify metabolite pages and relations in Notion

