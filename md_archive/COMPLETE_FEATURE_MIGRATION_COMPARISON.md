# Complete Feature Migration Comparison

## Overview
This document compares ALL functionality between the old Notion-based ingestion system and the new Postgres-only system to ensure nothing is missing.

---

## ✅ Dataset Ingestion

### Core Features

| Feature | Notion-Based | Postgres-Only | Status |
|---------|-------------|---------------|--------|
| **Text Content Extraction** | ✅ Extract from Notion blocks | ✅ Build from Postgres fields | ✅ Complete |
| **Text Chunking** | ✅ Uses `chunk_text()` | ✅ Uses `chunk_text()` | ✅ Complete |
| **Text Embedding** | ✅ Uses `embed_texts()` | ✅ Uses `embed_texts()` | ✅ Complete |
| **Pinecone Upsert** | ✅ Batch upsert with metadata | ✅ Batch upsert with metadata | ✅ Complete |
| **mwTab Data Extraction** | ✅ From page content OR MW API | ✅ From MW API (via external_ids) | ✅ Complete |
| **mwTab JSON in Text** | ✅ Included in text content | ✅ Included in text content | ✅ Complete |

### Metadata Extraction

| Feature | Notion-Based | Postgres-Only | Status |
|---------|-------------|---------------|--------|
| **Basic Metadata** | ✅ Title, omics_type | ✅ Name, omics_type | ✅ Complete |
| **Disease** | ✅ From "Disease" multi-select | ✅ From `dataset.disease` array | ✅ Complete |
| **Matrix** | ✅ From "Matrix" multi-select | ✅ From `dataset.sample_type` array | ✅ Complete |
| **Model Systems** | ✅ From "Model Systems" multi-select | ✅ From `dataset.organism` array | ✅ Complete |
| **External IDs** | ✅ From properties | ✅ From `dataset.external_ids` JSONB | ✅ Complete |
| **DOI** | ✅ From properties | ✅ From `dataset.external_ids['doi']` | ✅ Complete |
| **Study IDs** | ✅ From properties | ✅ From `dataset.external_ids` | ✅ Complete |
| **File URLs** | ✅ From properties | ✅ From `dataset.file_urls` array | ✅ Complete |
| **Dataset Source Type** | ✅ From "Dataset Source Type" select | ⚠️ **MISSING** | ⚠️ Gap |
| **Data Origin** | ✅ From "Data Origin" select | ⚠️ **MISSING** | ⚠️ Gap |

### Semantic Metadata (Rich Metadata)

| Feature | Notion-Based | Postgres-Only | Status |
|---------|-------------|---------------|--------|
| **Lipid Signatures (Relation)** | ✅ From "Related Signature(s)" relation | ⚠️ **MISSING** - Requires Notion page ID | ⚠️ Gap |
| **Lipid Signature Short IDs** | ✅ Extracted from signature pages | ⚠️ **MISSING** | ⚠️ Gap |
| **Lipid Signature Role** | ✅ From signature pages | ⚠️ **MISSING** | ⚠️ Gap |
| **Phenotype Axes** | ✅ From signature pages | ⚠️ **MISSING** | ⚠️ Gap |
| **Signature Ownership** | ✅ From signature pages | ⚠️ **MISSING** | ⚠️ Gap |
| **Program Links** | ✅ From "Related Programs" relation | ✅ From `dataset.programs` relationship | ✅ Complete |
| **Experiment Links** | ✅ From "Related Experiments" relation | ✅ From `dataset.experiments` relationship | ✅ Complete |

### Feature Linking

| Feature | Notion-Based | Postgres-Only | Status |
|---------|-------------|---------------|--------|
| **Feature Extraction from mwTab** | ✅ `extract_features_from_mwtab()` | ✅ `extract_features_from_mwtab()` | ✅ Complete |
| **Feature Linking** | ✅ `link_features_to_notion_items()` | ✅ `batch_link_features_to_dataset_in_postgres()` | ✅ Complete |
| **Multi-Omics Support** | ✅ Genes, proteins, metabolites, lipids | ✅ Genes, proteins, metabolites, lipids | ✅ Complete |
| **Feature Normalization** | ✅ Omics-specific normalization | ✅ Omics-specific normalization | ✅ Complete |

### Scientific Metadata (from mwTab)

| Feature | Notion-Based | Postgres-Only | Status |
|---------|-------------|---------------|--------|
| **Model Systems from mwTab** | ✅ Extracted and updated to Notion | ⚠️ **NOT EXTRACTED** - Only from Postgres fields | ⚠️ Gap |
| **Disease Terms from mwTab** | ✅ Extracted and updated to Notion | ⚠️ **NOT EXTRACTED** - Only from Postgres fields | ⚠️ Gap |
| **Matrix Terms from mwTab** | ✅ Extracted and updated to Notion | ⚠️ **NOT EXTRACTED** - Only from Postgres fields | ⚠️ Gap |
| **Methods from mwTab** | ✅ Extracted and updated to Notion | ⚠️ **NOT EXTRACTED** | ⚠️ Gap |
| **Summary from mwTab** | ✅ Extracted and updated to Notion | ⚠️ **NOT EXTRACTED** | ⚠️ Gap |
| **Results from mwTab** | ✅ Extracted and updated to Notion | ⚠️ **NOT EXTRACTED** | ⚠️ Gap |
| **Conclusions from mwTab** | ✅ Extracted and updated to Notion | ⚠️ **NOT EXTRACTED** | ⚠️ Gap |
| **Data Origin from mwTab** | ✅ Extracted and updated to Notion | ⚠️ **NOT EXTRACTED** | ⚠️ Gap |
| **Source URL from mwTab** | ✅ Extracted and updated to Notion | ⚠️ **NOT EXTRACTED** | ⚠️ Gap |

### Signature Systems

| Feature | Notion-Based | Postgres-Only | Status |
|---------|-------------|---------------|--------|
| **Signature Matching** | ✅ Matches dataset species against signatures | ✅ Matches dataset species (if Notion page ID available) | ⚠️ Partial |
| **Signature Detection** | ✅ Detects new signatures from content | ✅ Detects new signatures (if Notion page ID available) | ⚠️ Partial |
| **Signature Updates to Notion** | ✅ Updates dataset page with matches | ⚠️ **MISSING** - No Postgres storage for matches | ⚠️ Gap |
| **Multi-Omics Signature Matching** | ✅ Supports all omics types | ✅ Supports all omics types | ✅ Complete |

### Notion Integration

| Feature | Notion-Based | Postgres-Only | Status |
|---------|-------------|---------------|--------|
| **Embedding Metadata Update** | ✅ Updates "Embedding IDs" and "Last Embedded" | ⚠️ **OPTIONAL** - Only if `update_notion=True` | ⚠️ Optional |
| **Scientific Metadata Update** | ✅ Updates Notion page with mwTab metadata | ⚠️ **OPTIONAL** - Only if `update_notion=True` | ⚠️ Optional |
| **Force Re-ingestion** | ✅ Checks "Last Embedded" timestamp | ✅ No timestamp check (force flag) | ⚠️ Different |

---

## ✅ Experiment Ingestion

### Core Features

| Feature | Notion-Based | Postgres-Only | Status |
|---------|-------------|---------------|--------|
| **Text Content Extraction** | ✅ Extract from Notion blocks + properties | ✅ Build from Postgres fields | ✅ Complete |
| **Text Chunking** | ✅ Uses `chunk_text()` | ✅ Uses `chunk_text()` | ✅ Complete |
| **Text Embedding** | ✅ Uses `embed_texts()` | ✅ Uses `embed_texts()` | ✅ Complete |
| **Pinecone Upsert** | ✅ Batch upsert with metadata | ✅ Batch upsert with metadata | ✅ Complete |
| **Linked Datasets** | ✅ Included in text content | ✅ Included in text content | ✅ Complete |
| **Linked Programs** | ✅ Included in text content | ✅ Included in text content | ✅ Complete |

### Metadata Extraction

| Feature | Notion-Based | Postgres-Only | Status |
|---------|-------------|---------------|--------|
| **Basic Metadata** | ✅ Name, Type | ✅ Name, Type | ✅ Complete |
| **Disease** | ✅ From "Disease" multi-select | ✅ From `experiment.disease` array | ✅ Complete |
| **Matrix** | ✅ From "Matrix" multi-select | ✅ From `experiment.matrix` array | ✅ Complete |
| **Model Systems** | ✅ From "Model Systems" multi-select | ✅ From `experiment.model_systems` array | ✅ Complete |
| **Type** | ✅ From "Type" select | ✅ From `experiment.type` | ✅ Complete |
| **Targets** | ✅ From "Targets" multi-select | ⚠️ **MISSING** - Not in Postgres model | ⚠️ Gap |
| **Modality** | ✅ From "Modality" multi-select | ⚠️ **MISSING** - Not in Postgres model | ⚠️ Gap |
| **Stage** | ✅ From "Stage" select | ⚠️ **MISSING** - Not in Postgres model | ⚠️ Gap |
| **Biomarker Role** | ✅ From "Biomarker Role" multi-select | ⚠️ **MISSING** - Not in Postgres model | ⚠️ Gap |
| **Treatment Arms** | ✅ From "Treatment Arms" multi-select | ⚠️ **MISSING** - Not in Postgres model | ⚠️ Gap |

### Semantic Metadata (Rich Metadata)

| Feature | Notion-Based | Postgres-Only | Status |
|---------|-------------|---------------|--------|
| **Lipid Signatures (Relation)** | ✅ From "Lipid Signatures" relation | ⚠️ **MISSING** - Requires Notion page ID | ⚠️ Gap |
| **Lipid Signature Short IDs** | ✅ Extracted from signature pages | ⚠️ **MISSING** | ⚠️ Gap |
| **Lipid Signature Role** | ✅ From signature pages | ⚠️ **MISSING** | ⚠️ Gap |
| **Phenotype Axes** | ✅ From experiment + signature pages | ⚠️ **MISSING** | ⚠️ Gap |
| **Signature Ownership** | ✅ From signature pages | ⚠️ **MISSING** | ⚠️ Gap |
| **Program Links** | ✅ From "Related Programs" relation | ✅ From `experiment.programs` relationship | ✅ Complete |

### Feature Linking

| Feature | Notion-Based | Postgres-Only | Status |
|---------|-------------|---------------|--------|
| **Feature Extraction from Text** | ✅ `extract_features_from_text()` | ✅ `extract_features_from_text()` | ✅ Complete |
| **Feature Linking to Datasets** | ✅ Links to related datasets | ✅ Links to related datasets | ✅ Complete |
| **Multi-Omics Support** | ✅ All feature types | ✅ All feature types | ✅ Complete |

### Signature Systems

| Feature | Notion-Based | Postgres-Only | Status |
|---------|-------------|---------------|--------|
| **Signature Detection** | ✅ Detects new signatures from content | ✅ Detects new signatures (if Notion page ID available) | ⚠️ Partial |

---

## ✅ Email Ingestion

### Core Features

| Feature | Notion-Based | Postgres-Only | Status |
|---------|-------------|---------------|--------|
| **Text Content Extraction** | ✅ Extract from Notion blocks | ✅ Direct email content | ✅ Complete |
| **Text Chunking** | ✅ Uses `chunk_text()` | ✅ Uses `chunk_text()` | ✅ Complete |
| **Text Embedding** | ✅ Uses `embed_texts()` | ✅ Uses `embed_texts()` | ✅ Complete |
| **Pinecone Upsert** | ✅ Batch upsert with metadata | ✅ Batch upsert with metadata | ✅ Complete |
| **Email Metadata** | ✅ Title, From, Tags | ✅ Title, From, Tags | ✅ Complete |
| **Idempotency** | ✅ Checks "Embedding Status" | ✅ Content hash check | ✅ Complete (Better!) |

### Metadata Extraction

| Feature | Notion-Based | Postgres-Only | Status |
|---------|-------------|---------------|--------|
| **Basic Metadata** | ✅ Title, From, Tags | ✅ Title, From, Tags | ✅ Complete |
| **Email ID** | ✅ From Notion page ID | ✅ From `email_id` parameter | ✅ Complete |
| **Semantic Metadata** | ✅ `get_email_semantic_metadata()` | ⚠️ **MISSING** - No semantic extraction | ⚠️ Gap |

### Semantic Metadata (Rich Metadata)

| Feature | Notion-Based | Postgres-Only | Status |
|---------|-------------|---------------|--------|
| **Disease** | ✅ Extracted from email content | ⚠️ **MISSING** | ⚠️ Gap |
| **Targets** | ✅ Extracted from email content | ⚠️ **MISSING** | ⚠️ Gap |
| **Lipid Signatures** | ✅ Extracted from email content | ⚠️ **MISSING** | ⚠️ Gap |
| **Phenotype Axes** | ✅ Extracted from email content | ⚠️ **MISSING** | ⚠️ Gap |

### Signature Systems

| Feature | Notion-Based | Postgres-Only | Status |
|---------|-------------|---------------|--------|
| **Signature Detection** | ✅ Detects new signatures from content | ⚠️ **SKIPPED** - No Notion page ID | ⚠️ Gap |

### Notion Integration

| Feature | Notion-Based | Postgres-Only | Status |
|---------|-------------|---------------|--------|
| **RAG Chunk Pages** | ✅ Creates chunk pages in Notion | ⚠️ **NOT CREATED** - Direct to Pinecone only | ⚠️ Different |
| **Email Status Update** | ✅ Updates "Embedding Status" to "Embedded" | ⚠️ **NOT UPDATED** - No Notion storage | ⚠️ Different |

---

## ✅ Literature/Zotero Ingestion

### Core Features

| Feature | Notion-Based | Postgres-Only | Status |
|---------|-------------|---------------|--------|
| **Text Content Extraction** | ✅ Extract from Zotero + Notion blocks | ✅ Extract from Zotero + PDFs | ✅ Complete |
| **Text Chunking** | ✅ Uses `chunk_text()` | ✅ Uses `chunk_text()` | ✅ Complete |
| **Text Embedding** | ✅ Uses `embed_texts()` | ✅ Uses `embed_texts()` | ✅ Complete |
| **Pinecone Upsert** | ✅ Batch upsert with metadata | ✅ Batch upsert with metadata | ✅ Complete |
| **PDF Attachment Processing** | ✅ Extract text from PDFs | ✅ Extract text from PDFs | ✅ Complete |
| **Notes Processing** | ✅ Extract HTML notes | ✅ Extract HTML notes | ✅ Complete |
| **Idempotency** | ✅ Checks Notion page existence | ✅ Content hash check | ✅ Complete (Better!) |

### Metadata Extraction

| Feature | Notion-Based | Postgres-Only | Status |
|---------|-------------|---------------|--------|
| **Basic Metadata** | ✅ Title, Authors, DOI, Zotero Key | ✅ Title, Authors, DOI, Zotero Key | ✅ Complete |
| **Abstract** | ✅ From Zotero item | ✅ From Zotero item | ✅ Complete |
| **Semantic Metadata** | ✅ `get_literature_semantic_metadata()` | ⚠️ **MISSING** - No semantic extraction | ⚠️ Gap |

### Semantic Metadata (Rich Metadata)

| Feature | Notion-Based | Postgres-Only | Status |
|---------|-------------|---------------|--------|
| **Disease** | ✅ Extracted from literature content | ⚠️ **MISSING** | ⚠️ Gap |
| **Targets** | ✅ Extracted from literature content | ⚠️ **MISSING** | ⚠️ Gap |
| **Lipid Signatures** | ✅ Extracted from literature content | ⚠️ **MISSING** | ⚠️ Gap |
| **Phenotype Axes** | ✅ Extracted from literature content | ⚠️ **MISSING** | ⚠️ Gap |

### Signature Systems

| Feature | Notion-Based | Postgres-Only | Status |
|---------|-------------|---------------|--------|
| **Signature Detection** | ✅ Detects new signatures from content | ⚠️ **SKIPPED** - No Notion page ID | ⚠️ Gap |

### Notion Integration

| Feature | Notion-Based | Postgres-Only | Status |
|---------|-------------|---------------|--------|
| **Notion Page Creation** | ✅ Creates literature page in Notion | ⚠️ **NOT CREATED** - Direct to Pinecone only | ⚠️ Different |
| **RAG Chunk Pages** | ✅ Creates chunk pages in Notion | ⚠️ **NOT CREATED** - Direct to Pinecone only | ⚠️ Different |

---

## 🔍 Summary of Gaps

### Critical Gaps (Features that should be migrated)

1. **Scientific Metadata Extraction from mwTab** (Datasets)
   - Model systems, disease terms, matrix terms from mwTab
   - Methods, summary, results, conclusions from mwTab
   - Currently only extracted to Notion, not stored in Postgres

2. **Semantic Metadata Extraction** (All types)
   - Disease, targets, lipid signatures from content
   - Currently uses Notion-based extraction functions
   - Need Postgres-compatible extraction

3. **Signature Metadata from Relations** (Datasets & Experiments)
   - Lipid signature short IDs, roles, axes, ownership
   - Currently requires Notion page IDs to fetch signature pages
   - Need Postgres-based signature linking

4. **Experiment Rich Metadata** (Experiments)
   - Targets, Modality, Stage, Biomarker Role, Treatment Arms
   - Not currently in Postgres Experiment model

5. **Dataset Rich Metadata** (Datasets)
   - Dataset Source Type, Data Origin
   - Not currently in Postgres Dataset model

### Optional Gaps (Features that may not be needed)

1. **Notion Page Creation/Updates**
   - RAG chunk pages in Notion
   - Email/literature pages in Notion
   - Status updates to Notion
   - These are optional features for documentation/UI

2. **Signature Detection/Matching**
   - Currently requires Notion page IDs
   - Can be migrated to Postgres-based signatures later

---

## 🎯 Recommendations

### High Priority (Core Functionality)

1. **Add Scientific Metadata Extraction to Postgres**
   - Extract mwTab metadata and store in Postgres Dataset model
   - Update `ingest_dataset_from_postgres()` to extract and store

2. **Migrate Semantic Metadata Extraction**
   - Create Postgres-compatible semantic metadata extraction
   - Extract from text content directly (no Notion dependency)

3. **Add Missing Experiment Fields**
   - Add Targets, Modality, Stage, Biomarker Role, Treatment Arms to Experiment model
   - Update ingestion to extract and store these fields

4. **Add Missing Dataset Fields**
   - Add Dataset Source Type, Data Origin to Dataset model
   - Update ingestion to extract and store these fields

### Medium Priority (Enhanced Features)

1. **Migrate Signature Systems**
   - Store signatures in Postgres
   - Update signature matching/detection to use Postgres

2. **Content Hash-Based Idempotency**
   - Already implemented for emails/literature ✅
   - Add to datasets/experiments for better incremental updates

### Low Priority (Documentation/UI)

1. **Notion Integration** (Optional)
   - Keep as optional feature for backward compatibility
   - Can be enabled with `ENABLE_NOTION_SYNC=true`

---

## ✅ What's Already Complete

- ✅ Core ingestion pipeline (text extraction, chunking, embedding, Pinecone upsert)
- ✅ Basic metadata extraction (disease, matrix, model systems, external IDs)
- ✅ Feature linking (all omics types, Postgres-first)
- ✅ Program/Experiment linking
- ✅ mwTab data fetching and inclusion in text
- ✅ Idempotency for emails/literature (content hash)
- ✅ Direct-to-Pinecone ingestion (faster, no Notion overhead)
- ✅ Gmail integration (replaces Zapier)
- ✅ Zotero integration (direct API access)

---

## 📋 Next Steps

1. **Review this comparison** to confirm which gaps are critical
2. **Prioritize migrations** based on usage/importance
3. **Implement high-priority gaps** first
4. **Test thoroughly** to ensure feature parity
5. **Update documentation** once complete

