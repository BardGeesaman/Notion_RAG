# Optional Enhancements Implementation Plan

## Overview
This document outlines the implementation plan for the optional enhancements:
1. Signature Systems Migration (Postgres-first)
2. LLM-Based Semantic Metadata Extraction

---

## Phase 1: Signature Systems Migration

### Step 1: Database Schema Updates ✅ (Done)
- ✅ Added fields to Signature model:
  - `short_id` (String) - Short identifier
  - `biomarker_role` (Array[String]) - Biomarker roles
  - `phenotype_axes` (Array[String]) - Phenotype axes
  - `data_ownership` (String) - Data ownership classification

### Step 2: Create Migration (Pending)
- Create Alembic migration for Signature fields
- Apply migration: `alembic upgrade head`

### Step 3: Postgres-Based Signature Metadata Extraction
Create new module: `amprenta_rag/ingestion/metadata/postgres_signature_metadata.py`

**Functions to create:**
- `collect_signature_metadata_from_postgres(signature_ids: List[UUID]) -> Dict[str, Any]`
  - Fetch signatures from Postgres by IDs
  - Extract: short_id, biomarker_role, phenotype_axes, data_ownership
  - Return structured metadata dictionary

- `get_signature_metadata_from_postgres(signature: SignatureModel) -> Dict[str, Any]`
  - Extract metadata from a single Signature model instance

### Step 4: Postgres-Based Feature Extraction for Datasets
Create new module: `amprenta_rag/ingestion/postgres_feature_extraction.py`

**Functions to create:**
- `extract_dataset_features_by_type_from_postgres(dataset_id: UUID, omics_type: Optional[str] = None) -> Dict[str, Set[str]]`
  - Extract features from Postgres dataset
  - Group by feature type (gene, protein, metabolite, lipid)
  - Use existing `dataset.features` relationship

### Step 5: Update Signature Matching to Use Postgres
Update: `amprenta_rag/ingestion/signature_matching/matching.py`

**Changes:**
- Add parameter `dataset_id: Optional[UUID]` alongside `dataset_page_id`
- Use `extract_dataset_features_by_type_from_postgres()` when `dataset_id` is provided
- Load signatures from Postgres instead of Notion
- Fallback to Notion if Postgres dataset_id not available

### Step 6: Update Signature Detection to Use Postgres
Update: `amprenta_rag/ingestion/signature_integration.py`

**Changes:**
- Add parameter `source_dataset_id: Optional[UUID]` alongside `source_page_id`
- Store detected signatures in Postgres
- Link signatures to datasets in Postgres using `dataset_signature_assoc` table
- Update `link_signature_to_source()` to work with Postgres

### Step 7: Postgres-Based Signature Linking
Create new module: `amprenta_rag/ingestion/postgres_signature_linking.py`

**Functions to create:**
- `link_signature_to_dataset_in_postgres(signature_id: UUID, dataset_id: UUID, match_score: Optional[float] = None) -> None`
  - Link signature to dataset using `dataset_signature_assoc` table
  - Store match score if provided

- `get_dataset_signatures_from_postgres(dataset_id: UUID) -> List[SignatureModel]`
  - Get all signatures linked to a dataset

- `get_signature_datasets_from_postgres(signature_id: UUID) -> List[DatasetModel]`
  - Get all datasets linked to a signature

### Step 8: Update Dataset Ingestion
Update: `amprenta_rag/ingestion/postgres_dataset_ingestion.py`

**Changes:**
- Pass `dataset_id` instead of `notion_page_id` to signature functions
- Use Postgres signature linking after matching
- Store signature matches in Postgres

---

## Phase 2: LLM-Based Semantic Metadata Extraction

### Step 1: Create LLM Extraction Module
Create new module: `amprenta_rag/ingestion/metadata/llm_semantic_extraction.py`

**Functions to create:**
- `extract_semantic_metadata_with_llm(text: str, source_type: str) -> Dict[str, Any]`
  - Use OpenAI API to extract structured metadata
  - Extract: diseases, targets, signatures, phenotype axes, biomarker roles
  - Return structured dictionary

- `enhance_metadata_with_llm(metadata: Dict[str, Any], text: str) -> Dict[str, Any]`
  - Enhance existing metadata with LLM extraction
  - Merge results intelligently

### Step 2: Update Semantic Extraction Module
Update: `amprenta_rag/ingestion/metadata/postgres_semantic_extraction.py`

**Changes:**
- Add optional LLM-based extraction
- Use pattern matching as fallback
- Configurable LLM usage via config flag

### Step 3: Integration
Update ingestion functions to optionally use LLM extraction:
- `postgres_dataset_ingestion.py`
- `postgres_experiment_ingestion.py`
- `postgres_content_ingestion.py`

---

## Implementation Order

1. ✅ Add Signature fields to model
2. ⏳ Create and apply Signature migration
3. ⏳ Create Postgres signature metadata extraction
4. ⏳ Create Postgres feature extraction
5. ⏳ Update signature matching
6. ⏳ Update signature detection
7. ⏳ Create Postgres signature linking
8. ⏳ Update dataset ingestion
9. ⏳ Create LLM extraction module
10. ⏳ Integrate LLM extraction

---

## Files to Create

1. `amprenta_rag/ingestion/metadata/postgres_signature_metadata.py`
2. `amprenta_rag/ingestion/postgres_feature_extraction.py`
3. `amprenta_rag/ingestion/postgres_signature_linking.py`
4. `amprenta_rag/ingestion/metadata/llm_semantic_extraction.py`
5. `alembic/versions/XXXX_add_signature_metadata_fields.py`

## Files to Modify

1. `amprenta_rag/ingestion/signature_matching/matching.py`
2. `amprenta_rag/ingestion/signature_integration.py`
3. `amprenta_rag/ingestion/postgres_dataset_ingestion.py`
4. `amprenta_rag/ingestion/metadata/postgres_semantic_extraction.py`

---

## Testing Strategy

1. Test signature metadata extraction from Postgres
2. Test feature extraction from Postgres datasets
3. Test signature matching with Postgres dataset_id
4. Test signature detection with Postgres dataset_id
5. Test signature linking in Postgres
6. Test LLM semantic extraction (with rate limiting)

---

## Configuration

Add to config:
- `ENABLE_LLM_SEMANTIC_EXTRACTION: bool = False` (default: False)
- `ENABLE_POSTGRES_SIGNATURES: bool = True` (default: True)

