# Code Refactoring Analysis - 2025

**Date**: 2025-12-03

**Purpose**: Analyze large files and create a plan to break them into smaller, more maintainable modules.

---

## 📊 File Size Analysis

### Files Over 600 Lines (High Priority)

1. **`cross_omics_reasoning.py`** - 892 lines ⚠️ **LARGEST**
2. **`feature_extraction.py`** - 821 lines
3. **`signature_notion_crud.py`** - 801 lines
4. **`signature_matching.py`** - 670 lines
5. **`lipidomics_ingestion.py`** - 630 lines
6. **`transcriptomics_ingestion.py`** - 625 lines
7. **`metadata_semantic.py`** - 605 lines

### Files Over 500 Lines (Medium Priority)

8. **`metabolomics_ingestion.py`** - 548 lines
9. **`proteomics_ingestion.py`** - 547 lines
10. **`rag_engine.py`** - 545 lines
11. **`mwtab_extraction.py`** - 508 lines
12. **`zotero_ingest.py`** - 503 lines

---

## 🎯 Refactoring Strategy

### Principle: Single Responsibility

Each module should have one clear purpose. Break up files by:
- **Functional area** (CRUD, extraction, scoring, etc.)
- **Data type** (signatures, features, datasets)
- **Operation type** (queries, ingestion, matching)

---

## 📋 Detailed Refactoring Plan

### 1. `cross_omics_reasoning.py` (892 lines) ⚠️ **HIGHEST PRIORITY**

**Current Structure**: Single file with all cross-omics functions

**Split Into**:
```
amprenta_rag/query/cross_omics/
├── __init__.py                    # Exports
├── program_summary.py             # cross_omics_program_summary
├── signature_summary.py           # cross_omics_signature_summary
├── feature_summary.py             # cross_omics_feature_summary
├── dataset_summary.py             # cross_omics_dataset_summary
├── helpers.py                     # Shared helper functions
└── prompt_templates.py            # LLM prompt templates
```

**Functions to Extract**:
- `_fetch_notion_page()` → `helpers.py`
- `_extract_relation_ids()` → `helpers.py`
- `_extract_select_values()` → `helpers.py`
- `_extract_text_property()` → `helpers.py`
- `_get_chunk_text()` → `helpers.py`
- `_retrieve_chunks_for_objects()` → `helpers.py`
- `_group_chunks_by_omics_type()` → `helpers.py`
- `_synthesize_cross_omics_summary()` → `helpers.py`

**Estimated Reduction**: 892 → ~150 lines per file (6 files)

---

### 2. `feature_extraction.py` (821 lines)

**Current Structure**: Mixed feature extraction and linking functions

**Split Into**:
```
amprenta_rag/ingestion/features/
├── __init__.py                    # Exports
├── extraction.py                  # Feature extraction functions
├── normalization.py               # Normalization functions
├── linking.py                     # Feature linking to Notion
└── constants.py                   # Metabolite synonyms, lists
```

**Functions to Organize**:
- Normalization: `normalize_metabolite_name()`, etc. → `normalization.py`
- Extraction: `extract_features_from_mwtab()`, etc. → `extraction.py`
- Linking: `link_feature()`, `_find_or_create_feature_page()` → `linking.py`
- Constants: `METABOLITE_SYNONYMS`, `AMINO_ACIDS`, etc. → `constants.py`

**Estimated Reduction**: 821 → ~200 lines per file (4 files)

---

### 3. `signature_notion_crud.py` (801 lines)

**Current Structure**: All signature CRUD operations in one file

**Split Into**:
```
amprenta_rag/ingestion/signatures/
├── __init__.py                    # Exports
├── signature_crud.py              # Signature page CRUD
├── component_crud.py              # Component page CRUD
├── species_crud.py                # Lipid species CRUD
└── short_id.py                    # Short ID generation
```

**Functions to Organize**:
- `generate_signature_short_id()` → `short_id.py`
- `find_or_create_signature_page()` → `signature_crud.py`
- `update_signature_page_if_needed()` → `signature_crud.py`
- `find_or_create_component_page()` → `component_crud.py`
- `find_or_create_lipid_species_page()` → `species_crud.py`
- `update_lipid_species_synonyms()` → `species_crud.py`

**Estimated Reduction**: 801 → ~200 lines per file (4 files)

**Note**: This partially duplicates existing `signature_ingestion.py` structure. Consider consolidation.

---

### 4. `signature_matching.py` (670 lines)

**Current Structure**: Matching, scoring, and writeback functions

**Split Into**:
```
amprenta_rag/ingestion/signatures/
├── matching.py                    # Signature matching logic (keep)
├── scoring.py                     # Scoring functions
└── writeback.py                   # Notion writeback functions
```

**Functions to Organize**:
- Scoring: `score_signature_against_dataset()`, etc. → `scoring.py`
- Writeback: `update_dataset_with_signature_matches()` → `writeback.py`
- Matching: Core matching logic → `matching.py` (simplified)

**Estimated Reduction**: 670 → ~250 lines per file (3 files)

---

### 5. `lipidomics_ingestion.py` (630 lines)

**Current Structure**: Complete ingestion pipeline

**Analysis**: This is already reasonably modular. Consider:
- Extract species normalization to shared module
- Extract file parsing to shared utilities
- Keep orchestration logic in main file

**Potential Split**:
```
amprenta_rag/ingestion/lipidomics/
├── __init__.py
├── ingestion.py                   # Main orchestration (~300 lines)
├── normalization.py               # Species normalization (~150 lines)
└── file_parsing.py                # File parsing utilities (~180 lines)
```

---

### 6. `transcriptomics_ingestion.py` (625 lines)

**Similar to lipidomics** - consider same pattern.

---

### 7. `metadata_semantic.py` (605 lines)

**Current Structure**: Semantic metadata extraction

**Split Into**:
```
amprenta_rag/ingestion/metadata/
├── __init__.py
├── extraction.py                  # Main extraction function
├── disease_extraction.py          # Disease detection
├── matrix_extraction.py           # Matrix detection
├── signature_extraction.py        # Signature detection
└── patterns.py                    # Regex patterns
```

---

## 🎯 Recommended Priority Order

### Phase 1: High-Impact Refactoring (Start Here)

1. **`cross_omics_reasoning.py`** (892 lines)
   - Highest value - split into logical modules
   - Clear separation of concerns
   - Easy to test independently

2. **`feature_extraction.py`** (821 lines)
   - Already has clear functional boundaries
   - Extracting normalization helps other modules
   - Better code reuse

3. **`signature_notion_crud.py`** (801 lines)
   - Clear separation: signature/component/species
   - Easier to maintain
   - Better testability

### Phase 2: Medium-Impact Refactoring

4. **`signature_matching.py`** (670 lines)
5. **`metadata_semantic.py`** (605 lines)
6. **Omics ingestion files** (630, 625, 548, 547 lines)

---

## ✅ Refactoring Principles

1. **Maintain Backward Compatibility**
   - Use `__init__.py` to re-export functions
   - Add deprecation warnings if needed
   - Gradual migration path

2. **Single Responsibility**
   - Each module does one thing well
   - Clear boundaries between modules

3. **Testability**
   - Smaller files = easier to test
   - Clear interfaces between modules

4. **Idempotency**
   - All operations remain idempotent
   - No breaking changes to behavior

---

## 📋 Implementation Checklist

For each file to refactor:

- [ ] Analyze current structure
- [ ] Identify logical splits
- [ ] Create new module structure
- [ ] Move functions to new modules
- [ ] Update imports
- [ ] Add `__init__.py` exports
- [ ] Test all functionality
- [ ] Update documentation

---

## 🚀 Start With: `cross_omics_reasoning.py`

**Why Start Here**:
- Largest file (892 lines)
- Clear functional boundaries
- High impact on maintainability
- Easy to split logically

**Estimated Time**: 2-3 hours

---

**Ready to proceed with refactoring? Let's start with `cross_omics_reasoning.py`!**
