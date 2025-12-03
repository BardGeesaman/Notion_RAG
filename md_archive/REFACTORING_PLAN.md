# Refactoring Plan: Breaking Down Large Files

## Overview
This plan outlines how to refactor the three largest ingestion modules into smaller, more maintainable components while preserving functionality and backward compatibility.

## Current State

### File Sizes
- `dataset_ingestion.py`: **1,132 lines** (6 functions)
- `signature_ingestion.py`: **1,109 lines** (12 functions)
- `email_ingestion.py`: **676 lines** (6 functions)

### Goals
1. Reduce file sizes to <500 lines each
2. Improve separation of concerns
3. Enhance testability
4. Maintain backward compatibility
5. Improve code reusability

---

## Refactoring Strategy

### 1. `dataset_ingestion.py` → Split into 3 modules

#### **Current Structure:**
- `_extract_mwtab_from_page_content()` - 90 lines
- `_extract_metadata_from_mwtab()` - 213 lines
- `_update_experimental_data_asset_metadata()` - 120 lines
- `_fetch_notion_page()` - 22 lines
- `_update_notion_page_with_embeddings()` - 74 lines
- `ingest_dataset()` - 565 lines (main orchestration)

#### **Proposed Split:**

**A. `amprenta_rag/ingestion/mwtab_extraction.py`** (~300 lines)
- `extract_mwtab_from_page_content()` (public)
- `extract_metadata_from_mwtab()` (public)
- `extract_study_id_from_page_properties()` (helper)
- `fetch_mwtab_from_api()` (helper - MW API fallback)

**B. `amprenta_rag/ingestion/dataset_notion_utils.py`** (~200 lines)
- `update_dataset_embedding_metadata()` (public)
- `update_dataset_scientific_metadata()` (public)
- `fetch_dataset_page()` (public)

**C. `amprenta_rag/ingestion/dataset_ingestion.py`** (~450 lines)
- `ingest_dataset()` (main orchestration)
- Imports from `mwtab_extraction` and `dataset_notion_utils`

**Benefits:**
- mwTab logic isolated and testable
- Notion operations separated
- Main orchestration cleaner

---

### 2. `signature_ingestion.py` → Split into 4 modules

#### **Current Structure:**
- `_generate_short_id()` - 20 lines
- `_find_or_create_signature_page()` - 150 lines
- `_update_signature_page_if_needed()` - 95 lines
- `_find_or_create_component_page()` - 150 lines
- `_find_or_create_lipid_species_page()` - 140 lines
- `_update_lipid_species_synonyms()` - 45 lines
- `_link_component_to_lipid_species()` - 60 lines
- `ingest_signature_from_file()` - 170 lines
- `link_signature_to_source()` - 80 lines
- `link_component_to_metabolite_feature()` - 45 lines
- `embed_signature()` - 95 lines
- `_fetch_notion_page_helper()` - 7 lines

#### **Proposed Split:**

**A. `amprenta_rag/ingestion/signature_notion_crud.py`** (~400 lines)
- `find_or_create_signature_page()` (public)
- `update_signature_page_if_needed()` (public)
- `find_or_create_component_page()` (public)
- `find_or_create_lipid_species_page()` (public)
- `update_lipid_species_synonyms()` (public)
- `generate_signature_short_id()` (public)

**B. `amprenta_rag/ingestion/signature_linking.py`** (~200 lines)
- `link_component_to_lipid_species()` (public)
- `link_signature_to_source()` (public)
- `link_component_to_metabolite_feature()` (public)

**C. `amprenta_rag/ingestion/signature_embedding.py`** (~100 lines)
- `embed_signature()` (public)

**D. `amprenta_rag/ingestion/signature_ingestion.py`** (~400 lines)
- `ingest_signature_from_file()` (main orchestration)
- Imports from `signature_notion_crud`, `signature_linking`, `signature_embedding`

**Benefits:**
- CRUD operations isolated
- Linking logic separated
- Embedding logic isolated
- Main orchestration cleaner

---

### 3. `email_ingestion.py` → Split into 2 modules

#### **Current Structure:**
- `cleanup_orphaned_chunks()` - 150 lines
- `delete_email_and_chunks()` - 135 lines
- `ingest_email()` - 260 lines
- `batch_ingest_emails()` - 50 lines
- Helper functions: `_notion_base_url()`, `_rag_db_id()`

#### **Proposed Split:**

**A. `amprenta_rag/ingestion/email_cleanup.py`** (~200 lines)
- `cleanup_orphaned_chunks()` (public)
- `delete_email_and_chunks()` (public)
- Helper functions moved here

**B. `amprenta_rag/ingestion/email_ingestion.py`** (~450 lines)
- `ingest_email()` (main ingestion)
- `batch_ingest_emails()` (batch orchestration)
- Imports cleanup functions from `email_cleanup`

**Benefits:**
- Cleanup logic isolated
- Main ingestion focused
- Easier to test cleanup separately

---

## Implementation Plan

### Phase 1: Extract Utilities (Low Risk)
1. Create `mwtab_extraction.py` from `dataset_ingestion.py`
2. Create `signature_notion_crud.py` from `signature_ingestion.py`
3. Create `email_cleanup.py` from `email_ingestion.py`
4. Update imports in original files
5. Run tests to verify

### Phase 2: Extract Linking/Embedding (Medium Risk)
1. Create `signature_linking.py`
2. Create `signature_embedding.py`
3. Create `dataset_notion_utils.py`
4. Update imports
5. Run tests

### Phase 3: Refactor Main Orchestration (Low Risk)
1. Clean up main orchestration functions
2. Remove duplicate helpers
3. Update all imports across codebase
4. Run full test suite

### Phase 4: Documentation & Cleanup
1. Add module docstrings
2. Update README if needed
3. Verify all scripts still work
4. Commit changes

---

## Backward Compatibility Strategy

### Import Aliases
Keep original function names available via `__all__` exports:

```python
# In new modules
__all__ = [
    "extract_mwtab_from_page_content",
    "extract_metadata_from_mwtab",
    # ...
]

# In original files (temporary)
from amprenta_rag.ingestion.mwtab_extraction import (
    extract_mwtab_from_page_content as _extract_mwtab_from_page_content,
    extract_metadata_from_mwtab as _extract_metadata_from_mwtab,
)
```

### Deprecation Warnings (Optional)
Add deprecation warnings to old function locations pointing to new modules.

---

## Testing Strategy

### Unit Tests
- Test each extracted module independently
- Mock Notion/Pinecone API calls
- Test error handling paths

### Integration Tests
- Test full ingestion pipelines
- Verify Notion writebacks
- Verify Pinecone upserts

### Regression Tests
- Run existing test suite after each phase
- Verify all scripts still work
- Check for import errors

---

## File Size Targets

### After Refactoring:
- `dataset_ingestion.py`: ~450 lines ✅
- `signature_ingestion.py`: ~400 lines ✅
- `email_ingestion.py`: ~450 lines ✅
- `mwtab_extraction.py`: ~300 lines ✅
- `signature_notion_crud.py`: ~400 lines ✅
- `signature_linking.py`: ~200 lines ✅
- `signature_embedding.py`: ~100 lines ✅
- `dataset_notion_utils.py`: ~200 lines ✅
- `email_cleanup.py`: ~200 lines ✅

**Total:** ~2,700 lines (vs ~2,900 lines before)
**Benefit:** Better organization, easier maintenance, improved testability

---

## Risk Assessment

### Low Risk:
- Extracting utility functions
- Creating new modules
- Adding import statements

### Medium Risk:
- Changing function signatures
- Moving shared helpers
- Updating all import statements

### Mitigation:
- Implement in phases
- Keep original files working during transition
- Comprehensive testing after each phase
- Git commits after each successful phase

---

## Next Steps

1. **Review this plan** - Confirm approach and priorities
2. **Start with Phase 1** - Extract utilities (safest first step)
3. **Test thoroughly** - Ensure nothing breaks
4. **Iterate** - Continue with subsequent phases

---

## Questions to Consider

1. Should we maintain backward compatibility indefinitely or deprecate old imports?
2. Do we want to add type hints to all extracted functions?
3. Should we create a `__init__.py` that re-exports everything for convenience?
4. Are there any other large files (>500 lines) that should be refactored?
