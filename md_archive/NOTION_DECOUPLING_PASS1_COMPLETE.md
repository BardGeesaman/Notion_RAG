# Notion Decoupling Pass 1 - Complete

**Date**: 2025-12-05  
**Status**: ✅ **COMPLETE**

## Summary

All Notion-dependent modules have been properly isolated with runtime guards. Core library code and scripts now run cleanly with `ENABLE_NOTION_SYNC=false`.

---

## Changes Made

### 1. `amprenta_rag/ingestion/__init__.py`
- ✅ **Removed** Notion-dependent exports:
  - `dataset_ingestion.ingest_dataset` (legacy Notion-based)
  - `experiments_ingestion.ingest_experiment` (legacy Notion-based)
  - `zotero_collection.incremental_ingest_collection` (legacy Notion-based)
- ✅ **Kept** Postgres-first exports:
  - `email_ingestion.*` (Postgres-first)
  - `zotero_ingest.ingest_zotero_item` (Postgres-first)
- ✅ Added documentation explaining legacy modules are not exported

### 2. Runtime Guards Added

All legacy Notion-dependent modules now have top-of-file runtime guards:

#### ✅ `amprenta_rag/ingestion/dataset_ingestion.py`
```python
cfg = get_config()
if not getattr(cfg.pipeline, "enable_notion_sync", False):
    raise RuntimeError("Notion sync is disabled; this module is legacy/migration only and must not be imported in core flows.")
```

#### ✅ `amprenta_rag/ingestion/experiments_ingestion.py`
```python
cfg = get_config()
if not getattr(cfg.pipeline, "enable_notion_sync", False):
    raise RuntimeError("Notion sync is disabled; this module is legacy/migration only and must not be imported in core flows.")
```

#### ✅ `amprenta_rag/ingestion/zotero_collection.py`
```python
cfg = get_config()
if not getattr(cfg.pipeline, "enable_notion_sync", False):
    raise RuntimeError("Notion sync is disabled; this module is legacy/migration only and must not be imported in core flows.")
```

#### ✅ `amprenta_rag/ingestion/dataset_notion_utils.py`
- Already had guard (verified)

#### ✅ `amprenta_rag/ingestion/signature_matching/signature_loader.py`
- Already had guard (verified)

#### ✅ `amprenta_rag/ingestion/signature_matching/writeback.py`
- Already had guard (verified)

#### ✅ `amprenta_rag/chemistry/notion_integration.py`
- Already had guard (verified)

### 3. Conditional Imports

#### ✅ `amprenta_rag/ingestion/signature_matching/matching.py`
- Made `fetch_all_signatures_from_notion` and `load_signature_from_notion_page` conditionally imported
- Added runtime checks to handle `None` when Notion is disabled

#### ✅ `amprenta_rag/ingestion/signature_matching/__init__.py`
- Made `update_dataset_with_signature_matches` conditionally imported
- Handles `RuntimeError` when Notion sync is disabled

---

## Verification Tests

### ✅ Core Imports (No Notion Required)
```bash
# All pass with ENABLE_NOTION_SYNC=false
✅ from amprenta_rag.ingestion import ingest_email, batch_ingest_emails, ingest_zotero_item
✅ from amprenta_rag.ingestion.postgres_ingestion.feature_linking import extract_and_link_features
✅ from amprenta_rag.ingestion.signature_matching import map_raw_lipid_to_canonical_species
✅ from amprenta_rag.chemistry import create_or_update_compound_in_postgres
```

### ✅ Legacy Modules (Raise Errors When Notion Disabled)
```bash
# All correctly raise RuntimeError with ENABLE_NOTION_SYNC=false
✅ from amprenta_rag.ingestion.dataset_ingestion import ingest_dataset
✅ from amprenta_rag.ingestion.experiments_ingestion import ingest_experiment
✅ from amprenta_rag.ingestion.zotero_collection import incremental_ingest_collection
```

### ✅ Production Scripts (Run Without Notion)
```bash
# All run successfully with ENABLE_NOTION_SYNC=false
✅ python scripts/backfill_dataset_features.py
✅ python scripts/discover_signatures_from_postgres.py
✅ python scripts/validate_e2e.py
```

---

## Module Status

| Module | Status | Notes |
|--------|--------|-------|
| `ingestion/__init__.py` | ✅ Clean | No Notion exports |
| `ingestion/dataset_ingestion.py` | ✅ Guarded | Legacy, raises error if Notion disabled |
| `ingestion/experiments_ingestion.py` | ✅ Guarded | Legacy, raises error if Notion disabled |
| `ingestion/zotero_collection.py` | ✅ Guarded | Legacy, raises error if Notion disabled |
| `ingestion/dataset_notion_utils.py` | ✅ Guarded | Legacy, raises error if Notion disabled |
| `ingestion/signature_matching/signature_loader.py` | ✅ Guarded | Legacy, raises error if Notion disabled |
| `ingestion/signature_matching/writeback.py` | ✅ Guarded | Legacy, raises error if Notion disabled |
| `ingestion/signature_matching/matching.py` | ✅ Conditional | Handles Notion disabled gracefully |
| `ingestion/signature_matching/__init__.py` | ✅ Conditional | Handles Notion disabled gracefully |
| `chemistry/notion_integration.py` | ✅ Guarded | Legacy, raises error if Notion disabled |
| `chemistry/__init__.py` | ✅ Clean | No Notion exports |

---

## Postgres-First Alternatives

For all legacy Notion-dependent modules, Postgres-first alternatives exist:

| Legacy Module | Postgres Alternative |
|---------------|---------------------|
| `dataset_ingestion.ingest_dataset` | `postgres_dataset_ingestion.ingest_dataset_from_postgres` |
| `experiments_ingestion.ingest_experiment` | `postgres_experiment_ingestion.*` |
| `zotero_collection.*` | `zotero/ingestion.py` (Postgres-first) |
| `signature_matching.signature_loader.*` | `postgres_signature_matching.*` |

---

## Next Steps

1. ✅ **Complete**: All Notion-dependent modules isolated
2. ✅ **Complete**: Core imports work without Notion
3. ✅ **Complete**: Production scripts run without Notion
4. ⏸️ **Future**: Consider deprecation warnings for legacy modules
5. ⏸️ **Future**: Update documentation to emphasize Postgres-first paths

---

## Conclusion

**Notion Decoupling Pass 1 is complete.** All core library code and scripts run cleanly with `ENABLE_NOTION_SYNC=false`. Legacy Notion-dependent modules are properly isolated and will raise clear errors if accidentally imported when Notion sync is disabled.

The codebase is now ready for production runs with Postgres as the single source of truth, with Notion sync as an optional legacy/migration feature.
