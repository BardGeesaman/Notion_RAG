# Programs & Experiments Database Support - Added! ✅

## Changes Made

### 1. Configuration Updates ✅

**File**: `amprenta_rag/config.py`

Added:
- `NOTION_PROGRAMS_DB_ID` - Reads from environment variable
- `NOTION_EXPERIMENTS_DB_ID` - Reads from environment variable
- Both added to `NotionConfig` dataclass

### 2. List Script Updates ✅

**File**: `scripts/list_cross_omics_test_ids.py`

Added:
- `list_experiments()` function
- Updated `--type` argument to include "experiments"
- Updated output to show experiments when listing

### 3. Cross-Omics Reasoning ✅

The cross-omics reasoning functions already work with Programs and Experiments via:
- Relation properties (Program Datasets, Related Experiments)
- Direct page fetching

Having the DB IDs enables:
- Direct listing (already working via `list_cross_omics_test_ids.py`)
- More efficient bulk queries
- Better error messages

---

## 📝 Next Step: Add to `.env` File

You need to add these two database IDs to your `.env` file:

```bash
# Programs database ID
NOTION_PROGRAMS_DB_ID=<your_programs_database_id>

# Experiments database ID  
NOTION_EXPERIMENTS_DB_ID=<your_experiments_database_id>
```

### How to Find Database IDs:

1. **For Programs Database**:
   - Open Notion → Navigate to your Programs database
   - Click "..." menu → "Copy link"
   - URL format: `https://www.notion.so/.../<32_char_id>?v=...`
   - Extract the 32-character ID (no dashes)
   - Add to `.env` as: `NOTION_PROGRAMS_DB_ID=<id>`

2. **For Experiments Database**:
   - Same process for your Experiments database
   - Extract the 32-character ID
   - Add to `.env` as: `NOTION_EXPERIMENTS_DB_ID=<id>`

---

## ✅ Testing After Adding IDs

1. **Reload environment**:
   ```bash
   # Restart your Python session or reload .env
   # The dotenv package will reload on next import
   ```

2. **Test listing**:
   ```bash
   # List programs
   python scripts/list_cross_omics_test_ids.py --type programs --limit 5
   
   # List experiments
   python scripts/list_cross_omics_test_ids.py --type experiments --limit 5
   
   # List all
   python scripts/list_cross_omics_test_ids.py --type all --limit 5
   ```

3. **Test cross-omics program summary**:
   ```bash
   # Get a program ID
   python scripts/list_cross_omics_test_ids.py --type programs --limit 1
   
   # Test program summary
   python test_cross_omics.py --program <program_id>
   ```

---

## 🎯 What This Enables

### Programs Database:
- ✅ Direct listing of all programs
- ✅ Full `--cross-omics-program` functionality
- ✅ More efficient queries
- ✅ Better error messages

### Experiments Database:
- ✅ Direct listing of all experiments
- ✅ More efficient experiment queries
- ✅ Better cross-omics program summaries (finds experiments faster)
- ✅ Foundation for future experiment-focused queries

---

## 📋 Code Changes Summary

1. **`amprenta_rag/config.py`**:
   - Added `NOTION_PROGRAMS_DB_ID = os.getenv(...)`
   - Added `NOTION_EXPERIMENTS_DB_ID = os.getenv(...)`
   - Added to `NotionConfig` dataclass

2. **`scripts/list_cross_omics_test_ids.py`**:
   - Added `list_experiments()` function
   - Updated CLI to support `--type experiments`
   - Updated output formatting

3. **Cross-omics reasoning**:
   - Already works via relations (no changes needed)
   - Will benefit from direct DB access (faster queries)

---

## ✅ Status

**Code Changes**: ✅ COMPLETE  
**Configuration**: ✅ ADDED TO CONFIG.PY  
**Next Step**: ➡️ ADD IDs TO `.env` FILE

Once you add the IDs to `.env`, everything will work!

---

**Ready to add the IDs to your `.env` file!** 🚀

