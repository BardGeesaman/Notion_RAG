# Adding Programs & Experiments Database Support

## Changes Made

### 1. Configuration Updates ✅

**File**: `amprenta_rag/config.py`

Added support for:
- `NOTION_PROGRAMS_DB_ID` - Programs database
- `NOTION_EXPERIMENTS_DB_ID` - Experiments database

These are now:
- Read from environment variables (`.env` file)
- Added to `NotionConfig` dataclass
- Available throughout the codebase

### 2. List Script Updates ✅

**File**: `scripts/list_cross_omics_test_ids.py`

- Updated `list_programs()` to use `cfg.notion.programs_db_id`
- Clearer error message showing required env var

### 3. Cross-Omics Reasoning ✅

The cross-omics reasoning functions already support:
- Finding experiments via relations (doesn't need Experiments DB)
- Finding programs via relations (doesn't need Programs DB)

However, having the DB IDs enables:
- Direct listing of all programs/experiments
- More efficient queries
- Better error messages

---

## 📝 Required: Add to `.env` File

Add these two lines to your `.env` file:

```bash
# Programs database ID (get from Notion)
NOTION_PROGRAMS_DB_ID=<your_programs_database_id>

# Experiments database ID (get from Notion)
NOTION_EXPERIMENTS_DB_ID=<your_experiments_database_id>
```

### How to Find Database IDs:

1. Open your Notion workspace
2. Navigate to the Programs database
3. Click "..." menu → "Copy link"
4. The URL will look like: `https://www.notion.so/YOUR_WORKSPACE/<database_id>?v=...`
5. The `<database_id>` is a 32-character string (no dashes)
6. Copy that ID and add to `.env`

Repeat for Experiments database.

---

## ✅ Testing

After adding the IDs to `.env`:

1. **Reload environment**:
   ```bash
   # Restart your Python session or reload .env
   source .env  # Or just restart your terminal
   ```

2. **Test listing**:
   ```bash
   python scripts/list_cross_omics_test_ids.py --type programs --limit 5
   python scripts/list_cross_omics_test_ids.py --type all --limit 5
   ```

3. **Test cross-omics program summary**:
   ```bash
   # Get a program ID first
   python scripts/list_cross_omics_test_ids.py --type programs --limit 1
   
   # Then test
   python test_cross_omics.py --program <program_id_from_above>
   ```

---

## 🎯 Benefits

Once configured:

- ✅ Can list all programs directly
- ✅ Can list all experiments directly  
- ✅ Can test `--cross-omics-program` fully
- ✅ More efficient queries
- ✅ Better error messages

---

## ⚠️ Note

The cross-omics functions will still work without these DB IDs by using:
- Relation properties to find linked programs/experiments
- Indirect queries via datasets/signatures

But having direct DB access is more efficient and enables full functionality!

