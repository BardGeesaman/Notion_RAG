# Update Programs Database ID - Setup Guide

## ✅ New Programs Database Created

A fresh **Pipeline Programs** database with a single data source has been created. It's now API-queryable!

---

## 📋 Required Updates

### Step 1: Share Database with Integration

**CRITICAL**: The new database must be shared with your Notion integration!

1. Open the new Programs database in Notion
   - URL: `https://www.notion.so/bde04fc4ed6640a1b392445d7e1a08ae`

2. Click the "..." menu (top right)
3. Select "Connections" or "Add connections"
4. Add your Notion integration
5. Ensure it has at least "Read" access

### Step 2: Verify Database ID

From the URL: `https://www.notion.so/bde04fc4ed6640a1b392445d7e1a08ae`

The database ID is: **`bde04fc4ed6640a1b392445d7e1a08ae`** (32 characters, no dashes)

**Note**: The ID provided (`2bdadf6142ab80518c77e41a60f48e31`) might be different. We should verify which one is correct.

### Step 3: Update .env File

Update your `.env` file with the correct database ID:

```bash
# OLD (multi-source, not queryable)
NOTION_PROGRAMS_DB_ID=2b5adf6142ab802c8100c772dd41c650

# NEW (single source, API-queryable)
# Extract from URL: bde04fc4ed6640a1b392445d7e1a08ae
NOTION_PROGRAMS_DB_ID=bde04fc4ed6640a1b392445d7e1a08ae
```

Or if the provided ID is correct:
```bash
NOTION_PROGRAMS_DB_ID=2bdadf6142ab80518c77e41a60f48e31
```

---

## 🔍 How to Extract Database ID from Notion URL

1. Open the database in Notion
2. Look at the URL: `https://www.notion.so/<database_id>?v=...`
3. The `<database_id>` is 32 characters (no dashes)
4. Copy that ID to `.env`

**Example**:
- URL: `https://www.notion.so/bde04fc4ed6640a1b392445d7e1a08ae?pvs=21`
- Database ID: `bde04fc4ed6640a1b392445d7e1a08ae`

---

## ✅ Verification Steps

### Step 1: Share Database

First, ensure the database is shared with your integration (see Step 1 above).

### Step 2: Test Database Access

Run this test after updating `.env`:

```bash
python3 -c "
import sys
sys.path.insert(0, '.')
from amprenta_rag.config import get_config
from amprenta_rag.clients.notion_client import notion_headers
import requests

# Try both possible IDs
test_ids = [
    'bde04fc4ed6640a1b392445d7e1a08ae',  # From URL
    '2bdadf6142ab80518c77e41a60f48e31',  # Provided ID
]

cfg = get_config()

for db_id in test_ids:
    print(f'Testing ID: {db_id[:16]}...')
    try:
        url = f'{cfg.notion.base_url}/databases/{db_id}/query'
        payload = {'page_size': 1}
        resp = requests.post(url, headers=notion_headers(), json=payload, timeout=30)
        
        if resp.status_code == 200:
            print(f'✅ SUCCESS with ID: {db_id}')
            print(f'   Use this ID in .env: NOTION_PROGRAMS_DB_ID={db_id}')
            break
        elif resp.status_code == 404:
            print(f'   ❌ 404 - Database not found or not shared with integration')
        else:
            print(f'   ❌ Error {resp.status_code}: {resp.text[:100]}')
    except Exception as e:
        print(f'   ❌ Error: {e}')
"
```

### Step 3: Update .env with Working ID

Once you identify the correct ID, update `.env`:

```bash
NOTION_PROGRAMS_DB_ID=<working_id_from_test>
```

### Step 4: Test Listing Programs

```bash
python scripts/list_cross_omics_test_ids.py --type programs --limit 5
```

Should show:
- ✅ Programs listed with names and IDs
- ✅ No errors

---

## 🎯 Expected Results

After completing all steps:

1. ✅ Database is shared with integration
2. ✅ `.env` has correct database ID
3. ✅ `list_cross_omics_test_ids.py --type programs` works
4. ✅ Programs can be listed and queried via API
5. ✅ Cross-omics program summaries will work

---

## 📝 Migration Notes

The new database is fresh (empty). You mentioned:

- **AMP-001** and **AMP-002** need to be migrated from old database
- Old database: `2b5adf6142ab802c8100c772dd41c650`
- New database: (will be determined after verification)

After migration, the cross-omics system will work perfectly!

---

**Priority Actions**:
1. ✅ Share new database with integration
2. ✅ Verify correct database ID
3. ✅ Update `.env` file
4. ✅ Test API access
5. ✅ Migrate AMP-001 and AMP-002

