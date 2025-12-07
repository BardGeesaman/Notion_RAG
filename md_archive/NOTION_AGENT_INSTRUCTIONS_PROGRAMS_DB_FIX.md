# Notion Agent Instructions: Fix Programs Database Multiple Data Sources Issue

## 🎯 Objective

Fix the **Pipeline / Programs** database to use a single data source instead of multiple data sources (connected/synced databases), so that it can be queried via the Notion API.

---

## 🔍 Problem Statement

The Programs database currently has **multiple data sources** (connected/synced databases), which prevents direct API queries. The Notion API error indicates:

```
"Databases with multiple data sources are not supported in this API version."
```

**Database ID**: `2b5adf6142ab802c8100c772dd41c650`

**Child data sources identified**:
- `2b5adf61-42ab-8061-84c8-000b48092bcd`
- `2b5adf61-42ab-8018-8d60-000b98bb9297`

---

## ✅ Solution: Consolidate to Single Data Source

### Option A: Remove Sync Connections (Recommended)

If the Programs database is synced/connected to other databases, remove those connections:

1. **Open the Programs database** in Notion
   - Navigate to: Pipeline / Programs database

2. **Check for sync connections**:
   - Look for a sync icon or "Connected to" indicator at the top of the database
   - Check database properties for any that indicate syncing

3. **Remove sync connections**:
   - Click the sync/connection indicator
   - Select "Disconnect" or "Remove connection"
   - Confirm the disconnection

4. **Verify single data source**:
   - The database should now be standalone
   - All data should remain intact
   - Database should now support API queries

### Option B: Consolidate Synced Databases

If you need to keep data from multiple sources:

1. **Export data from child data sources** (if needed)
2. **Import/merge into main Programs database**
3. **Remove sync connections**
4. **Verify all programs are in single database**

### Option C: Use Primary Database Only

If the synced databases are redundant:

1. **Identify the primary Programs database**
2. **Ensure all programs are in the primary database**
3. **Remove all sync connections**
4. **Delete or archive the synced database copies**

---

## 📋 Step-by-Step Instructions

### Step 1: Identify the Issue

1. Open Notion and navigate to the **Pipeline / Programs** database
2. Check the database header for:
   - Sync icons
   - "Connected to" indicators
   - Multiple database references

### Step 2: Check Database Properties

Look at the database properties for any that indicate syncing or connections:

- Properties that show sync status
- Relations to databases that might be synced
- Any "Source" or "Linked from" indicators

### Step 3: Remove Sync Connections

**If you see a sync/connection indicator:**

1. Click on the sync icon or "Connected to" text
2. A menu should appear with connection options
3. Select "Disconnect" or "Remove connection"
4. Confirm the action
5. Wait for Notion to process the change

**If sync is via a property:**

1. Identify which property creates the sync
2. Remove or modify that property
3. Or convert it to a regular relation property

### Step 4: Verify the Fix

After removing sync connections:

1. The database should no longer show sync indicators
2. All programs should still be visible
3. Database structure should remain the same
4. API queries should now work

---

## 🔧 Technical Verification

After making changes, verify with this test:

```bash
python3 -c "
import sys
sys.path.insert(0, '.')
from amprenta_rag.config import get_config
from amprenta_rag.clients.notion_client import notion_headers
import requests

cfg = get_config()
db_id = cfg.notion.programs_db_id
url = f'{cfg.notion.base_url}/databases/{db_id}/query'
payload = {'page_size': 1}

resp = requests.post(url, headers=notion_headers(), json=payload, timeout=30)
if resp.status_code == 200:
    print('✅ SUCCESS: Programs database can now be queried!')
else:
    print(f'❌ Still has issues: {resp.status_code}')
    print(resp.text[:200])
"
```

**Expected result**: Status code 200 (success)

---

## ⚠️ Important Notes

1. **Data Preservation**:
   - Removing sync connections does NOT delete data
   - All programs remain in the database
   - Only the connection/sync is removed

2. **Backup Recommendation**:
   - Consider exporting the Programs database before making changes
   - This ensures you have a backup if needed

3. **Relation Properties**:
   - Relations to Experiments, Tasks, Questions, etc. remain intact
   - Only sync connections are removed, not relations

4. **Time Required**:
   - Removing sync connections is quick (seconds)
   - Notion will process the change immediately

---

## ✅ Success Criteria

After completing the fix:

- ✅ Database no longer shows sync/connection indicators
- ✅ All programs remain visible in the database
- ✅ Database properties unchanged
- ✅ Relations to Experiments, Tasks, etc. still work
- ✅ API queries return status 200 (success)
- ✅ `list_cross_omics_test_ids.py --type programs` works

---

## 📝 Alternative: Manual Verification Steps

If automated verification isn't possible:

1. Open Programs database in Notion
2. Verify no sync indicators visible
3. Check that all programs are still there
4. Verify relations (Experiments, Tasks) still work
5. Run the list script to test API access

---

## 🎯 Expected Outcome

Once fixed, you'll be able to:

- ✅ List all programs via API
- ✅ Query Programs database directly
- ✅ Use `--cross-omics-program` functionality fully
- ✅ Automatically discover programs for cross-omics summaries

---

**Priority**: Medium (workarounds exist, but direct querying is preferred)

**Impact**: Enables full Programs database API access for cross-omics reasoning

**Time Estimate**: 2-5 minutes to remove sync connections

---

**Ready to proceed? Follow the steps above and verify with the test script!**

