# Notion Agent Prompt: Fix Programs Database Multiple Data Sources

## 🎯 Task

Fix the **Pipeline / Programs** database to use a single data source by removing sync/connection relationships, so it can be queried via the Notion API.

---

## 📋 Instructions for Notion Agent

**Database to Fix**: Pipeline / Programs database  
**Database ID**: `2b5adf6142ab802c8100c772dd41c650`

### Problem

The database currently has multiple data sources (connected/synced databases), causing this API error:
```
"Databases with multiple data sources are not supported in this API version."
```

### Solution Steps

1. **Open the Pipeline / Programs database** in Notion

2. **Identify sync connections**:
   - Look for sync icons (🔄) or "Connected to" indicators in the database header
   - Check if the database is synced/connected to other databases
   - Note which databases are connected

3. **Remove sync connections**:
   - Click on any sync/connection indicator
   - Select "Disconnect" or "Remove connection"
   - Confirm the disconnection
   - Repeat for all connected data sources until only one source remains

4. **Verify**:
   - All program entries should still be visible
   - Database properties (Program, Target, Indication, Stage, etc.) unchanged
   - Relations (Experiments, Tasks, Questions, etc.) still intact
   - No sync indicators visible

### Important Notes

- **Data is preserved**: Removing sync connections does NOT delete any programs or data
- **Relations remain**: All relation properties (Experiments, Tasks, etc.) continue working
- **Structure unchanged**: Database schema and properties remain the same

### Expected Outcome

After removing sync connections:
- Database uses a single data source
- All programs remain visible
- API queries will work (status 200)
- Can list programs via `list_cross_omics_test_ids.py`

---

## ✅ Verification

After making changes, verify success with:

```bash
python scripts/list_cross_omics_test_ids.py --type programs --limit 3
```

**Success indicators**:
- ✅ Returns program names and IDs (not an error)
- ✅ Shows "Programs (showing up to 3):" with results
- ✅ No "400 Bad Request" errors

---

## 📝 Summary

**Action**: Remove sync/connection relationships from Programs database  
**Time**: 2-5 minutes  
**Risk**: Low (data preserved, only connections removed)  
**Impact**: Enables full API access for cross-omics reasoning

---

**Ready to proceed? Remove sync connections and verify with the test script!**

