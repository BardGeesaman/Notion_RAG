# Programs & Experiments Database - Verification Results

## ✅ Status Check

### Experiments Database ✅ WORKING!

Successfully connected and listing experiments:

- ✅ Found **7 experiments** in the database
- ✅ Successfully extracting experiment names and IDs
- ✅ Ready for cross-omics testing

**Example experiments found:**
1. 🧬 Multi-Omics Experiment
   - ID: `288fa0da-a8a5-49ae-a5ee-2d5b94eaf661`

2. (And 6 more experiments...)

### Programs Database ⚠️ Needs Attention

Encountered an error when querying:
- ❌ HTTP 400 Bad Request error
- Possible causes:
  1. Database ID might be incorrect
  2. Database permissions issue
  3. Database structure mismatch

**Action needed:**
- Verify the `NOTION_PROGRAMS_DB_ID` in `.env` matches the actual Programs database
- Check that the Notion integration has access to the Programs database
- Confirm the database ID is 32 characters (no dashes)

---

## ✅ What's Working

### Experiments ✅
```bash
# List experiments
python scripts/list_cross_omics_test_ids.py --type experiments --limit 10

# This works perfectly!
```

### Other Databases ✅
All other databases working:
- ✅ Signatures
- ✅ Datasets (Experimental Data Assets)
- ✅ Gene Features
- ✅ Protein Features
- ✅ Metabolite Features
- ✅ Lipid Species

---

## 🔧 Troubleshooting Programs Database

### Option 1: Verify Database ID

1. Open Notion → Navigate to Programs database
2. Click "..." menu → "Copy link"
3. Extract the 32-character ID from the URL
4. Make sure it matches exactly what's in `.env`

The ID should be:
- 32 characters long
- No dashes (if you see dashes in URL, remove them)
- Example: `2b5adf6142ab802c8100c772dd41c650`

### Option 2: Check Database Permissions

1. In Notion, go to the Programs database
2. Click "..." menu → "Connections"
3. Ensure your Notion integration is connected
4. The integration needs "Read" access at minimum

### Option 3: Test with Simple Query

Try a simpler test to verify the database ID works:

```bash
# This will show more detailed error info
python3 -c "
import sys
sys.path.insert(0, '.')
from amprenta_rag.config import get_config
from amprenta_rag.clients.notion_client import notion_headers
import requests

cfg = get_config()
db_id = cfg.notion.programs_db_id
print(f'Testing Programs DB ID: {db_id[:16]}...')
print(f'Full ID: {db_id}')

if db_id:
    url = f'{cfg.notion.base_url}/databases/{db_id}'
    resp = requests.get(url, headers=notion_headers(), timeout=30)
    print(f'Status: {resp.status_code}')
    if resp.status_code != 200:
        print(f'Error: {resp.text}')
"
```

---

## 🎯 Current Capabilities

### ✅ Can Test Right Now:

1. **Experiments** - Fully working!
   ```bash
   python scripts/list_cross_omics_test_ids.py --type experiments
   ```

2. **All Feature Types** - Working!
   ```bash
   python test_cross_omics.py --feature "gene:GFAP"
   python test_cross_omics.py --feature "metabolite:Glutamate"
   ```

3. **Signatures** - Working!
   ```bash
   python test_cross_omics.py --signature <signature_id>
   ```

4. **Datasets** - Working!
   ```bash
   python test_cross_omics.py --dataset <dataset_id>
   ```

### ⚠️ Programs - Needs Fix:

- Program listing: Error (needs DB ID verification)
- Program summaries: Will work once listing is fixed

---

## ✅ Summary

- **Experiments Database**: ✅ **FULLY OPERATIONAL**
- **Programs Database**: ⚠️ **Needs ID verification/permissions check**
- **All Other Features**: ✅ **Working perfectly**

You can proceed with testing experiments and all other features while troubleshooting the Programs database!

