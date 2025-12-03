# Programs Database - Multiple Data Sources Issue

## 🔍 Root Cause Identified

The Programs database has **multiple data sources** (connected/synced databases), which the current Notion API version doesn't support for direct queries.

**Error Message:**
```
Databases with multiple data sources are not supported in this API version.
```

**Details:**
- Database ID: `2b5adf6142ab802c8100c772dd41c650`
- Has 2 child data sources (synced/connected databases)
- Current API version doesn't support querying these

---

## ✅ Solutions

### Option 1: Query via Relations (Recommended)

Instead of querying the Programs database directly, we can find programs through their relations from other databases:

- **From Experiments**: Programs are linked via "Experiments" relation
- **From Datasets**: Programs might be linked via "Related Programs"
- **From Signatures**: If signatures link to programs

This approach works around the API limitation!

### Option 2: Query Individual Data Sources

If we have the child data source IDs, we could query them individually:
- Child source 1: `2b5adf61-42ab-8061-84c8-000b48092bcd`
- Child source 2: `2b5adf61-42ab-8018-8d60-000b98bb9297`

### Option 3: Update Notion API Version (Future)

The error mentions minimum API version `2025-09-03`, but this might require:
- Updating the Notion API client
- Using a newer API endpoint format

---

## 🔧 Implementation

I've already updated the code to:
1. Use the correct property names ("Program" instead of "Name")
2. Use the correct relation names ("Experiments" instead of "Related Experiments")

The cross-omics program summary function will still work by:
- Fetching the program page directly (this works)
- Finding linked experiments via the "Experiments" relation
- Finding datasets via experiments (datasets link to experiments)

---

## ✅ What Still Works

### Direct Program Page Access ✅
- You can still fetch individual program pages by ID
- Cross-omics program summaries will work if you provide a program page ID

### Via Relations ✅
- Find programs through experiments
- Find programs through datasets
- Cross-omics reasoning via linked entities

### Program Summary Function ✅
- `cross_omics_program_summary(program_page_id)` will work
- It fetches the page directly (not via database query)
- Finds linked experiments and datasets via relations

---

## 📝 Updated Code

The code has been updated to match your actual schema:

1. **Property names**:
   - Uses "Program" as title property (with "Name" fallback)
   - Uses "Experiments" relation (with "Related Experiments" fallback)

2. **Error handling**:
   - Gracefully handles database query failures
   - Falls back to relation-based discovery

---

## 🎯 Testing

You can still test program summaries if you have a program page ID:

```bash
# If you have a program page ID from elsewhere
python test_cross_omics.py --program <program_page_id>

# The function will:
# 1. Fetch the program page directly (works!)
# 2. Find linked experiments via "Experiments" relation
# 3. Find datasets via experiments
# 4. Generate cross-omics summary
```

---

## ✅ Summary

**Status**: Programs database direct querying is blocked by Notion API limitation (multiple data sources)

**Workaround**: 
- Program summaries work via direct page access + relations
- Can find programs through their relations from other databases

**Impact**: 
- Can't list all programs directly
- Can still use program summaries if you have the page ID
- Cross-omics reasoning via relations still works

This is a Notion API limitation, not a code issue! ✅

