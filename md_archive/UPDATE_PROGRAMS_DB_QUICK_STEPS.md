# Quick Steps: Update Programs Database ID

## ✅ New Database Created

A fresh Programs database with single data source is ready!

---

## 🎯 Quick Setup (3 Steps)

### Step 1: Share Database with Integration ⚠️ CRITICAL

1. Open: `https://www.notion.so/bde04fc4ed6640a1b392445d7e1a08ae`
2. Click "..." menu → "Connections"
3. Add your Notion integration
4. Grant "Read" access (minimum)

### Step 2: Get Correct Database ID

From the URL: `https://www.notion.so/bde04fc4ed6640a1b392445d7e1a08ae`

**Database ID**: `bde04fc4ed6640a1b392445d7e1a08ae` (from URL)

OR use the provided ID: `2bdadf6142ab80518c77e41a60f48e31`

### Step 3: Update .env File

```bash
# Edit .env file, change this line:
NOTION_PROGRAMS_DB_ID=bde04fc4ed6640a1b392445d7e1a08ae

# OR if the provided ID works:
NOTION_PROGRAMS_DB_ID=2bdadf6142ab80518c77e41a60f48e31
```

---

## ✅ Test

```bash
# Test database access
python scripts/list_cross_omics_test_ids.py --type programs --limit 3
```

**Expected**: Lists programs (not an error)

---

## 📝 Database Info

- **URL**: `https://www.notion.so/bde04fc4ed6640a1b392445d7e1a08ae`
- **Schema**: Identical to original
- **Status**: Single data source ✅
- **API-queryable**: Yes (after sharing)

---

**Ready to update? Share the database first, then update .env!**

