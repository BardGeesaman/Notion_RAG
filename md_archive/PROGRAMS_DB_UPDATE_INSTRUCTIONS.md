# Update Programs Database ID - Instructions

## ✅ New Database Verified!

The new Programs database is ready and verified to work with the API.

---

## 📝 Update Required in .env

**Change this line in your `.env` file:**

```bash
# OLD (multi-source, not queryable)
NOTION_PROGRAMS_DB_ID=2b5adf6142ab802c8100c772dd41c650

# NEW (single source, API-queryable) ✅
NOTION_PROGRAMS_DB_ID=bde04fc4ed6640a1b392445d7e1a08ae
```

**Database Details:**
- **Name**: Pipeline Programs
- **URL**: `https://www.notion.so/bde04fc4ed6640a1b392445d7e1a08ae`
- **Status**: ✅ Single data source
- **API-queryable**: ✅ Yes

---

## ✅ Verification

The database has been tested and confirmed:
- ✅ Database exists and is accessible
- ✅ Query successful (status 200)
- ✅ Ready for API queries

---

## 🚀 After Updating .env

Once you update the `.env` file, test with:

```bash
# List programs
python scripts/list_cross_omics_test_ids.py --type programs --limit 5

# Test program summary (once you have a program page ID)
python test_cross_omics.py --program <program_page_id>
```

---

## 📋 Migration Notes

**Important**: The new database is currently empty. You'll need to:

1. **Migrate programs** (AMP-001, AMP-002) from old database to new
2. **Update .env** with new database ID (see above)
3. **Test** that everything works
4. **Archive** the old database once migration is complete

---

## ✅ Summary

- ✅ New database created with single data source
- ✅ Database ID verified and working
- ✅ Schema matches original
- ⏳ **Next step**: Update `.env` file with new ID

**Ready to update! Just change the database ID in your `.env` file!**

