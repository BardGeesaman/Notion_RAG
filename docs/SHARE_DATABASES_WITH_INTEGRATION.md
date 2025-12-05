# Sharing Notion Databases with Integration

The chemistry and HTS databases need to be shared with your Notion integration to enable programmatic access.

## Quick Steps

1. **Open each database in Notion**
   - 🧪 Compound Features
   - 🧪 HTS Campaigns
   - 🧪 Biochemical Hits
   - 🧬 Pathways

2. **For each database:**
   - Click the "..." menu (top right)
   - Select "Share" or "Add connections"
   - Find your integration (the one with your `NOTION_API_KEY`)
   - Click "Invite" or "Add"

3. **Verify permissions:**
   - Ensure the integration has "Update" capability
   - Some databases may need "Full access" for relations

## Alternative: Share via URL

1. Open the database
2. Click "Share" in the top right
3. Click "Add people, emails, groups, or integrations"
4. Search for your integration name
5. Add it with "Can edit" permissions

## Verification

After sharing, run:

```bash
python scripts/verify_notion_setup.py
```

All databases should show as "✅ Accessible".

Then test compound creation:

```bash
python scripts/test_chemistry_notion.py
```

This should successfully create a test compound page.

## Troubleshooting

**404 Error:**
- Database not shared with integration
- Integration doesn't have access
- Database ID is incorrect

**403 Error:**
- Integration lacks permissions
- Need "Update" or "Full access"

**Database ID Format:**
- Should be 36 characters with dashes: `xxxxxxxx-xxxx-xxxx-xxxx-xxxxxxxxxxxx`
- Can be found in the database URL: `notion.so/.../db=xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx`

