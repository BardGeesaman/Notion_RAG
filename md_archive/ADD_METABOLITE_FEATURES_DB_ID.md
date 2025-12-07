# Add Metabolite Features Database ID to .env

## What to Add

Add this line to your `.env` file:

```bash
NOTION_METABOLITE_FEATURES_DB_ID=<your-database-id>
```

## How to Get the Database ID

1. Open the **Metabolite Features** database in Notion
2. Click the "..." menu → "Copy link"
3. The URL will look like:
   ```
   https://www.notion.so/<workspace>/<database-id>?v=...
   ```
4. Copy the database ID (32 hex characters, no dashes)
5. Add it to your `.env` file

## Example

After adding, your `.env` should have:

```bash
NOTION_EXP_DATA_DB_ID=2b6adf6142ab80e9820cee8778c95f91
NOTION_GENE_FEATURES_DB_ID=7e9e2eff9a57454b9798d243d31d2474
NOTION_LIPID_SPECIES_DB_ID=22fcb28946854dfdb6d0a89a3a665e12
NOTION_PROTEIN_FEATURES_DB_ID=57dd7eb883e644dba8d202071d1265af
NOTION_METABOLITE_FEATURES_DB_ID=<your-metabolite-features-id>  # ADD THIS
NOTION_SIGNATURE_COMPONENT_DB_ID=ba5657beee034bc69a234a7b9e6667b5
NOTION_SIGNATURE_DB_ID=18d9e6a95b64463994ceee9078a5250b
```

## After Adding

Reload your environment:

```bash
# If using source to reload
source .env

# Or restart your shell/environment
```

## Benefits

Adding this will enable:
- ✅ Full feature linking for metabolomics datasets
- ✅ Automatic metabolite feature page creation
- ✅ Dataset → metabolite feature relations
- ✅ Complete feature-level knowledge graph for metabolomics

## Current Status

**Without this DB ID:**
- Metabolomics ingestion works ✅
- Feature linking is skipped (gracefully) ⚠️

**With this DB ID:**
- Metabolomics ingestion works ✅
- Feature linking works ✅✅✅

