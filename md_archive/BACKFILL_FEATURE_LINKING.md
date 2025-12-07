# Backfill Feature Linking for Existing Datasets

## Overview

Re-ingest existing lipidomics datasets to backfill feature linking. This will:
- ✅ Extract lipid species from the original file
- ✅ Create/find lipid species pages in Notion
- ✅ Link dataset to all lipid species (backfill feature linking)
- ✅ Re-score against signatures and update matches
- ✅ Update embeddings in Pinecone

## Command Template

```bash
python scripts/ingest_lipidomics.py \
  --file <original_file_path_or_re-export> \
  --dataset-page-id <existing_dataset_page_id>
```

## Example: Re-Ingest Test Dataset

```bash
python scripts/ingest_lipidomics.py \
  --file test_data/test_lipidomics_features.csv \
  --dataset-page-id 2beadf61-42ab-81fc-9e3c-f55fb0ccd975
```

## How to Find Dataset Page IDs

### Option 1: From Notion UI
1. Open the dataset page in Notion
2. Click "..." menu → "Copy link"
3. The page ID is in the URL: `https://www.notion.so/YourPageName-<page_id>`
4. Remove dashes: `<page_id>` should be 32 hex characters

### Option 2: From Previous Ingestion Logs
- Check ingestion logs for lines like:
  ```
  Created new dataset page 2beadf61-42ab-XXXX-XXXX-XXXXXXXXXXXX
  ```
- Or:
  ```
  Completed ingestion -> dataset page 2beadf61-42ab-XXXX-XXXX-XXXXXXXXXXXX
  ```

### Option 3: Query Notion API
Use the `list_lipidomics_datasets.py` script (may need schema adjustments).

## What Happens During Re-Ingestion

1. **File Parsing**: Re-extracts lipid species from the CSV/TSV file
2. **Feature Linking**: 
   - Creates/finds lipid species pages for each species
   - Links dataset to all species via "Experimental Data Assets" relation
3. **Signature Matching**:
   - Re-scores dataset against all active signatures
   - Updates "Related Signature(s)" relation
   - Updates "Signature Match Score" numeric property
4. **Embedding Update**:
   - Rebuilds text representation with latest species and matches
   - Updates Pinecone embeddings
   - Updates "Embedding IDs" and "Last Embedded" on dataset page

## Important Notes

- ✅ **Idempotent**: Re-running is safe - won't create duplicates
- ✅ **Non-destructive**: Only updates/creates missing links
- ✅ **Incremental**: Adds new features without removing existing ones

## Finding Your Dataset File

If you don't have the original file:
1. Check if it was attached to the Notion dataset page
2. Check your original data directory
3. Re-export from your data analysis tool (same format)
4. The file format just needs a column with lipid species names

## Expected Output

```
[INGEST][LIPIDOMICS] Starting ingestion of lipidomics file: <file>
[INGEST][LIPIDOMICS] Using existing dataset page <page_id>
[INGEST][LIPIDOMICS] Extracted N unique species from M rows
[INGEST][LIPIDOMICS] Linking N lipid species to Lipid Species DB
[INGEST][FEATURE] Added dataset <id> to lipid feature <id> (via property 'Experimental Data Assets')
...
[INGEST][LIPIDOMICS] Linked N/N lipid species to Lipid Species DB
[INGEST][LIPIDOMICS] Scoring dataset against signatures...
[INGEST][LIPIDOMICS] Found X matching signature(s)
[INGEST][LIPIDOMICS] Embedded dataset to Pinecone
[INGEST][LIPIDOMICS] Completed ingestion -> dataset page <page_id>
```

## Troubleshooting

**File not found**: Make sure the file path is correct and accessible

**Page ID not found**: Verify the page ID is correct (32 hex characters, with or without dashes)

**No species extracted**: Check that your file has a lipid species column (common names: "species", "lipid", "Lipid", "Name")

**Relations not created**: Verify the "Experimental Data Assets" relation property exists in the Lipid Species database

