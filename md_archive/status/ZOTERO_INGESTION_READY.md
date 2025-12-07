# Zotero Ingestion - Postgres-Only Script Ready! ✅

## Summary

A new Postgres-only Zotero ingestion script has been created that works **without Notion** dependencies.

## Script Location

`scripts/ingest_zotero_postgres.py`

## Features

✅ **Postgres-only** - No Notion API calls  
✅ **Direct to Pinecone** - Fast, scalable ingestion  
✅ **Processes**:
   - Item metadata (title, authors, DOI, etc.)
   - Abstracts
   - Notes
   - Attachments (PDFs, text files)
✅ **Automatic text extraction** from PDFs and other attachments

## Usage

### Prerequisites

You need a Zotero collection key. You can find it by:

1. **From Zotero Desktop**:
   - Right-click on a collection
   - Select "Properties" or check the collection URL
   - The key is usually 8 characters (e.g., `3RGXZTAY`)

2. **From Zotero API**:
   ```bash
   # List all collections (requires Zotero API access)
   curl -H "Zotero-API-Key: YOUR_KEY" \
        https://api.zotero.org/groups/YOUR_GROUP_ID/collections
   ```

### Basic Commands

#### 1. Dry Run (Recommended First Step)

```bash
python scripts/ingest_zotero_postgres.py \
    --collection-key YOUR_COLLECTION_KEY \
    --dry-run \
    --max-items 2
```

This will show what would be ingested without actually doing it.

#### 2. Ingest a Few Items (Test)

```bash
python scripts/ingest_zotero_postgres.py \
    --collection-key YOUR_COLLECTION_KEY \
    --max-items 2
```

This will ingest 2 items from the collection.

#### 3. Ingest All Items

```bash
python scripts/ingest_zotero_postgres.py \
    --collection-key YOUR_COLLECTION_KEY
```

This will ingest all literature items in the collection.

## What Gets Ingested

For each Zotero item, the script:

1. ✅ Fetches item metadata (title, authors, DOI, journal, etc.)
2. ✅ Extracts abstract
3. ✅ Processes notes (converts HTML to text)
4. ✅ Downloads and extracts text from attachments (PDFs, etc.)
5. ✅ Combines all content
6. ✅ Chunks and embeds to Pinecone
7. ✅ Includes rich metadata (authors, DOI, tags, etc.)

## Supported Item Types

The script processes these Zotero item types:
- `journalArticle`
- `book`
- `bookSection`
- `report`
- `preprint`
- `conferencePaper`

Other types are automatically skipped.

## Architecture

- ✅ **No Notion** - Completely Postgres-only
- ✅ **Direct to Pinecone** - Fast ingestion
- ✅ **Feature extraction** - Automatically extracts biological features
- ✅ **Rich metadata** - Authors, DOI, tags, journal, etc.

## Performance

- **Fast**: Direct database operations
- **Scalable**: Handles large collections
- **Efficient**: Only processes literature items

## Next Steps

1. **Get a collection key** from your Zotero library
2. **Test with dry-run** to see what would be ingested
3. **Ingest a few items** to verify it works
4. **Ingest full collection** when ready

## Example

```bash
# First, do a dry run to see what's there
python scripts/ingest_zotero_postgres.py \
    --collection-key 3RGXZTAY \
    --dry-run

# Then ingest 2 items as a test
python scripts/ingest_zotero_postgres.py \
    --collection-key 3RGXZTAY \
    --max-items 2

# Finally, ingest everything
python scripts/ingest_zotero_postgres.py \
    --collection-key 3RGXZTAY
```

## Troubleshooting

### "Collection not found"
- Verify the collection key is correct
- Check that the collection is in your Zotero library
- Ensure your Zotero API key has access

### "No items found"
- Check that the collection has literature items (journal articles, books, etc.)
- Other item types (notes, standalone attachments) are skipped

### "Insufficient content"
- Items with very short abstracts/notes may be skipped
- Try items with PDF attachments for more content

## Ready to Test!

The script is ready to use. Do you have a Zotero collection key you'd like to test with?

