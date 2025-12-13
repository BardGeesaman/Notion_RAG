# Idempotent Ingestion - Incremental Updates

## Overview

All ingestion pipelines now support **idempotent ingestion** - they only process items that are new or have changed. This matches the behavior of the old Notion-based system and enables efficient incremental updates.

## How It Works

### Content Hashing

Each piece of content gets a unique identifier (`content_id`) and a hash (`content_hash`):
- **content_id**: Unique identifier (e.g., Gmail message ID, Zotero item key)
- **content_hash**: SHA256 hash of the content text

### Idempotency Check

Before ingesting, the system:

1. ✅ **Checks if content exists** in Pinecone using `content_id`
2. ✅ **Compares content hash** if content exists
3. ✅ **Skips if unchanged** (same hash = no changes)
4. ✅ **Re-ingests if changed** (deletes old vectors, ingests new content)
5. ✅ **Ingests if new** (doesn't exist in Pinecone)

### Change Detection

- **Same hash** → Content unchanged → Skip ingestion
- **Different hash** → Content changed → Delete old, ingest new
- **Not found** → New content → Ingest

## Supported Ingestion Types

### ✅ Email Ingestion (Gmail)

**Idempotency Key**: Gmail message ID (`gmail_{message_id}`)

```bash
# Only ingest new/changed emails
python scripts/ingest_gmail.py --query "is:unread"

# Force re-ingest all emails
python scripts/ingest_gmail.py --query "in:inbox" --force
```

### ✅ Zotero/Literature Ingestion

**Idempotency Key**: Zotero item key

```bash
# Only ingest new/changed items
python scripts/ingest_zotero_postgres.py --collection-key 3RGXZTAY

# Force re-ingest all items
python scripts/ingest_zotero_postgres.py --collection-key 3RGXZTAY --force
```

## Command-Line Options

### `--force` Flag

Use `--force` to re-ingest everything, even if already ingested:

```bash
# Gmail
python scripts/ingest_gmail.py --force

# Zotero
python scripts/ingest_zotero_postgres.py --collection-key KEY --force
```

**When to use `--force`:**
- Testing/debugging
- After fixing ingestion bugs
- When you want to refresh all content

## Example Output

### First Run (New Content)

```
[GMAIL] [1/23] Processing: Important email
[GMAIL] ✅ Ingested email: Important email (5 chunks)
```

### Second Run (Unchanged Content)

```
[GMAIL] [1/23] Processing: Important email
[GMAIL] ⏭️  Skipped email: Important email (already ingested)
```

### Third Run (Changed Content)

```
[GMAIL] [1/23] Processing: Important email (updated)
[INGEST] Content has changed; deleting old vectors and re-ingesting
[GMAIL] ✅ Ingested email: Important email (updated) (6 chunks)
```

## Implementation Details

### Content Hash Storage

Content hash is stored in Pinecone metadata:
```python
{
    "content_id": "gmail_19",
    "content_hash": "abc123...",  # SHA256 hash
    "source_type": "email",
    "title": "...",
    ...
}
```

### Query Performance

Idempotency checks use Pinecone metadata filters:
- Fast queries by `content_id`
- No full-text search required
- Efficient change detection

## Benefits

### Performance

- ✅ **Faster ingestion** - Only processes new/changed items
- ✅ **Reduced API calls** - Skips unchanged content
- ✅ **Lower costs** - Less embedding generation

### Reliability

- ✅ **No duplicates** - Same content only ingested once
- ✅ **Change tracking** - Automatically detects updates
- ✅ **Safe re-runs** - Can re-run scripts without duplicating

### Efficiency

- ✅ **Incremental updates** - Only process what's new
- ✅ **Batch processing** - Handle large collections efficiently
- ✅ **Automation-friendly** - Safe to run on schedule

## Usage Patterns

### Daily Incremental Updates

```bash
# Cron job: Only process new emails
0 * * * * python scripts/ingest_gmail.py --query "is:unread"

# Daily: Only process new Zotero items
0 2 * * * python scripts/ingest_zotero_postgres.py --collection-key KEY
```

### Manual Full Refresh

```bash
# Re-ingest everything (use sparingly)
python scripts/ingest_gmail.py --all --force
python scripts/ingest_zotero_postgres.py --collection-key KEY --force
```

## Summary

✅ **Idempotent ingestion is now enabled** for:
- Email ingestion (Gmail)
- Literature ingestion (Zotero)

✅ **Only new/changed content** is processed  
✅ **Automatic change detection** via content hashing  
✅ **Force flag available** for full re-ingestion  

This matches the behavior of the old Notion-based system while providing better performance and scalability!

