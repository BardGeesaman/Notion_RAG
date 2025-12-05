# Gmail Ingestion - Quick Start

## Replace Zapier with Direct Gmail Integration

This replaces your Zapier workflow (Gmail → Notion) with direct Gmail API integration that ingests emails using Postgres-only ingestion (no Notion required).

## Quick Setup (5 minutes)

### 1. Install Dependencies
```bash
pip install google-auth google-auth-oauthlib google-auth-httplib2 google-api-python-client
```

Or use the requirements file:
```bash
pip install -r requirements_gmail.txt
```

### 2. Get Gmail API Credentials

1. Go to https://console.cloud.google.com/
2. Create/select a project
3. Enable **Gmail API**:
   - APIs & Services > Library > Search "Gmail API" > Enable
4. Create OAuth 2.0 credentials:
   - APIs & Services > Credentials > Create Credentials > OAuth client ID
   - Application type: **Desktop app**
   - Download the JSON file
5. Save as: `credentials/gmail_credentials.json`

### 3. First Run (Authentication)

```bash
python scripts/ingest_gmail.py --dry-run
```

This will:
- Open a browser window
- Ask you to sign in with `amprenta.email@gmail.com`
- Request permission to read emails
- Save the token for future use

### 4. Ingest Emails

```bash
# Ingest recent inbox emails
python scripts/ingest_gmail.py

# Ingest unread emails
python scripts/ingest_gmail.py --query "is:unread"

# Ingest emails from last 7 days
python scripts/ingest_gmail.py --days 7
```

## Automation

Set up a cron job or scheduled task:

```bash
# Runs every hour, ingests emails from last day
0 * * * * cd /path/to/project && python scripts/ingest_gmail.py --days 1
```

## Benefits vs Zapier

✅ **Faster** - Direct API access (no Zapier middleman)  
✅ **Free** - No Zapier subscription needed  
✅ **No Notion Required** - Uses Postgres-only ingestion  
✅ **More Control** - Custom queries, filtering, etc.  
✅ **Better Logging** - See exactly what's happening  

## What Happens

1. **Fetches emails** from Gmail API
2. **Ingests to Pinecone** using Postgres-only ingestion
3. **No Notion** - Completely bypassed (much faster!)
4. **Searchable** - Emails are immediately available for RAG queries

## Full Documentation

See `docs/GMAIL_INGESTION_SETUP.md` for complete setup instructions.

