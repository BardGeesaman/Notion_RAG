# Gmail Ingestion Setup Guide

## Overview

This guide explains how to set up direct Gmail ingestion, replacing the Zapier → Notion workflow. Emails are fetched directly from Gmail and ingested using Postgres-only ingestion (no Notion required).

## Prerequisites

1. **Python packages**:
   ```bash
   pip install google-auth google-auth-oauthlib google-auth-httplib2 google-api-python-client
   ```

2. **Google Cloud Project** with Gmail API enabled

## Setup Steps

### Step 1: Create Google Cloud Project

1. Go to [Google Cloud Console](https://console.cloud.google.com/)
2. Create a new project (or select existing)
3. Enable the **Gmail API**:
   - Go to "APIs & Services" > "Library"
   - Search for "Gmail API"
   - Click "Enable"

### Step 2: Create OAuth 2.0 Credentials

1. Go to "APIs & Services" > "Credentials"
2. Click "Create Credentials" > "OAuth client ID"
3. If prompted, configure OAuth consent screen:
   - User Type: "External" (or "Internal" if using Google Workspace)
   - App name: "Amprenta Email Ingestor"
   - User support email: your email
   - Scopes: Add `https://www.googleapis.com/auth/gmail.readonly`
   - Save
4. Application type: **"Desktop app"**
5. Name: "Amprenta Gmail Client"
6. Click "Create"
7. Download the credentials JSON file

### Step 3: Save Credentials

1. Create a `credentials` directory in the project root:
   ```bash
   mkdir -p credentials
   ```

2. Move the downloaded credentials file to:
   ```
   credentials/gmail_credentials.json
   ```

### Step 4: Authenticate (First Time)

Run the ingestion script. On first run, it will:
1. Open a browser window
2. Ask you to sign in with your Gmail account (`amprenta.email@gmail.com`)
3. Ask for permission to read emails
4. Save the token for future use

```bash
python scripts/ingest_gmail.py --dry-run
```

The token will be saved to `credentials/gmail_token.json` for future use.

### Step 5: Configure (Optional)

Add to your `.env` file (optional - defaults work):

```bash
# Gmail API credentials (optional - defaults shown)
GMAIL_CREDENTIALS_FILE=credentials/gmail_credentials.json
GMAIL_TOKEN_FILE=credentials/gmail_token.json
```

## Usage

### Basic Usage

```bash
# Ingest recent inbox emails (default: last 100)
python scripts/ingest_gmail.py

# Ingest unread emails
python scripts/ingest_gmail.py --query "is:unread"

# Ingest emails from specific sender
python scripts/ingest_gmail.py --query "from:sender@example.com"

# Ingest emails from last 7 days
python scripts/ingest_gmail.py --days 7

# Ingest all inbox emails (up to 1000)
python scripts/ingest_gmail.py --all

# Dry run (see what would be ingested)
python scripts/ingest_gmail.py --dry-run
```

### Advanced Queries

Use any Gmail search query:
- `in:inbox` - All inbox emails
- `is:unread` - Unread emails
- `from:example@gmail.com` - From specific sender
- `subject:important` - Subject contains "important"
- `after:2024/01/01` - After date
- `has:attachment` - With attachments
- Combine: `in:inbox is:unread from:example@gmail.com`

### Automation (Recommended)

Set up a cron job or scheduled task to run automatically:

```bash
# Add to crontab (runs every hour)
0 * * * * cd /path/to/project && python scripts/ingest_gmail.py --days 1
```

Or use systemd timer, GitHub Actions, etc.

## How It Works

1. **Gmail API** - Fetches emails directly from Gmail
2. **Postgres-Only Ingestion** - Uses `ingest_email_content()` (no Notion)
3. **Pinecone** - Emails are embedded and stored for RAG queries
4. **No Notion** - Completely bypasses Notion (much faster!)

## Comparison: Zapier vs Direct Gmail

| Feature | Zapier → Notion | Direct Gmail |
|---------|----------------|--------------|
| Speed | Slow (multiple API calls) | Fast (direct) |
| Cost | Zapier subscription | Free |
| Notion Required | ✅ Yes | ❌ No |
| Reliability | Dependent on Zapier | Direct connection |
| Setup | Zapier webhook | OAuth2 (one-time) |

## Troubleshooting

### "Gmail API libraries not installed"
```bash
pip install google-auth google-auth-oauthlib google-auth-httplib2 google-api-python-client
```

### "Credentials file not found"
- Download OAuth2 credentials from Google Cloud Console
- Save as `credentials/gmail_credentials.json`

### "Token expired"
- Delete `credentials/gmail_token.json`
- Run script again to re-authenticate

### "Permission denied"
- Make sure Gmail API is enabled in Google Cloud Console
- Check OAuth consent screen is configured
- Verify scopes include `gmail.readonly`

## Security Notes

- Credentials files contain sensitive data - **do not commit to git**
- Add to `.gitignore`:
  ```
  credentials/gmail_credentials.json
  credentials/gmail_token.json
  ```
- Token file is user-specific and should be kept secure
- OAuth2 only grants read-only access to emails

## Next Steps

1. Set up authentication (Step 4)
2. Test with `--dry-run`
3. Run first ingestion
4. Set up automated scheduling
5. Monitor logs for any issues

## Support

For issues or questions:
- Check logs for detailed error messages
- Verify Gmail API is enabled
- Ensure OAuth2 credentials are correct
- Check network connectivity

