# Gmail Setup Steps - Quick Guide

## Important: No Password Needed! 🎉

**Gmail API uses OAuth2 (not passwords)**. You just need to authenticate once, and then it saves a token for future use.

## Step-by-Step Setup

### Step 1: Get OAuth2 Credentials from Google Cloud Console

1. **Go to Google Cloud Console**:
   - Open: https://console.cloud.google.com/
   - Sign in with your Google account

2. **Create or Select a Project**:
   - Click project dropdown at top
   - Create new project: "Amprenta Email Ingestor"
   - Or select existing project

3. **Enable Gmail API**:
   - Go to: APIs & Services > Library
   - Search for: "Gmail API"
   - Click "Enable"

4. **Create OAuth 2.0 Credentials**:
   - Go to: APIs & Services > Credentials
   - Click: "+ CREATE CREDENTIALS" > "OAuth client ID"
   - **First time?** Configure OAuth consent screen:
     - User Type: "External" (or "Internal" if using Google Workspace)
     - App name: "Amprenta Email Ingestor"
     - User support email: your email (amprenta.email@gmail.com)
     - Developer contact: your email
     - Click "Save and Continue"
     - Scopes: Click "+ ADD OR REMOVE SCOPES"
       - Search and select: `.../auth/gmail.readonly`
       - Click "Update" > "Save and Continue"
     - Test users: Add `amprenta.email@gmail.com`
     - Click "Save and Continue"
   - Application type: **"Desktop app"**
   - Name: "Amprenta Gmail Client"
   - Click "CREATE"
   - **Download the JSON file** (click download icon)

5. **Save Credentials**:
   - Save the downloaded file as: `credentials/gmail_credentials.json`
   - Make sure it's in the project root, in a `credentials/` folder

### Step 2: Authenticate (First Time Only)

Run the setup script - it will:
1. Open a browser window
2. Ask you to sign in with `amprenta.email@gmail.com`
3. Request permission to read emails
4. Save the token automatically

```bash
python scripts/setup_gmail.py
```

Or directly test:

```bash
python scripts/ingest_gmail.py --dry-run
```

### Step 3: Test Ingestion

Once authenticated, you can ingest emails:

```bash
# Dry run (see what would be ingested)
python scripts/ingest_gmail.py --dry-run

# Ingest recent inbox emails
python scripts/ingest_gmail.py

# Ingest unread emails
python scripts/ingest_gmail.py --query "is:unread"

# Ingest last 7 days
python scripts/ingest_gmail.py --days 7
```

## What Gets Created

- `credentials/gmail_credentials.json` - OAuth2 credentials (from Google Cloud Console)
- `credentials/gmail_token.json` - Authentication token (auto-created on first run)

Both files are automatically ignored by git (already in .gitignore).

## Optional: Environment Variables

You can optionally add to `.env` (but defaults work fine):

```bash
# Optional - defaults shown
GMAIL_CREDENTIALS_FILE=credentials/gmail_credentials.json
GMAIL_TOKEN_FILE=credentials/gmail_token.json
```

**No password needed!** OAuth2 handles authentication securely.

## Troubleshooting

- **"Credentials file not found"**: Make sure you downloaded and saved the OAuth2 credentials JSON file
- **"Permission denied"**: Make sure Gmail API is enabled and OAuth consent screen is configured
- **"Token expired"**: Delete `credentials/gmail_token.json` and run again to re-authenticate

## Next Steps

Once setup is complete, you can:
1. Set up automated ingestion (cron job, scheduled task)
2. Configure custom email filters/queries
3. Monitor ingestion logs

See `docs/GMAIL_INGESTION_SETUP.md` for full documentation.

