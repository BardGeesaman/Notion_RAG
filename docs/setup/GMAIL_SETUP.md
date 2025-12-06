# Gmail API Setup Guide

Complete guide for setting up Gmail API integration for email ingestion.

## Overview

Gmail API uses OAuth2 authentication (no passwords required). You authenticate once, and the system saves a token for future use.

## Quick Start (5 minutes)

### 1. Install Dependencies

```bash
pip install google-auth google-auth-oauthlib google-auth-httplib2 google-api-python-client
```

Or use the requirements file:
```bash
pip install -r requirements_gmail.txt
```

### 2. Get Gmail API Credentials

1. **Go to Google Cloud Console**: https://console.cloud.google.com/
2. **Create or Select a Project**:
   - Click project dropdown at top
   - Create new project: "Amprenta Email Ingestor"
   - Or select existing project

3. **Enable Gmail API**:
   - Go to: APIs & Services > Library
   - Search for: "Gmail API"
   - Click "Enable"

4. **Configure OAuth Consent Screen** (first time only):
   - Go to: APIs & Services > OAuth consent screen
   - User Type: "External" (or "Internal" if using Google Workspace)
   - App name: "Amprenta Email Ingestor"
   - User support email: your email
   - Developer contact: your email
   - Click "Save and Continue"
   - Scopes: Click "+ ADD OR REMOVE SCOPES"
     - Search and select: `.../auth/gmail.readonly`
     - Click "Update" > "Save and Continue"
   - Test users: Add your Gmail address
   - Click "Save and Continue"

5. **Create OAuth 2.0 Credentials**:
   - Go to: APIs & Services > Credentials
   - Click: "+ CREATE CREDENTIALS" > "OAuth client ID"
   - Application type: **"Desktop app"**
   - Name: "Amprenta Gmail Client"
   - Click "CREATE"
   - **Download the JSON file** (click download icon)

6. **Save Credentials**:
   - Save the downloaded file as: `credentials/gmail_credentials.json`
   - Make sure it's in the project root, in a `credentials/` folder

### 3. First Run (Authentication)

```bash
python scripts/ingest_gmail.py --dry-run
```

This will:
- Open a browser window
- Ask you to sign in with your Gmail account
- Request permission to read emails
- Save the token for future use (`credentials/gmail_token.json`)

### 4. Ingest Emails

```bash
# Ingest recent inbox emails
python scripts/ingest_gmail.py

# Ingest unread emails
python scripts/ingest_gmail.py --query "is:unread"

# Ingest emails from specific sender
python scripts/ingest_gmail.py --query "from:example@domain.com"
```

## Troubleshooting

### 403 Error: Access Denied

- Make sure you've added your email as a test user in OAuth consent screen
- Check that Gmail API is enabled in your project
- Verify the credentials file is in the correct location

### Token Expired

Delete `credentials/gmail_token.json` and run the setup again to re-authenticate.

### Can't Find Credentials

Make sure `credentials/gmail_credentials.json` exists in the project root directory.

## Additional Resources

For detailed troubleshooting, see:
- `docs/setup/GMAIL_CREDENTIALS_GUIDE.md` - Detailed credential setup
- `docs/setup/GMAIL_SETUP_WALKTHROUGH.md` - Step-by-step walkthrough
- `docs/GMAIL_INGESTION_SETUP.md` - Full ingestion guide

## Related Documentation

- [Email Ingestion Guide](../docs/EMAIL_INGESTION.md) (if exists)
- [Configuration Guide](../docs/CONFIGURATION.md)

