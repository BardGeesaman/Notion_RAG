# Get Gmail Credentials - Quick Start

## 🚀 Fast Track (5 minutes)

### Step 1: Open Google Cloud Console
**Click this link**: https://console.cloud.google.com/

### Step 2: Quick Links (Copy/Paste)

1. **Enable Gmail API**: https://console.cloud.google.com/apis/library/gmail.googleapis.com
   - Click "ENABLE" button

2. **OAuth Consent Screen**: https://console.cloud.google.com/apis/credentials/consent
   - User Type: "External"
   - App name: "Amprenta Email Ingestor"
   - Scopes: Add `gmail.readonly`
   - Test users: Add `amprenta.email@gmail.com`

3. **Create Credentials**: https://console.cloud.google.com/apis/credentials
   - Click "+ CREATE CREDENTIALS" > "OAuth client ID"
   - Application type: **"Desktop app"**
   - Name: "Amprenta Gmail Client"
   - Click "CREATE"
   - **Download the JSON file** (click 📥 icon)

### Step 3: Save the File

1. The downloaded file has a name like `client_secret_XXXXXX.json`
2. Rename it to: `gmail_credentials.json`
3. Move it to: `/Users/bard/Documents/Notion RAG/credentials/`

### Step 4: Verify

Run this command to check:
```bash
python scripts/verify_gmail_credentials.py
```

### Step 5: Authenticate

Once verified, run:
```bash
python scripts/setup_gmail.py
```

This will open your browser to authenticate with `amprenta.email@gmail.com`.

---

## 📖 Detailed Instructions

For step-by-step screenshots and detailed instructions, see:
- **`GMAIL_CREDENTIALS_GUIDE.md`** - Complete walkthrough

## ⚡ What You Need

- Google account (for `amprenta.email@gmail.com`)
- 5 minutes
- Browser access

## 🎯 Goal

Get a JSON file from Google Cloud Console and save it as:
```
credentials/gmail_credentials.json
```

That's it! No passwords, no API keys in .env - just this one file.

