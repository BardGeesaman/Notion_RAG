# Get Gmail API Credentials - Step-by-Step Guide

Follow these steps exactly to get your Gmail API credentials.

## Step 1: Open Google Cloud Console

1. Open your web browser
2. Go to: **https://console.cloud.google.com/**
3. Sign in with your Google account (the one that has access to `amprenta.email@gmail.com`)

## Step 2: Create or Select a Project

**Option A: Create New Project (Recommended)**
1. Click the project dropdown at the top of the page (next to "Google Cloud")
2. Click **"+ Create Project"**
3. Project name: `Amprenta Email Ingestor`
4. Click **"Create"**
5. Wait for project to be created (few seconds)
6. Make sure the new project is selected (check dropdown at top)

**Option B: Use Existing Project**
1. Click the project dropdown at the top
2. Select an existing project you want to use

## Step 3: Enable Gmail API

1. In the left sidebar, click **"APIs & Services"** > **"Library"**
   - Or go directly to: https://console.cloud.google.com/apis/library
2. In the search box at the top, type: **"Gmail API"**
3. Click on **"Gmail API"** in the search results
4. Click the blue **"ENABLE"** button
5. Wait for it to enable (may take a few seconds)
6. You should see a green checkmark or "API enabled" message

## Step 4: Configure OAuth Consent Screen (First Time Only)

**Skip this step if you've already configured OAuth consent screen before.**

1. In the left sidebar, click **"APIs & Services"** > **"OAuth consent screen"**
   - Or go directly to: https://console.cloud.google.com/apis/credentials/consent
2. Select **"External"** (unless you're using Google Workspace, then use "Internal")
3. Click **"CREATE"**
4. Fill in the form:
   - **App name**: `Amprenta Email Ingestor`
   - **User support email**: Select your email (or enter `amprenta.email@gmail.com`)
   - **Developer contact information**: Enter your email
   - Click **"SAVE AND CONTINUE"**
5. **Scopes** page:
   - Click **"+ ADD OR REMOVE SCOPES"**
   - Scroll down and find: **`https://www.googleapis.com/auth/gmail.readonly`**
   - Check the box next to it
   - Click **"UPDATE"** at the bottom
   - Click **"SAVE AND CONTINUE"**
6. **Test users** page (if using External):
   - Click **"+ ADD USERS"**
   - Enter: `amprenta.email@gmail.com`
   - Click **"ADD"**
   - Click **"SAVE AND CONTINUE"**
7. **Summary** page:
   - Review and click **"BACK TO DASHBOARD"**

## Step 5: Create OAuth 2.0 Credentials

1. In the left sidebar, click **"APIs & Services"** > **"Credentials"**
   - Or go directly to: https://console.cloud.google.com/apis/credentials
2. At the top, click **"+ CREATE CREDENTIALS"**
3. Click **"OAuth client ID"**
4. **Application type**: Select **"Desktop app"**
5. **Name**: `Amprenta Gmail Client`
6. Click **"CREATE"**
7. A popup will appear showing your credentials:
   - **Client ID** (long string)
   - **Client secret** (long string)
   - **Download JSON** button (📥 icon)
8. **Click the download icon** (📥) to download the JSON file
   - Or click **"OK"** and download it from the credentials list

## Step 6: Save the Credentials File

1. The downloaded file will be named something like: `client_secret_XXXXXX.json`
2. **Rename it to**: `gmail_credentials.json`
3. **Move it to**: `/Users/bard/Documents/Notion RAG/credentials/gmail_credentials.json`
   - Create the `credentials` folder if it doesn't exist
   - Make sure the file is named exactly: `gmail_credentials.json`

## Step 7: Verify the File

The file should be located at:
```
/Users/bard/Documents/Notion RAG/credentials/gmail_credentials.json
```

You can verify it exists:
```bash
ls -la credentials/gmail_credentials.json
```

The file should contain JSON with fields like:
- `installed.client_id`
- `installed.client_secret`
- `installed.auth_uri`
- etc.

## Quick Checklist

- [ ] Signed in to Google Cloud Console
- [ ] Created/selected a project
- [ ] Enabled Gmail API
- [ ] Configured OAuth consent screen (if first time)
- [ ] Created OAuth 2.0 Desktop app credentials
- [ ] Downloaded the credentials JSON file
- [ ] Renamed it to `gmail_credentials.json`
- [ ] Saved it to `credentials/gmail_credentials.json`

## Next Steps

Once you've saved the credentials file, run:

```bash
python scripts/setup_gmail.py
```

This will test the connection and authenticate you (opens browser).

## Common Issues

**"Can't find APIs & Services"**
- Make sure you're signed in
- Check that you have a project selected

**"Can't enable Gmail API"**
- Make sure billing is enabled (Gmail API requires billing, but it's free for this use case)
- Or try refreshing the page

**"OAuth consent screen shows errors"**
- Make sure you filled all required fields
- For External apps, you may need to verify your domain (can skip for testing)

**"Credentials file not found"**
- Check the exact path: `credentials/gmail_credentials.json`
- Make sure the `credentials` folder exists
- Check file name is exactly `gmail_credentials.json` (not `.txt` or anything else)

## Need Help?

If you get stuck at any step, let me know which step and what error you see!

