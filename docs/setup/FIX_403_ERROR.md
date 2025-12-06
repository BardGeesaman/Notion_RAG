# Fix Error 403: Access Denied

## The Problem

You're getting "Error 403: access_denied" when trying to authenticate. This means the OAuth consent screen needs to be configured properly.

## Solution: Configure OAuth Consent Screen

### Step 1: Go to OAuth Consent Screen

**Direct link**: https://console.cloud.google.com/apis/credentials/consent

Or navigate:
1. Go to Google Cloud Console
2. APIs & Services > OAuth consent screen

### Step 2: Configure the Consent Screen

If you see "CONFIGURE CONSENT SCREEN" or "EDIT APP":
1. Click it
2. Select **"External"** as user type (unless using Google Workspace)
3. Fill in:
   - **App name**: `Amprenta Email Ingestor`
   - **User support email**: Your email
   - **Developer contact**: Your email
4. Click "SAVE AND CONTINUE"

### Step 3: Add Scopes

1. Click "**Scopes**" tab (or "SAVE AND CONTINUE" to get to Scopes)
2. Click "**+ ADD OR REMOVE SCOPES**"
3. Search for: `gmail.readonly`
4. Check the box for: `https://www.googleapis.com/auth/gmail.readonly`
5. Click "**UPDATE**"
6. Click "**SAVE AND CONTINUE**"

### Step 4: Add Test User (IMPORTANT!)

1. Click "**Test users**" tab
2. Click "**+ ADD USERS**" button
3. Enter: `amprenta.email@gmail.com`
4. Click "**ADD**"
5. You should see `amprenta.email@gmail.com` in the test users list
6. Click "**SAVE AND CONTINUE**"

### Step 5: Complete

1. Review the summary
2. Click "**BACK TO DASHBOARD**"

## Try Again

After configuring the consent screen and adding the test user, try running:

```bash
python scripts/ingest_gmail.py --dry-run
```

This should now work!

## Common Issues

**"I don't see Test users tab"**
- Make sure you selected "External" as user type
- Test users are only for External apps

**"Can't add test user"**
- Make sure you're on the "Test users" tab
- Try refreshing the page

**"Still getting 403 error"**
- Make sure you added `amprenta.email@gmail.com` as a test user
- Try signing out and signing in again with that account

