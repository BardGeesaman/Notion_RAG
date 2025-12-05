# Fix Gmail OAuth 403 Error - Do This Now

## The Problem
You're getting `Error 403: access_denied` because the OAuth consent screen isn't configured yet. The "Test users" tab won't appear until you configure the consent screen first.

## What To Do Right Now

### Step 1: Open OAuth Consent Screen
Go to this exact URL:
```
https://console.cloud.google.com/apis/credentials/consent
```

### Step 2: What Do You See?
Look at the **center of the page**. Tell me which of these you see:

**Option A:** A big button that says **"CONFIGURE CONSENT SCREEN"**
- If you see this → Click it and go to Step 3

**Option B:** A button that says **"EDIT APP"**
- If you see this → You've started configuring. Go to Step 3

**Option C:** Tabs at the top like:
  - App information
  - Scopes
  - Test users
- If you see tabs → Go to Step 4

**Option D:** Something else entirely?
- Tell me what you see

### Step 3: Initial Setup (If you clicked "CONFIGURE CONSENT SCREEN")
1. Select **"External"** as user type
2. Click **"CREATE"**
3. You should now see tabs at the top

### Step 4: Fill Out "App information" Tab
1. Click the **"App information"** tab (if not already selected)
2. Fill in:
   - **App name**: `Amprenta Email Ingestor`
   - **User support email**: Select `amprenta.email@gmail.com` or your email
   - **App logo**: (Skip - optional)
   - **App domain**: (Skip - optional)
   - **Application home page**: (Skip - optional)
   - **Application privacy policy link**: (Skip - optional)
   - **Application terms of service link**: (Skip - optional)
   - **Authorized domains**: (Skip - optional)
   - **Developer contact information**: Enter your email
3. Click **"SAVE AND CONTINUE"** at the bottom

### Step 5: Add Scopes
1. Click the **"Scopes"** tab
2. Click **"+ ADD OR REMOVE SCOPES"** button
3. In the filter/search box, type: `gmail.readonly`
4. Check the box next to: `https://www.googleapis.com/auth/gmail.readonly`
5. Click **"UPDATE"** button
6. Click **"SAVE AND CONTINUE"** at the bottom

### Step 6: Add Test Users (THIS IS THE KEY STEP!)
1. Click the **"Test users"** tab (it should appear now!)
2. Click **"+ ADD USERS"** button
3. Type: `amprenta.email@gmail.com`
4. Click **"ADD"**
5. You should see `amprenta.email@gmail.com` in the list
6. Click **"SAVE AND CONTINUE"** at the bottom

### Step 7: Complete
1. Review the summary
2. Click **"BACK TO DASHBOARD"** or just close the tab

## Important Notes

- **You MUST add `amprenta.email@gmail.com` as a test user** or you'll keep getting the 403 error
- The "Test users" tab only appears after you've configured at least the App information
- For External apps, you can only use test users until you publish the app (we don't need to publish)

## After Configuration

Once you've added the test user, try running the Gmail ingestion again:
```bash
python scripts/ingest_gmail.py --dry-run
```

The authentication should work now!

## Still Stuck?

Tell me:
1. What do you see when you go to https://console.cloud.google.com/apis/credentials/consent?
2. Which step are you on?
3. What's blocking you?

