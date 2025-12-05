# Gmail Setup Walkthrough - Step by Step

## After Logging In to Google Cloud Console

### Step 1: Create or Select a Project

1. Look at the top of the page - you'll see a project dropdown (next to "Google Cloud" logo)
2. **Click on the project dropdown**
3. Choose one:
   - **Option A**: Click "**+ NEW PROJECT**" 
     - Project name: `Amprenta Email Ingestor`
     - Click "**CREATE**"
     - Wait a few seconds for it to create
   - **Option B**: Select an existing project
4. Make sure your project is selected (check the dropdown at the top)

### Step 2: Enable Gmail API

**Option 1: Using the Search Box**
1. At the top of the page, you'll see a search box that says "Search products and resources"
2. Type: **"Gmail API"**
3. Click on "Gmail API" in the search results
4. Click the blue **"ENABLE"** button
5. Wait a few seconds - you'll see "API enabled" or a checkmark

**Option 2: Direct Navigation**
1. In the left sidebar, click **"APIs & Services"**
2. Click **"Library"** 
3. Search for "Gmail API"
4. Click "Gmail API"
5. Click **"ENABLE"**

### Step 3: Configure OAuth Consent Screen (First Time Only)

1. In the left sidebar, click **"APIs & Services"**
2. Click **"OAuth consent screen"**
3. If you see "OAuth consent screen is not configured", click the button to configure it
4. Fill in:
   - **User Type**: Select "**External**" (unless you're using Google Workspace)
   - Click "**CREATE**"
5. **App information** tab:
   - **App name**: `Amprenta Email Ingestor`
   - **User support email**: Select your email from dropdown (or enter `amprenta.email@gmail.com`)
   - **Developer contact information**: Enter your email
   - Click "**SAVE AND CONTINUE**" at the bottom
6. **Scopes** tab:
   - Click "**+ ADD OR REMOVE SCOPES**"
   - In the search box, type: `gmail.readonly`
   - Find: `https://www.googleapis.com/auth/gmail.readonly`
   - Check the box next to it
   - Click "**UPDATE**" at the bottom
   - Click "**SAVE AND CONTINUE**" at the bottom
7. **Test users** tab:
   - Click "**+ ADD USERS**"
   - Enter: `amprenta.email@gmail.com`
   - Click "**ADD**"
   - Click "**SAVE AND CONTINUE**" at the bottom
8. **Summary** tab:
   - Review the information
   - Click "**BACK TO DASHBOARD**"

### Step 4: Create OAuth 2.0 Credentials

1. In the left sidebar, click **"APIs & Services"**
2. Click **"Credentials"**
3. At the top of the page, click **"+ CREATE CREDENTIALS"**
4. From the dropdown, click **"OAuth client ID"**
5. **Application type**: Select **"Desktop app"** from the dropdown
6. **Name**: Type `Amprenta Gmail Client`
7. Click "**CREATE**" button at the bottom
8. A popup window will appear showing:
   - **Client ID** (a long string)
   - **Client secret** (a long string)
   - A **download icon** (📥)
9. **Click the download icon** (📥) to download the JSON file
   - The file will download to your Downloads folder
   - The file name will be something like: `client_secret_123456789-abcdefg.json`

### Step 5: Save the Credentials File

1. **Find the downloaded file** in your Downloads folder
   - It will have a name like: `client_secret_123456789-abcdefg.json`
2. **Rename it** to: `gmail_credentials.json`
3. **Move it** to: `/Users/bard/Documents/Notion RAG/credentials/`
   - You can drag and drop it into the credentials folder
   - Or use Finder to navigate and move it

### Step 6: Verify

Once you've saved the file, let me know and I'll verify it's set up correctly!

Or run this command to check:
```bash
python scripts/verify_gmail_credentials.py
```

## Visual Guide

### What You'll See:

1. **Top Bar**: Project dropdown, search box
2. **Left Sidebar**: Navigation menu (APIs & Services, etc.)
3. **Main Area**: Content that changes based on what you click

### Key Buttons to Look For:

- **"ENABLE"** - Blue button to enable APIs
- **"+ CREATE CREDENTIALS"** - Button to create OAuth credentials
- **"CREATE"** - Button to create projects/credentials
- **Download icon (📥)** - To download the credentials file

## Troubleshooting

**Can't find "APIs & Services"?**
- Look in the left sidebar menu (hamburger menu ≡ if it's collapsed)
- Or use the search box at the top to search for "APIs"

**Gmail API already enabled?**
- That's fine! Skip to Step 3 or 4

**OAuth consent screen already configured?**
- Great! Skip to Step 4

**Stuck?**
- Tell me which step you're on and what you see on your screen
- I can help you find the right buttons/options

