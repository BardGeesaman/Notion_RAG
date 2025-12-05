# Configure OAuth Consent Screen - Step by Step

## The Issue

The "Test users" tab only appears **after** you start configuring the OAuth consent screen. Let's configure it first!

## Step-by-Step Configuration

### Step 1: Go to OAuth Consent Screen

**Direct link**: https://console.cloud.google.com/apis/credentials/consent

### Step 2: Start Configuration

Look at the **main content area** (center of the page):

**If you see "CONFIGURE CONSENT SCREEN" button:**
1. Click it
2. Select **"External"** as user type
3. Click **"CREATE"**

**If you see "EDIT APP" button:**
1. Click it
2. You're already in configuration mode

**If you see configuration tabs already:**
- You're already configured! Skip to Step 3

### Step 3: Fill Out App Information Tab

1. **App name**: `Amprenta Email Ingestor`
2. **User support email**: Select your email or enter `amprenta.email@gmail.com`
3. **Developer contact information**: Enter your email
4. Click **"SAVE AND CONTINUE"** at the bottom

### Step 4: Add Scopes Tab

1. Click **"Scopes"** tab (or it will appear after Step 3)
2. Click **"+ ADD OR REMOVE SCOPES"** button
3. In the search box, type: `gmail.readonly`
4. Find and check: `https://www.googleapis.com/auth/gmail.readonly`
5. Click **"UPDATE"** button
6. Click **"SAVE AND CONTINUE"**

### Step 5: Add Test Users Tab (This Will Appear Now!)

1. Click **"Test users"** tab (should be visible now!)
2. Click **"+ ADD USERS"** button
3. Enter: `amprenta.email@gmail.com`
4. Click **"ADD"**
5. You should see the email in the list
6. Click **"SAVE AND CONTINUE"**

### Step 6: Complete

1. Review the summary
2. Click **"BACK TO DASHBOARD"**

## What You Should See

After clicking "CONFIGURE CONSENT SCREEN" and selecting "External", you should see tabs at the top:
- App information
- Scopes  
- Test users
- Summary

## If You Don't See "CONFIGURE CONSENT SCREEN" Button

You might already be on a different page. Try:
1. Using the direct link: https://console.cloud.google.com/apis/credentials/consent
2. Or search for "OAuth consent screen" in the search box
3. Or navigate: APIs & Services > OAuth consent screen

## Tell Me

What do you see when you go to the OAuth consent screen page?
- A button to configure?
- Configuration tabs?
- Something else?

