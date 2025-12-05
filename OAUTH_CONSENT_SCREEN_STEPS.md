# OAuth Consent Screen Configuration Steps

You're on the OAuth overview page. Here's what to do:

## Step-by-Step Configuration

### 1. Select User Type
- You should see "User Type" section
- Select: **"External"** 
  - (Choose "Internal" only if you're using Google Workspace/Enterprise)
- Click **"CREATE"** button

### 2. App Information Tab

Fill in these fields:

- **App name**: `Amprenta Email Ingestor`
- **User support email**: 
  - Select your email from the dropdown, OR
  - Enter: `amprenta.email@gmail.com`
- **App logo**: Leave blank (optional)
- **Application home page**: Leave blank or enter a URL (optional)
- **Application privacy policy link**: Leave blank (optional)
- **Application terms of service link**: Leave blank (optional)
- **Authorized domains**: Leave blank (optional)
- **Developer contact information**: Enter your email address

Then click **"SAVE AND CONTINUE"** at the bottom

### 3. Scopes Tab

- Click **"+ ADD OR REMOVE SCOPES"** button
- In the filter/search box, type: `gmail.readonly`
- Find and check the box for: **`https://www.googleapis.com/auth/gmail.readonly`**
- Click **"UPDATE"** button at the bottom
- Click **"SAVE AND CONTINUE"** at the bottom

### 4. Test Users Tab

- Click **"+ ADD USERS"** button
- Enter: `amprenta.email@gmail.com`
- Click **"ADD"** button
- The email should appear in the list
- Click **"SAVE AND CONTINUE"** at the bottom

### 5. Summary Tab

- Review the information you entered
- Click **"BACK TO DASHBOARD"** button

## What You Should See

After completing all tabs, you'll be back at the OAuth consent screen dashboard, and you should see:
- Your app name: "Amprenta Email Ingestor"
- Status: "Testing" (that's fine - it's in test mode)

## Next Step

Once you've completed the OAuth consent screen configuration, the next step is to create the OAuth 2.0 credentials!

