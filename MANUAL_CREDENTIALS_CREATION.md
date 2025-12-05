# Manual Credentials File Creation

If the download button didn't work, you can create the file manually!

## Step 1: Copy the Credentials from Popup

In the popup window that shows your OAuth client, you should see:
- **Client ID**: A long string (looks like: `123456789-abcdefg.apps.googleusercontent.com`)
- **Client secret**: Another long string (looks like: `GOCSPX-abcdefghijklmnopqrstuvwxyz`)

**Copy both of these** - you'll need them!

## Step 2: Create the JSON File

I'll create a template file for you. Tell me:
1. What's your Client ID?
2. What's your Client secret?

Or, if you want to create it yourself:

Create a file at: `/Users/bard/Documents/Notion RAG/credentials/gmail_credentials.json`

With this content (replace YOUR_CLIENT_ID and YOUR_CLIENT_SECRET):

```json
{
  "installed": {
    "client_id": "YOUR_CLIENT_ID",
    "project_id": "your-project-id",
    "auth_uri": "https://accounts.google.com/o/oauth2/auth",
    "token_uri": "https://oauth2.googleapis.com/token",
    "auth_provider_x509_cert_url": "https://www.googleapis.com/oauth2/v1/certs",
    "client_secret": "YOUR_CLIENT_SECRET",
    "redirect_uris": ["http://localhost"]
  }
}
```

## Alternative: Right-Click Download

If there's a download icon, try:
- Right-click on it → "Save link as..."
- Or look for a "Save" option in the popup

## Check Browser Downloads

1. Check your browser's download bar (usually at bottom)
2. Or press Cmd+Shift+J (Chrome) or Cmd+Option+L (Firefox) to see downloads
3. Look for any .json file

## Tell Me

Can you see the Client ID and Client secret in the popup? 
- If yes, tell me what they are and I'll create the file for you!
- Or try right-clicking the download button to save

