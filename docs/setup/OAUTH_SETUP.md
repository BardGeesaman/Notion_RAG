# OAuth Setup Guide

Complete guide for setting up OAuth2 authentication for various services (Gmail, Google APIs, etc.).

## Overview

OAuth2 is a secure authentication method that allows applications to access user data without storing passwords. The user authenticates once, and the application receives a token for future use.

## Common OAuth Setup Steps

### 1. Create OAuth Client

1. Go to the service's developer console (e.g., Google Cloud Console)
2. Navigate to APIs & Services > Credentials
3. Click "+ CREATE CREDENTIALS" > "OAuth client ID"
4. Configure the OAuth consent screen if prompted
5. Select application type (usually "Desktop app" for local development)
6. Download the credentials JSON file

### 2. Configure OAuth Consent Screen

**First-time setup required:**

1. Go to APIs & Services > OAuth consent screen
2. Select user type (External or Internal)
3. Fill in app information:
   - App name
   - User support email
   - Developer contact information
4. Add required scopes (permissions)
5. Add test users (for external apps in testing)
6. Save and continue through all steps

### 3. Save Credentials

- Save the downloaded JSON file to `credentials/` directory
- Name it appropriately (e.g., `gmail_credentials.json`, `oauth_credentials.json`)
- **Never commit credentials to git!** (should be in `.gitignore`)

### 4. Authenticate

Run the relevant script which will:
- Open a browser window
- Prompt you to sign in
- Request permissions
- Save the token for future use

## Service-Specific Guides

### Gmail API

See [Gmail Setup Guide](./GMAIL_SETUP.md) for detailed Gmail-specific instructions.

### Other Google APIs

The process is similar:
1. Enable the specific API in Google Cloud Console
2. Create OAuth credentials
3. Configure consent screen with appropriate scopes
4. Authenticate and use

## Troubleshooting

### Common Issues

**403 Error: Access Denied**
- Verify test users are added (for external apps)
- Check that required APIs are enabled
- Ensure scopes are correctly configured

**Token Expired**
- Delete the token file and re-authenticate
- Tokens typically expire after a set period

**Can't Find Credentials**
- Verify file path is correct
- Check file permissions
- Ensure credentials directory exists

## Security Best Practices

1. **Never commit credentials to version control**
2. **Use environment variables for sensitive data**
3. **Rotate credentials periodically**
4. **Use least-privilege scopes** (only request what you need)
5. **Store tokens securely** (encrypted if possible)

## Additional Resources

For specific setup instructions, see:
- `docs/setup/GMAIL_SETUP.md` - Gmail-specific setup
- `docs/setup/CHECK_OAUTH_PAGE.md` - OAuth page verification
- `docs/setup/CONFIGURE_CONSENT_SCREEN_STEPS.md` - Detailed consent screen setup

