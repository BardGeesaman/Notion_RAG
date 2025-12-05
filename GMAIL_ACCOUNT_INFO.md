# Gmail Account Setup - Which Account to Use?

## Short Answer

✅ **Use your MAIN Google account** to log into Google Cloud Console
✅ **Use `amprenta.email@gmail.com`** later when the script authenticates

## Detailed Explanation

### Google Cloud Console Login
- **Use any Google account** - Your personal/main account is fine
- This is just for setting up the API credentials
- The account you use here doesn't matter for reading emails

### Email Account Access
- **`amprenta.email@gmail.com`** is the account whose emails you want to read
- You'll sign in with this account **later** when running the script
- The script will ask for permission to read emails from this account

### How It Works

1. **Now (Setup)**: 
   - Log into Google Cloud Console with **your main account**
   - Create project, enable API, create credentials
   - Download credentials file

2. **Later (First Run)**:
   - Run: `python scripts/ingest_gmail.py`
   - Browser opens for authentication
   - Sign in with **`amprenta.email@gmail.com`**
   - Grant permission to read emails
   - Token saved for future use

3. **OAuth Consent Screen**:
   - When configuring OAuth, add `amprenta.email@gmail.com` as a test user
   - This allows that account to use your OAuth app

## Recommendation

✅ **Use your main Google account** for Google Cloud Console
- Easier access
- Likely already has permissions
- Can manage projects better

✅ **Use `amprenta.email@gmail.com`** for email access
- This is the account whose emails you want to read
- You'll authenticate with this account when the script runs

## Summary

| Step | Account to Use |
|------|---------------|
| Google Cloud Console login | Your main/personal account ✅ |
| Create project/credentials | Your main account ✅ |
| OAuth test users | Add `amprenta.email@gmail.com` ✅ |
| Script authentication (later) | `amprenta.email@gmail.com` ✅ |

