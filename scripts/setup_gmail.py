#!/usr/bin/env python3
"""
Interactive Gmail setup script.

Helps set up Gmail API credentials and test the connection.
"""

import sys
from pathlib import Path

# Add project root to path
project_root = Path(__file__).parent.parent
sys.path.insert(0, str(project_root))

def check_dependencies():
    """Check if Gmail API dependencies are installed."""
    try:
        print("✅ Gmail API dependencies are installed")
        return True
    except ImportError as e:
        print(f"❌ Missing Gmail API dependencies: {e}")
        print("\nInstall with:")
        print("  pip install google-auth google-auth-oauthlib google-auth-httplib2 google-api-python-client")
        print("\nOr:")
        print("  pip install -r requirements_gmail.txt")
        return False

def check_credentials():
    """Check if credentials file exists."""
    creds_path = project_root / "credentials" / "gmail_credentials.json"
    token_path = project_root / "credentials" / "gmail_token.json"

    has_creds = creds_path.exists()
    has_token = token_path.exists()

    print("\n" + "=" * 60)
    print("Gmail API Setup Status")
    print("=" * 60)

    if has_creds:
        print("✅ Credentials file found: credentials/gmail_credentials.json")
    else:
        print("❌ Credentials file NOT found: credentials/gmail_credentials.json")
        print("\n   You need to:")
        print("   1. Go to https://console.cloud.google.com/")
        print("   2. Create/select a project")
        print("   3. Enable Gmail API")
        print("   4. Create OAuth 2.0 Desktop app credentials")
        print("   5. Download and save as: credentials/gmail_credentials.json")

    if has_token:
        print("✅ Token file found: credentials/gmail_token.json")
        print("   (Already authenticated)")
    else:
        print("❌ Token file NOT found: credentials/gmail_token.json")
        print("   (Need to authenticate on first run)")

    print("=" * 60)

    return has_creds, has_token

def test_connection():
    """Test Gmail API connection."""
    print("\n" + "=" * 60)
    print("Testing Gmail API Connection")
    print("=" * 60)

    try:
        from amprenta_rag.clients.gmail_client import GmailClient

        print("\nInitializing Gmail client...")
        client = GmailClient()

        print("Fetching inbox emails (limit 5)...")
        emails = client.fetch_emails(query="in:inbox", max_results=5)

        print(f"✅ Success! Found {len(emails)} email(s)")

        if emails:
            print("\nSample emails:")
            for i, email in enumerate(emails[:3], 1):
                subject = email.get("subject", "(no subject)")[:50]
                from_sender = email.get("from", "unknown")[:30]
                print(f"  {i}. {subject} (from: {from_sender})")

        return True

    except FileNotFoundError as e:
        print(f"\n❌ Error: {e}")
        print("\nPlease download OAuth2 credentials from Google Cloud Console")
        return False
    except Exception as e:
        print(f"\n❌ Connection test failed: {e}")
        print("\nTroubleshooting:")
        print("  - Check credentials file exists")
        print("  - Verify Gmail API is enabled in Google Cloud Console")
        print("  - Try re-authenticating (delete credentials/gmail_token.json)")
        return False

def main():
    print("Gmail API Setup and Test")
    print("=" * 60)

    # Check dependencies
    if not check_dependencies():
        sys.exit(1)

    # Check credentials
    has_creds, has_token = check_credentials()

    if not has_creds:
        print("\n⚠️  Please set up credentials first (see instructions above)")
        sys.exit(1)

    # Test connection
    if not has_token:
        print("\n🔐 First-time authentication will happen during connection test...")
        input("Press Enter to continue (will open browser for OAuth)...")

    success = test_connection()

    if success:
        print("\n" + "=" * 60)
        print("✅ Setup Complete!")
        print("=" * 60)
        print("\nYou can now ingest emails with:")
        print("  python scripts/ingest_gmail.py")
        print("\nOptions:")
        print("  --query 'is:unread'     # Unread emails")
        print("  --days 7                # Last 7 days")
        print("  --dry-run               # Test without ingesting")
    else:
        print("\n⚠️  Setup incomplete. Please check the errors above.")
        sys.exit(1)

if __name__ == "__main__":
    main()

