#!/usr/bin/env python3
"""
Verify Gmail credentials file is correctly placed and formatted.
"""

import json
import sys
from pathlib import Path

# Add project root to path
project_root = Path(__file__).parent.parent
sys.path.insert(0, str(project_root))

def main():
    print("=" * 60)
    print("Gmail Credentials Verification")
    print("=" * 60)

    creds_path = project_root / "credentials" / "gmail_credentials.json"

    # Check if file exists
    if not creds_path.exists():
        print("\n❌ Credentials file NOT found at:")
        print(f"   {creds_path}")
        print("\nPlease:")
        print("1. Download OAuth2 credentials from Google Cloud Console")
        print("2. Save as: credentials/gmail_credentials.json")
        print("\nSee GMAIL_CREDENTIALS_GUIDE.md for detailed steps")
        sys.exit(1)

    print("\n✅ Credentials file found at:")
    print(f"   {creds_path}")

    # Check if it's valid JSON
    try:
        with open(creds_path, 'r') as f:
            creds = json.load(f)

        print("\n✅ File is valid JSON")

        # Check structure
        if "installed" in creds:
            client_id = creds["installed"].get("client_id", "")
            client_secret = creds["installed"].get("client_secret", "")

            if client_id and client_secret:
                print("\n✅ Credentials structure is correct")
                print(f"   Client ID: {client_id[:30]}...")
                print(f"   Client Secret: {'*' * 20} (hidden)")
                print("\n✅ Credentials file looks good!")
                print("\nNext step: Run authentication")
                print("   python scripts/setup_gmail.py")
                return 0
            else:
                print("\n⚠️  Credentials file missing required fields")
                print("   Expected: installed.client_id and installed.client_secret")
                sys.exit(1)
        else:
            print("\n⚠️  Credentials file has unexpected structure")
            print("   Expected 'installed' key with client_id and client_secret")
            print(f"   File keys: {list(creds.keys())}")
            sys.exit(1)

    except json.JSONDecodeError as e:
        print(f"\n❌ File is not valid JSON: {e}")
        print("   Please check the file contents")
        sys.exit(1)
    except Exception as e:
        print(f"\n❌ Error reading file: {e}")
        sys.exit(1)

if __name__ == "__main__":
    sys.exit(main())

