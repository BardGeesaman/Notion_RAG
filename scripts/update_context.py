#!/usr/bin/env python3
"""
Update Context and Archive Script

This script:
1. Updates the master context file (context/MASTER_CONTEXT_FOR_NEW_CHAT.md)
2. Archives any extraneous .md files to md_archive/
3. Keeps only README.md in the root directory

Usage:
    python scripts/update_context.py

Or from root directory:
    python scripts/update_context.py
"""

from __future__ import annotations

import shutil
from datetime import datetime
from pathlib import Path

# Get project root (assuming script is in scripts/ folder)
PROJECT_ROOT = Path(__file__).parent.parent
CONTEXT_DIR = PROJECT_ROOT / "context"
ARCHIVE_DIR = PROJECT_ROOT / "md_archive"


def get_current_timestamp() -> str:
    """Get current timestamp for context updates."""
    return datetime.now().strftime("%Y-%m-%d")


def archive_extraneous_md_files() -> list[str]:
    """
    Move extraneous .md files from root to md_archive/.

    Keeps only README.md in root.

    Returns:
        List of files that were archived
    """
    archived = []
    root_md_files = list(PROJECT_ROOT.glob("*.md"))

    # Files to keep in root
    keep_files = {"README.md"}

    # Ensure archive directory exists
    ARCHIVE_DIR.mkdir(exist_ok=True)

    for md_file in root_md_files:
        if md_file.name in keep_files:
            continue  # Skip files we want to keep

        # Move to archive
        archive_path = ARCHIVE_DIR / md_file.name
        if archive_path.exists():
            # If file already exists in archive, append timestamp
            stem = md_file.stem
            suffix = md_file.suffix
            timestamp = get_current_timestamp()
            archive_path = ARCHIVE_DIR / f"{stem}_{timestamp}{suffix}"

        shutil.move(str(md_file), str(archive_path))
        archived.append(md_file.name)
        print(f"  📦 Archived: {md_file.name}")

    return archived


def update_master_context_timestamp() -> None:
    """Update the timestamp in MASTER_CONTEXT_FOR_NEW_CHAT.md."""
    master_context_path = CONTEXT_DIR / "MASTER_CONTEXT_FOR_NEW_CHAT.md"

    if not master_context_path.exists():
        print(f"  ⚠️  Warning: {master_context_path} not found")
        return

    # Read current content
    content = master_context_path.read_text(encoding="utf-8")

    # Update timestamp
    timestamp = get_current_timestamp()
    old_pattern = "**Last Updated**: [Auto-updated on each significant change]"
    new_pattern = f"**Last Updated**: {timestamp}"

    if old_pattern in content:
        content = content.replace(old_pattern, new_pattern)
        master_context_path.write_text(content, encoding="utf-8")
        print("  ✅ Updated timestamp in MASTER_CONTEXT_FOR_NEW_CHAT.md")
    else:
        # Try alternative pattern
        import re
        pattern = r"\*\*Last Updated\*\*:.*"
        if re.search(pattern, content):
            content = re.sub(pattern, f"**Last Updated**: {timestamp}", content)
            master_context_path.write_text(content, encoding="utf-8")
            print("  ✅ Updated timestamp in MASTER_CONTEXT_FOR_NEW_CHAT.md")
        else:
            print("  ⚠️  Could not find timestamp pattern in MASTER_CONTEXT_FOR_NEW_CHAT.md")


def check_context_files() -> dict[str, bool]:
    """Check if all expected context files exist."""
    expected_files = [
        "MASTER_CONTEXT_FOR_NEW_CHAT.md",
        "NEW_CHAT_QUICK_START.md",
        "HIGH_LEVEL_STRATEGY_CONTEXT.md",
        "UNIFIED_STRATEGIC_ROADMAP.md",
        "README.md",
    ]

    status = {}
    for filename in expected_files:
        path = CONTEXT_DIR / filename
        status[filename] = path.exists()

    return status


def main() -> None:
    """Main function to update context and archive files."""
    print("=" * 70)
    print("🔄 Updating Context and Archiving Files")
    print("=" * 70)
    print()

    # Step 1: Archive extraneous .md files
    print("📦 Step 1: Archiving extraneous .md files...")
    archived = archive_extraneous_md_files()
    if archived:
        print(f"  ✅ Archived {len(archived)} file(s)")
    else:
        print("  ✅ No files to archive (root directory is clean)")
    print()

    # Step 2: Update master context timestamp
    print("📝 Step 2: Updating master context timestamp...")
    update_master_context_timestamp()
    print()

    # Step 3: Verify context files
    print("✅ Step 3: Verifying context files...")
    status = check_context_files()
    all_present = all(status.values())

    if all_present:
        print("  ✅ All expected context files are present")
    else:
        print("  ⚠️  Some context files are missing:")
        for filename, exists in status.items():
            if not exists:
                print(f"     - {filename}")
    print()

    # Summary
    print("=" * 70)
    print("✅ Context Update Complete!")
    print("=" * 70)
    print(f"  • Files archived: {len(archived)}")
    print(f"  • Context files verified: {sum(status.values())}/{len(status)}")
    print()
    print("📋 Next steps:")
    print("  • Review context/MASTER_CONTEXT_FOR_NEW_CHAT.md for any manual updates")
    print("  • Check md_archive/ for archived files")
    print()


if __name__ == "__main__":
    main()

