#!/usr/bin/env python3
"""
Helper script to apply Alembic migrations.

This script applies pending migrations to bring the database up to date.
"""

import subprocess
import sys
from pathlib import Path

# Add parent directory to path
sys.path.insert(0, str(Path(__file__).parent.parent))


def main() -> None:
    """Apply migrations."""
    import argparse
    
    parser = argparse.ArgumentParser(description="Apply Alembic migrations")
    parser.add_argument(
        "--revision",
        default="head",
        help="Target revision (default: head)",
    )
    parser.add_argument(
        "--dry-run",
        action="store_true",
        help="Show what would be done without applying",
    )
    
    args = parser.parse_args()
    
    print("=" * 80)
    print("Applying Alembic migrations...")
    print("=" * 80)
    print()
    
    if args.dry_run:
        print("DRY RUN MODE - No changes will be made")
        print()
        cmd = [sys.executable, "-m", "alembic", "history"]
    else:
        cmd = [sys.executable, "-m", "alembic", "upgrade", args.revision]
    
    print(f"Running: {' '.join(cmd)}")
    print()
    
    try:
        result = subprocess.run(cmd, check=True, cwd=Path(__file__).parent.parent)
        if not args.dry_run:
            print()
            print(f"✅ Migrations applied successfully to revision: {args.revision}")
        sys.exit(0)
    except subprocess.CalledProcessError as e:
        print(f"❌ Error applying migrations: {e}")
        sys.exit(1)
    except FileNotFoundError:
        print("❌ Error: Alembic not found. Make sure it's installed:")
        print("  pip install alembic")
        sys.exit(1)


if __name__ == "__main__":
    main()

