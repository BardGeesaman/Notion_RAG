#!/usr/bin/env python3
"""
Helper script to create the initial Alembic migration.

This script creates the initial migration based on the SQLAlchemy models.
"""

import subprocess
import sys
from pathlib import Path

# Add parent directory to path
sys.path.insert(0, str(Path(__file__).parent.parent))


def main() -> None:
    """Create the initial migration."""
    print("=" * 80)
    print("Creating initial Alembic migration...")
    print("=" * 80)
    print()
    
    # Run alembic revision with autogenerate
    cmd = [
        sys.executable,
        "-m",
        "alembic",
        "revision",
        "--autogenerate",
        "-m",
        "Initial schema: programs, experiments, datasets, features, signatures",
    ]
    
    print(f"Running: {' '.join(cmd)}")
    print()
    
    try:
        result = subprocess.run(cmd, check=True, cwd=Path(__file__).parent.parent)
        print()
        print("✅ Initial migration created successfully!")
        print()
        print("Next steps:")
        print("  1. Review the migration file in alembic/versions/")
        print("  2. Apply the migration: alembic upgrade head")
        print()
    except subprocess.CalledProcessError as e:
        print(f"❌ Error creating migration: {e}")
        sys.exit(1)
    except FileNotFoundError:
        print("❌ Error: Alembic not found. Make sure it's installed:")
        print("  pip install alembic")
        sys.exit(1)


if __name__ == "__main__":
    main()

