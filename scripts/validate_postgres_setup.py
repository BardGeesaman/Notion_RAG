#!/usr/bin/env python3
"""
Validate Postgres database setup and connection.

This script checks:
- Postgres connection configuration
- Database connectivity
- Schema tables exist
- Migrations are applied
"""

import sys
from pathlib import Path

# Add parent directory to path
sys.path.insert(0, str(Path(__file__).parent.parent))

from amprenta_rag.config import get_config
from amprenta_rag.database.base import get_engine, get_session_local
from amprenta_rag.database.models import (
    Program,
    Experiment,
    Dataset,
    Feature,
    Signature,
    SignatureComponent,
)
from amprenta_rag.logging_utils import get_logger

logger = get_logger(__name__)


def check_config() -> bool:
    """Check Postgres configuration is present."""
    print("=" * 80)
    print("1. Checking Postgres Configuration...")
    print("=" * 80)

    try:
        cfg = get_config()
        postgres_cfg = cfg.postgres

        print(f"  Host: {postgres_cfg.host}")
        print(f"  Port: {postgres_cfg.port}")
        print(f"  Database: {postgres_cfg.db}")
        print(f"  User: {postgres_cfg.user}")
        print(f"  Echo (SQL logging): {postgres_cfg.echo}")

        if postgres_cfg.url:
            print(f"  Connection URL: {postgres_cfg.url[:50]}...")
        else:
            print("  Connection URL: (constructed from components)")

        if not postgres_cfg.user or not postgres_cfg.db:
            print("  ❌ ERROR: Postgres user or database not configured")
            return False

        print("  ✅ Configuration looks good")
        return True

    except Exception as e:
        print(f"  ❌ ERROR: Failed to load configuration: {e}")
        return False


def check_connection() -> bool:
    """Check database connection."""
    print("\n" + "=" * 80)
    print("2. Testing Database Connection...")
    print("=" * 80)

    try:
        engine = get_engine()

        # Test connection
        from sqlalchemy import text
        with engine.connect() as conn:
            result = conn.execute(text("SELECT version();"))
            version = result.scalar()
            print(f"  ✅ Connected successfully")
            print(f"  PostgreSQL version: {version[:50]}...")

        return True

    except Exception as e:
        print(f"  ❌ ERROR: Failed to connect to database: {e}")
        print("\n  Troubleshooting:")
        print("    - Check Postgres is running")
        print("    - Verify connection settings in .env")
        print("    - Check database exists: createdb amprenta_rag")
        return False


def check_schema() -> bool:
    """Check schema tables exist."""
    print("\n" + "=" * 80)
    print("3. Checking Schema Tables...")
    print("=" * 80)

    try:
        engine = get_engine()

        # Expected tables
        expected_tables = [
            "programs",
            "experiments",
            "datasets",
            "features",
            "signatures",
            "signature_components",
            "dataset_features",  # Association table
            "dataset_programs",  # Association table
            "dataset_experiments",  # Association table
            "experiment_programs",  # Association table
        ]

        # Query information_schema for existing tables
        from sqlalchemy import text
        with engine.connect() as conn:
            result = conn.execute(
                text("""
                    SELECT table_name
                    FROM information_schema.tables
                    WHERE table_schema = 'public'
                    ORDER BY table_name
                """)
            )
            existing_tables = {row[0] for row in result}

        print(f"  Found {len(existing_tables)} tables in database")

        missing_tables = []
        for table in expected_tables:
            if table in existing_tables:
                print(f"  ✅ {table}")
            else:
                print(f"  ❌ {table} (missing)")
                missing_tables.append(table)

        if missing_tables:
            print(f"\n  ❌ ERROR: {len(missing_tables)} tables missing")
            print("\n  Run migrations to create tables:")
            print("    python scripts/migrate_database.py")
            return False

        print(f"\n  ✅ All {len(expected_tables)} expected tables exist")
        return True

    except Exception as e:
        print(f"  ❌ ERROR: Failed to check schema: {e}")
        return False


def check_migrations() -> bool:
    """Check Alembic migrations status."""
    print("\n" + "=" * 80)
    print("4. Checking Migration Status...")
    print("=" * 80)

    try:
        import subprocess

        # Check if alembic_version table exists
        from sqlalchemy import text
        engine = get_engine()
        with engine.connect() as conn:
            result = conn.execute(
                text("""
                    SELECT EXISTS (
                        SELECT FROM information_schema.tables
                        WHERE table_schema = 'public'
                        AND table_name = 'alembic_version'
                    )
                """)
            )
            alembic_table_exists = result.scalar()

        if not alembic_table_exists:
            print("  ⚠️  Alembic version table not found")
            print("     This means migrations haven't been applied yet")
            print("\n     Run migrations:")
            print("       python scripts/migrate_database.py")
            return False

        # Get current revision
        from sqlalchemy import text
        with engine.connect() as conn:
            result = conn.execute(text("SELECT version_num FROM alembic_version;"))
            current_revision = result.scalar()

        print(f"  ✅ Migrations applied")
        print(f"  Current revision: {current_revision}")

        # Check if up to date
        try:
            result = subprocess.run(
                ["alembic", "current"],
                capture_output=True,
                text=True,
                cwd=Path(__file__).parent.parent,
            )
            if result.returncode == 0:
                print(f"  Alembic status: {result.stdout.strip()}")
        except Exception:
            print("  ⚠️  Could not run 'alembic current' (alembic may not be in PATH)")

        return True

    except Exception as e:
        print(f"  ❌ ERROR: Failed to check migrations: {e}")
        return False


def check_models() -> bool:
    """Check SQLAlchemy models can be queried."""
    print("\n" + "=" * 80)
    print("5. Testing SQLAlchemy Models...")
    print("=" * 80)

    try:
        SessionLocal = get_session_local()
        db = SessionLocal()

        try:
            # Try to query each model (even if empty)
            model_checks = [
                ("Programs", Program),
                ("Experiments", Experiment),
                ("Datasets", Dataset),
                ("Features", Feature),
                ("Signatures", Signature),
                ("Signature Components", SignatureComponent),
            ]

            all_passed = True
            for name, model in model_checks:
                try:
                    count = db.query(model).count()
                    print(f"  ✅ {name}: {count} records")
                except Exception as e:
                    print(f"  ❌ {name}: Error - {e}")
                    all_passed = False

            if all_passed:
                print("\n  ✅ All models are queryable")
                return True
            else:
                print("\n  ❌ Some models failed")
                return False

        finally:
            db.close()

    except Exception as e:
        print(f"  ❌ ERROR: Failed to test models: {e}")
        return False


def main() -> None:
    """Run all validation checks."""
    print("\n" + "=" * 80)
    print("POSTGRES DATABASE VALIDATION")
    print("=" * 80)
    print()

    checks = [
        ("Configuration", check_config),
        ("Connection", check_connection),
        ("Schema", check_schema),
        ("Migrations", check_migrations),
        ("Models", check_models),
    ]

    results = []
    for name, check_func in checks:
        try:
            passed = check_func()
            results.append((name, passed))
        except Exception as e:
            logger.exception("Error during %s check", name)
            print(f"\n  ❌ Fatal error: {e}")
            results.append((name, False))

    # Summary
    print("\n" + "=" * 80)
    print("SUMMARY")
    print("=" * 80)

    passed = sum(1 for _, result in results if result)
    total = len(results)

    for name, result in results:
        status = "✅ PASS" if result else "❌ FAIL"
        print(f"  {status}: {name}")

    print(f"\n  Total: {passed}/{total} checks passed")
    print()

    if passed == total:
        print("✅ All checks passed! Postgres is ready to use.")
        sys.exit(0)
    else:
        print("⚠️  Some checks failed. Please review the errors above.")
        sys.exit(1)


if __name__ == "__main__":
    main()

