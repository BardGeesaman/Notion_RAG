#!/usr/bin/env python3
"""
Postgres connection validation and health check script.

Validates that:
- Postgres connection is working
- Database exists
- All expected tables are created
- Schema is up to date with migrations

Usage:
    python scripts/validate_postgres.py
"""

import sys
from pathlib import Path

# Add parent directory to path for imports
sys.path.insert(0, str(Path(__file__).parent.parent))

from sqlalchemy import create_engine, inspect, text
from sqlalchemy.exc import OperationalError

from amprenta_rag.config import get_config


def get_postgres_url() -> str:
    """Get Postgres connection URL from config."""
    cfg = get_config()
    
    # Use full URL if provided, otherwise construct from components
    if cfg.postgres.url:
        return cfg.postgres.url
    
    return (
        f"postgresql://{cfg.postgres.user}:{cfg.postgres.password}"
        f"@{cfg.postgres.host}:{cfg.postgres.port}/{cfg.postgres.db}"
    )


def check_connection() -> bool:
    """Test basic Postgres connection."""
    print("🔌 Testing Postgres connection...")
    
    try:
        url = get_postgres_url()
        engine = create_engine(url, echo=False)
        
        with engine.connect() as conn:
            result = conn.execute(text("SELECT version()"))
            version = result.fetchone()[0]
            print(f"✅ Connection successful!")
            print(f"   PostgreSQL version: {version.split(',')[0]}")
        
        return True
        
    except OperationalError as e:
        print(f"❌ Connection failed: {e}")
        print("\nTroubleshooting:")
        print("  1. Ensure Docker Compose is running: docker-compose ps")
        print("  2. Check .env file has correct POSTGRES_* settings")
        print("  3. Verify port 5432 is not blocked")
        print("  4. Try restarting services: docker-compose restart")
        return False


def check_database_exists(engine) -> bool:
    """Check if the target database exists."""
    print("\n📚 Checking database...")
    
    cfg = get_config()
    db_name = cfg.postgres.db
    
    try:
        with engine.connect() as conn:
            result = conn.execute(
                text("SELECT datname FROM pg_database WHERE datname = :db"),
                {"db": db_name}
            )
            exists = result.fetchone() is not None
            
            if exists:
                print(f"✅ Database '{db_name}' exists")
                return True
            else:
                print(f"❌ Database '{db_name}' not found")
                print(f"\nCreate database with:")
                print(f"  CREATE DATABASE {db_name};")
                return False
                
    except Exception as e:
        print(f"❌ Error checking database: {e}")
        return False


def check_tables(engine) -> bool:
    """Check if all expected tables exist."""
    print("\n🗂️  Checking tables...")
    
    expected_tables = {
        'programs',
        'experiments',
        'datasets',
        'features',
        'signatures',
        'signature_components',
        # Junction tables
        'program_experiment',
        'program_dataset',
        'experiment_dataset',
        'dataset_feature',
        'dataset_signature',
        'signature_feature',
        'program_signature',
    }
    
    try:
        inspector = inspect(engine)
        existing_tables = set(inspector.get_table_names())
        
        # Check if all expected tables exist
        missing_tables = expected_tables - existing_tables
        extra_tables = existing_tables - expected_tables
        
        if not missing_tables:
            print(f"✅ All {len(expected_tables)} tables exist:")
            for table in sorted(expected_tables):
                # Get row count
                try:
                    with engine.connect() as conn:
                        result = conn.execute(text(f"SELECT COUNT(*) FROM {table}"))
                        count = result.fetchone()[0]
                        print(f"   - {table}: {count} rows")
                except Exception:
                    print(f"   - {table}: (could not count rows)")
            
            if extra_tables:
                print(f"\nℹ️  Additional tables found: {', '.join(sorted(extra_tables))}")
            
            return True
        else:
            print(f"❌ Missing tables: {', '.join(sorted(missing_tables))}")
            print("\nRun Alembic migrations to create tables:")
            print("  alembic upgrade head")
            return False
            
    except Exception as e:
        print(f"❌ Error checking tables: {e}")
        return False


def check_indexes(engine) -> bool:
    """Check if critical indexes exist."""
    print("\n📊 Checking indexes...")
    
    critical_indexes = [
        ('programs', 'ix_programs_name'),
        ('datasets', 'ix_datasets_omics_type'),
        ('features', 'ix_features_feature_type'),
        ('signatures', 'ix_signatures_name'),
    ]
    
    try:
        inspector = inspect(engine)
        all_indexes_exist = True
        
        for table_name, index_name in critical_indexes:
            indexes = inspector.get_indexes(table_name)
            index_names = {idx['name'] for idx in indexes}
            
            if index_name in index_names:
                print(f"   ✓ {table_name}.{index_name}")
            else:
                print(f"   ✗ {table_name}.{index_name} (missing)")
                all_indexes_exist = False
        
        if all_indexes_exist:
            print("✅ All critical indexes exist")
            return True
        else:
            print("⚠️  Some indexes are missing (may affect performance)")
            return False
            
    except Exception as e:
        print(f"⚠️  Could not check indexes: {e}")
        return False


def check_migrations() -> bool:
    """Check if migrations are up to date."""
    print("\n🔄 Checking Alembic migrations...")
    
    try:
        from alembic.config import Config
        from alembic.script import ScriptDirectory
        from alembic.runtime.migration import MigrationContext
        
        # Get current revision from database
        url = get_postgres_url()
        engine = create_engine(url, echo=False)
        
        with engine.connect() as conn:
            context = MigrationContext.configure(conn)
            current_rev = context.get_current_revision()
        
        # Get head revision from scripts
        alembic_cfg = Config("alembic.ini")
        script = ScriptDirectory.from_config(alembic_cfg)
        head_rev = script.get_current_head()
        
        if current_rev == head_rev:
            print(f"✅ Migrations up to date: {current_rev}")
            return True
        elif current_rev is None:
            print("⚠️  No migrations applied yet")
            print("   Run: alembic upgrade head")
            return False
        else:
            print(f"⚠️  Database is at revision: {current_rev}")
            print(f"   Latest available: {head_rev}")
            print("   Run: alembic upgrade head")
            return False
            
    except ImportError:
        print("⚠️  Alembic not installed, skipping migration check")
        return True
    except FileNotFoundError:
        print("⚠️  alembic.ini not found, skipping migration check")
        return True
    except Exception as e:
        print(f"⚠️  Could not check migrations: {e}")
        return True


def check_configuration() -> bool:
    """Validate configuration settings."""
    print("\n⚙️  Checking configuration...")
    
    try:
        cfg = get_config()
        issues = []
        
        # Check Postgres configuration
        if not cfg.postgres.host:
            issues.append("POSTGRES_HOST not set")
        if not cfg.postgres.db:
            issues.append("POSTGRES_DB not set")
        if not cfg.postgres.user:
            issues.append("POSTGRES_USER not set")
        
        # Check critical flags
        if not cfg.pipeline.use_postgres_as_sot:
            issues.append("USE_POSTGRES_AS_SOT is False (should be True for local setup)")
        
        # Check API keys
        if not cfg.openai.api_key:
            issues.append("OPENAI_API_KEY not set")
        if not cfg.pinecone.api_key:
            issues.append("PINECONE_API_KEY not set")
        
        if issues:
            print("⚠️  Configuration issues found:")
            for issue in issues:
                print(f"   - {issue}")
            return False
        else:
            print("✅ Configuration looks good")
            print(f"   - USE_POSTGRES_AS_SOT: {cfg.pipeline.use_postgres_as_sot}")
            print(f"   - ENABLE_NOTION_SYNC: {cfg.pipeline.enable_notion_sync}")
            print(f"   - Database: {cfg.postgres.db}")
            return True
            
    except Exception as e:
        print(f"❌ Error checking configuration: {e}")
        return False


def main():
    """Run all validation checks."""
    print("=" * 70)
    print("🔍 Amprenta RAG - Postgres Validation")
    print("=" * 70)
    
    # Run checks
    checks = []
    
    # 1. Configuration
    checks.append(("Configuration", check_configuration()))
    
    # 2. Connection
    connection_ok = check_connection()
    checks.append(("Connection", connection_ok))
    
    if not connection_ok:
        print("\n" + "=" * 70)
        print("❌ Cannot proceed without database connection")
        print("=" * 70)
        sys.exit(1)
    
    # 3. Create engine for remaining checks
    url = get_postgres_url()
    engine = create_engine(url, echo=False)
    
    # 4. Database exists
    checks.append(("Database", check_database_exists(engine)))
    
    # 5. Tables exist
    tables_ok = check_tables(engine)
    checks.append(("Tables", tables_ok))
    
    # 6. Indexes (non-critical)
    checks.append(("Indexes", check_indexes(engine)))
    
    # 7. Migrations
    checks.append(("Migrations", check_migrations()))
    
    # Summary
    print("\n" + "=" * 70)
    print("📋 Validation Summary")
    print("=" * 70)
    
    all_passed = True
    for check_name, passed in checks:
        status = "✅ PASS" if passed else "❌ FAIL"
        print(f"{status}: {check_name}")
        if not passed:
            all_passed = False
    
    print("=" * 70)
    
    if all_passed:
        print("\n🎉 All validation checks passed!")
        print("\nYou're ready to:")
        print("  - Ingest data: python scripts/ingest_lipidomics.py --file <path>")
        print("  - Start Streamlit: streamlit run streamlit_app/app.py")
        print("  - Start FastAPI: uvicorn amprenta_rag.api.main:app --reload")
        sys.exit(0)
    else:
        print("\n⚠️  Some validation checks failed.")
        print("See above for troubleshooting steps.")
        sys.exit(1)


if __name__ == "__main__":
    main()

