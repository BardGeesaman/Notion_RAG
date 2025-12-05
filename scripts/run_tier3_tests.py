#!/usr/bin/env python3
"""
Test runner for TIER 3 tests with dependency checking.

This script:
- Checks for required dependencies
- Runs tests that can run
- Skips tests that require missing dependencies
- Provides clear feedback
"""

import sys
from pathlib import Path

# Add parent directory to path
sys.path.insert(0, str(Path(__file__).parent.parent))


def check_dependency(name: str, module_name: str = None) -> tuple[bool, str]:
    """Check if a dependency is available."""
    if module_name is None:
        module_name = name.lower()
    
    try:
        __import__(module_name)
        return True, f"✅ {name} installed"
    except ImportError:
        return False, f"❌ {name} not installed"


def main() -> None:
    """Run TIER 3 tests with dependency checking."""
    print("=" * 80)
    print("TIER 3 TEST RUNNER")
    print("=" * 80)
    print()
    
    # Check dependencies
    print("Checking dependencies...")
    print("-" * 80)
    
    dependencies = [
        ("pytest", "pytest"),
        ("SQLAlchemy", "sqlalchemy"),
        ("FastAPI", "fastapi"),
        ("psycopg2", "psycopg2"),
        ("pydantic", "pydantic"),
    ]
    
    available = {}
    for name, module in dependencies:
        is_available, message = check_dependency(name, module)
        print(f"  {message}")
        available[name] = is_available
    
    print()
    
    # Determine what we can test
    can_test_database = available.get("SQLAlchemy", False)
    can_test_api = available.get("FastAPI", False)
    can_test_postgres = available.get("psycopg2", False)
    
    print("=" * 80)
    print("Test Plan")
    print("=" * 80)
    print()
    
    if can_test_database:
        print("  ✅ Can test database models (SQLAlchemy available)")
    else:
        print("  ❌ Cannot test database models (SQLAlchemy missing)")
    
    if can_test_api:
        print("  ✅ Can test API endpoints (FastAPI available)")
    else:
        print("  ❌ Cannot test API endpoints (FastAPI missing)")
        print("     Install: pip install fastapi uvicorn")
    
    if can_test_postgres:
        print("  ✅ Can test Postgres connection (psycopg2 available)")
    else:
        print("  ❌ Cannot test Postgres connection (psycopg2 missing)")
        print("     Install: pip install psycopg2-binary")
    
    print()
    
    # Run tests
    import subprocess
    
    tests_to_run = []
    
    if can_test_api:
        tests_to_run.append("amprenta_rag/tests/api/")
    
    if can_test_database and can_test_postgres:
        tests_to_run.append("amprenta_rag/tests/database/")
        tests_to_run.append("amprenta_rag/tests/integration/")
    
    if not tests_to_run:
        print("=" * 80)
        print("❌ No tests can run - missing dependencies")
        print("=" * 80)
        print()
        print("Install missing dependencies:")
        print("  pip install -r requirements.txt")
        sys.exit(1)
    
    print("=" * 80)
    print("Running Tests")
    print("=" * 80)
    print()
    
    for test_path in tests_to_run:
        print(f"Running: {test_path}")
        print("-" * 80)
        
        try:
            result = subprocess.run(
                [sys.executable, "-m", "pytest", test_path, "-v"],
                cwd=Path(__file__).parent.parent,
                timeout=120,
            )
            
            if result.returncode == 0:
                print(f"  ✅ {test_path} passed")
            else:
                print(f"  ❌ {test_path} failed")
        except subprocess.TimeoutExpired:
            print(f"  ⏱️  {test_path} timed out")
        except Exception as e:
            print(f"  ❌ Error running {test_path}: {e}")
        
        print()
    
    print("=" * 80)
    print("Test Run Complete")
    print("=" * 80)


if __name__ == "__main__":
    main()

