#!/usr/bin/env python3
"""
Comprehensive validation script for TIER 3 infrastructure.

This script runs all validation checks:
- Postgres setup
- API endpoints
- Integration tests

Usage:
    python scripts/validate_tier3_complete.py
    python scripts/validate_tier3_complete.py --skip-api  # Skip API tests (if server not running)
"""

import argparse
import subprocess
import sys
from pathlib import Path

# Add parent directory to path
sys.path.insert(0, str(Path(__file__).parent.parent))


def run_validation_script(script_name: str, description: str) -> bool:
    """Run a validation script and return success status."""
    print("\n" + "=" * 80)
    print(f"Running: {description}")
    print("=" * 80)
    print()

    script_path = Path(__file__).parent / script_name

    if not script_path.exists():
        print(f"  ⚠️  Script not found: {script_path}")
        return False

    try:
        result = subprocess.run(
            [sys.executable, str(script_path)],
            cwd=Path(__file__).parent.parent,
            timeout=60,
        )
        return result.returncode == 0
    except subprocess.TimeoutExpired:
        print(f"  ❌ Script timed out after 60 seconds")
        return False
    except Exception as e:
        print(f"  ❌ Error running script: {e}")
        return False


def run_pytest_tests(test_path: str, description: str) -> bool:
    """Run pytest tests and return success status."""
    print("\n" + "=" * 80)
    print(f"Running: {description}")
    print("=" * 80)
    print()

    test_dir = Path(__file__).parent.parent / test_path

    if not test_dir.exists():
        print(f"  ⚠️  Test directory not found: {test_dir}")
        return False

    try:
        result = subprocess.run(
            [sys.executable, "-m", "pytest", str(test_dir), "-v"],
            cwd=Path(__file__).parent.parent,
            timeout=120,
        )
        return result.returncode == 0
    except subprocess.TimeoutExpired:
        print(f"  ❌ Tests timed out after 120 seconds")
        return False
    except FileNotFoundError:
        print(f"  ⚠️  pytest not installed. Install with: pip install pytest")
        return False
    except Exception as e:
        print(f"  ❌ Error running tests: {e}")
        return False


def main() -> None:
    """Run all TIER 3 validation checks."""
    parser = argparse.ArgumentParser(
        description="Comprehensive TIER 3 infrastructure validation"
    )
    parser.add_argument(
        "--skip-api",
        action="store_true",
        help="Skip API endpoint tests (useful if server is not running)",
    )
    parser.add_argument(
        "--skip-tests",
        action="store_true",
        help="Skip pytest unit/integration tests",
    )

    args = parser.parse_args()

    print("\n" + "=" * 80)
    print("TIER 3 COMPREHENSIVE VALIDATION")
    print("=" * 80)
    print()

    results = []

    # 1. Postgres validation
    results.append(
        (
            "Postgres Setup",
            run_validation_script("validate_postgres_setup.py", "Postgres Database Validation"),
        )
    )

    # 2. Database tests
    if not args.skip_tests:
        results.append(
            (
                "Database Tests",
                run_pytest_tests(
                    "amprenta_rag/tests/database",
                    "Database Connection & Model Tests",
                ),
            )
        )

    # 3. API endpoint validation
    if not args.skip_api:
        results.append(
            (
                "API Endpoints",
                run_validation_script("test_api_endpoints.py", "FastAPI Endpoint Validation"),
            )
        )
    else:
        print("\n" + "=" * 80)
        print("Skipping API Endpoint Tests (--skip-api flag)")
        print("=" * 80)
        results.append(("API Endpoints", None))  # Skipped

    # 4. API tests
    if not args.skip_tests and not args.skip_api:
        results.append(
            (
                "API Tests",
                run_pytest_tests("amprenta_rag/tests/api", "FastAPI Endpoint Tests"),
            )
        )

    # 5. Integration tests
    if not args.skip_tests:
        results.append(
            (
                "Integration Tests",
                run_pytest_tests(
                    "amprenta_rag/tests/integration",
                    "Postgres + API Integration Tests",
                ),
            )
        )

    # Summary
    print("\n" + "=" * 80)
    print("VALIDATION SUMMARY")
    print("=" * 80)
    print()

    passed = sum(1 for _, result in results if result is True)
    failed = sum(1 for _, result in results if result is False)
    skipped = sum(1 for _, result in results if result is None)

    for name, result in results:
        if result is None:
            status = "⏭️  SKIPPED"
        elif result:
            status = "✅ PASS"
        else:
            status = "❌ FAIL"
        print(f"  {status}: {name}")

    print()
    print(f"  Total: {passed} passed, {failed} failed, {skipped} skipped")
    print()

    if failed == 0 and passed > 0:
        print("✅ All validations passed! TIER 3 infrastructure is ready.")
        sys.exit(0)
    elif failed > 0:
        print("⚠️  Some validations failed. Please review the errors above.")
        sys.exit(1)
    else:
        print("⚠️  No validations were run.")
        sys.exit(1)


if __name__ == "__main__":
    main()

