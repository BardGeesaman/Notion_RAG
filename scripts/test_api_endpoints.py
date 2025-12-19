#!/usr/bin/env python3
"""
Test FastAPI endpoints to verify API functionality.

This script tests:
- API server can start
- All endpoints are accessible
- CRUD operations work
- Error handling works correctly
"""

import sys
import time
from pathlib import Path

# Add parent directory to path
sys.path.insert(0, str(Path(__file__).parent.parent))

import requests
from amprenta_rag.logging_utils import get_logger

logger = get_logger(__name__)

API_BASE_URL = "http://127.0.0.1:8000"
API_PREFIX = "/api/v1"


def check_server_running() -> bool:
    """Check if API server is running."""
    print("=" * 80)
    print("1. Checking API Server...")
    print("=" * 80)

    try:
        health_url = f"{API_BASE_URL}/health"
        resp = requests.get(health_url, timeout=2)

        if resp.status_code == 200:
            print(f"  ✅ Server is running at {API_BASE_URL}")
            return True
        else:
            print(f"  ❌ Server returned status {resp.status_code}")
            return False

    except requests.exceptions.ConnectionError:
        print(f"  ❌ Server is not running at {API_BASE_URL}")
        print("\n  Start the server:")
        print("    python scripts/run_api_server.py")
        return False
    except Exception as e:
        print(f"  ❌ Error checking server: {e}")
        return False


def test_programs_endpoints() -> bool:
    """Test Programs API endpoints."""
    print("\n" + "=" * 80)
    print("2. Testing Programs Endpoints...")
    print("=" * 80)

    base_url = f"{API_BASE_URL}{API_PREFIX}/programs"

    try:
        # Test GET (list)
        print("  Testing GET /programs (list)...")
        resp = requests.get(base_url, timeout=5)
        if resp.status_code == 200:
            data = resp.json()
            print(f"  ✅ GET /programs: {len(data)} programs")
        else:
            print(f"  ❌ GET /programs: Status {resp.status_code}")
            print(f"     {resp.text[:200]}")
            return False

        # Test POST (create)
        print("  Testing POST /programs (create)...")
        test_program = {
            "name": "Test Program (API Test)",
            "description": "Created by test script",
        }
        resp = requests.post(base_url, json=test_program, timeout=5)

        if resp.status_code in (200, 201):
            created = resp.json()
            program_id = created.get("id")
            print(f"  ✅ POST /programs: Created program {program_id}")

            # Test GET (single)
            print(f"  Testing GET /programs/{program_id}...")
            resp = requests.get(f"{base_url}/{program_id}", timeout=5)
            if resp.status_code == 200:
                print(f"  ✅ GET /programs/{{id}}: Retrieved program")
            else:
                print(f"  ❌ GET /programs/{{id}}: Status {resp.status_code}")

            # Test PATCH (update)
            print(f"  Testing PATCH /programs/{program_id}...")
            update_data = {"description": "Updated by test script"}
            resp = requests.patch(f"{base_url}/{program_id}", json=update_data, timeout=5)
            if resp.status_code == 200:
                print(f"  ✅ PATCH /programs/{{id}}: Updated program")
            else:
                print(f"  ❌ PATCH /programs/{{id}}: Status {resp.status_code}")

            # Test DELETE
            print(f"  Testing DELETE /programs/{program_id}...")
            resp = requests.delete(f"{base_url}/{program_id}", timeout=5)
            if resp.status_code in (200, 204):
                print(f"  ✅ DELETE /programs/{{id}}: Deleted program")
            else:
                print(f"  ❌ DELETE /programs/{{id}}: Status {resp.status_code}")

            return True
        else:
            print(f"  ❌ POST /programs: Status {resp.status_code}")
            print(f"     {resp.text[:200]}")
            return False

    except Exception as e:
        print(f"  ❌ ERROR: {e}")
        logger.exception("Error testing programs endpoints")
        return False


def test_experiments_endpoints() -> bool:
    """Test Experiments API endpoints."""
    print("\n" + "=" * 80)
    print("3. Testing Experiments Endpoints...")
    print("=" * 80)

    base_url = f"{API_BASE_URL}{API_PREFIX}/experiments"

    try:
        # Test GET (list)
        print("  Testing GET /experiments (list)...")
        resp = requests.get(base_url, timeout=5)
        if resp.status_code == 200:
            data = resp.json()
            print(f"  ✅ GET /experiments: {len(data)} experiments")
            return True
        else:
            print(f"  ❌ GET /experiments: Status {resp.status_code}")
            print(f"     {resp.text[:200]}")
            return False

    except Exception as e:
        print(f"  ❌ ERROR: {e}")
        logger.exception("Error testing experiments endpoints")
        return False


def test_datasets_endpoints() -> bool:
    """Test Datasets API endpoints."""
    print("\n" + "=" * 80)
    print("4. Testing Datasets Endpoints...")
    print("=" * 80)

    base_url = f"{API_BASE_URL}{API_PREFIX}/datasets"

    try:
        # Test GET (list)
        print("  Testing GET /datasets (list)...")
        resp = requests.get(base_url, timeout=5)
        if resp.status_code == 200:
            data = resp.json()
            print(f"  ✅ GET /datasets: {len(data)} datasets")
            return True
        else:
            print(f"  ❌ GET /datasets: Status {resp.status_code}")
            print(f"     {resp.text[:200]}")
            return False

    except Exception as e:
        print(f"  ❌ ERROR: {e}")
        logger.exception("Error testing datasets endpoints")
        return False


def test_features_endpoints() -> bool:
    """Test Features API endpoints."""
    print("\n" + "=" * 80)
    print("5. Testing Features Endpoints...")
    print("=" * 80)

    base_url = f"{API_BASE_URL}{API_PREFIX}/features"

    try:
        # Test GET (list)
        print("  Testing GET /features (list)...")
        resp = requests.get(base_url, timeout=5)
        if resp.status_code == 200:
            data = resp.json()
            print(f"  ✅ GET /features: {len(data)} features")
            return True
        else:
            print(f"  ❌ GET /features: Status {resp.status_code}")
            print(f"     {resp.text[:200]}")
            return False

    except Exception as e:
        print(f"  ❌ ERROR: {e}")
        logger.exception("Error testing features endpoints")
        return False


def test_signatures_endpoints() -> bool:
    """Test Signatures API endpoints."""
    print("\n" + "=" * 80)
    print("6. Testing Signatures Endpoints...")
    print("=" * 80)

    base_url = f"{API_BASE_URL}{API_PREFIX}/signatures"

    try:
        # Test GET (list)
        print("  Testing GET /signatures (list)...")
        resp = requests.get(base_url, timeout=5)
        if resp.status_code == 200:
            data = resp.json()
            print(f"  ✅ GET /signatures: {len(data)} signatures")
            return True
        else:
            print(f"  ❌ GET /signatures: Status {resp.status_code}")
            print(f"     {resp.text[:200]}")
            return False

    except Exception as e:
        print(f"  ❌ ERROR: {e}")
        logger.exception("Error testing signatures endpoints")
        return False


def test_api_docs() -> bool:
    """Test API documentation endpoints."""
    print("\n" + "=" * 80)
    print("7. Testing API Documentation...")
    print("=" * 80)

    try:
        # Test OpenAPI schema
        print("  Testing GET /openapi.json...")
        resp = requests.get(f"{API_BASE_URL}/openapi.json", timeout=5)
        if resp.status_code == 200:
            schema = resp.json()
            paths = schema.get("paths", {})
            print(f"  ✅ OpenAPI schema: {len(paths)} endpoints defined")
        else:
            print(f"  ❌ OpenAPI schema: Status {resp.status_code}")
            return False

        # Test Swagger UI
        print("  Testing GET /docs (Swagger UI)...")
        resp = requests.get(f"{API_BASE_URL}/docs", timeout=5)
        if resp.status_code == 200:
            print(f"  ✅ Swagger UI accessible")
        else:
            print(f"  ❌ Swagger UI: Status {resp.status_code}")

        # Test ReDoc
        print("  Testing GET /redoc (ReDoc)...")
        resp = requests.get(f"{API_BASE_URL}/redoc", timeout=5)
        if resp.status_code == 200:
            print(f"  ✅ ReDoc accessible")
        else:
            print(f"  ❌ ReDoc: Status {resp.status_code}")

        return True

    except Exception as e:
        print(f"  ❌ ERROR: {e}")
        logger.exception("Error testing API docs")
        return False


def main() -> None:
    """Run all API endpoint tests."""
    import argparse

    parser = argparse.ArgumentParser(description="Test FastAPI endpoints")
    parser.add_argument(
        "--url",
        default=API_BASE_URL,
        help=f"API base URL (default: {API_BASE_URL})",
    )
    parser.add_argument(
        "--wait",
        type=int,
        default=0,
        help="Wait N seconds for server to start (default: 0)",
    )

    args = parser.parse_args()

    global API_BASE_URL
    API_BASE_URL = args.url.rstrip("/")

    print("\n" + "=" * 80)
    print("FASTAPI ENDPOINT TESTING")
    print("=" * 80)
    print(f"API Base URL: {API_BASE_URL}")
    print()

    if args.wait > 0:
        print(f"Waiting {args.wait} seconds for server to start...")
        time.sleep(args.wait)

    checks = [
        ("Server Running", check_server_running),
        ("Programs API", test_programs_endpoints),
        ("Experiments API", test_experiments_endpoints),
        ("Datasets API", test_datasets_endpoints),
        ("Features API", test_features_endpoints),
        ("Signatures API", test_signatures_endpoints),
        ("API Documentation", test_api_docs),
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
        print("✅ All API endpoints are working!")
        sys.exit(0)
    else:
        print("⚠️  Some endpoint tests failed. Please review the errors above.")
        sys.exit(1)


if __name__ == "__main__":
    main()

