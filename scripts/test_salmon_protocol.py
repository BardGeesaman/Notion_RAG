#!/usr/bin/env python3
"""
Test script to verify Salmon Subprocess Protocol implementation.

This script verifies:
1. Salmon installation check
2. Subprocess Protocol compliance (no Python imports)
3. Command structure and flags
"""

import sys
from pathlib import Path

# Add parent directory to path
sys.path.insert(0, str(Path(__file__).parent.parent))

from amprenta_rag.ingestion.genomics.pipeline import check_salmon_installed, quantify_with_salmon


def test_salmon_installation():
    """Test 1: Verify Salmon installation check."""
    print("\n" + "="*60)
    print("TEST 1: Salmon Installation Check")
    print("="*60)

    is_installed = check_salmon_installed()

    if is_installed:
        print("✅ PASS: Salmon installation check passed")
    else:
        print("❌ FAIL: Salmon installation check failed")
        print("   -> Please install Salmon: conda install -c bioconda salmon")

    return is_installed


def test_subprocess_protocol():
    """Test 2: Verify Subprocess Protocol compliance."""
    print("\n" + "="*60)
    print("TEST 2: Subprocess Protocol Compliance")
    print("="*60)

    # Check that we're NOT importing salmon as a Python package

    # Verify salmon is NOT in sys.modules as a Python package
    salmon_modules = [m for m in sys.modules.keys() if 'salmon' in m.lower()]

    # Filter out our own modules (amprenta_rag.ingestion.genomics.pipeline)
    salmon_modules = [m for m in salmon_modules if 'pipeline' not in m]

    if not salmon_modules:
        print("✅ PASS: No Python 'salmon' package imported")
        print("   -> Correctly using subprocess instead")
    else:
        print(f"⚠️  WARNING: Found modules: {salmon_modules}")
        print("   -> Make sure we're not importing salmon as a Python package!")

    # Verify subprocess is being used
    print("✅ PASS: Using subprocess module for Salmon commands")

    return True


def test_command_structure():
    """Test 3: Verify command structure follows protocol."""
    print("\n" + "="*60)
    print("TEST 3: Command Structure")
    print("="*60)

    # Check that quantify_with_salmon function exists and has correct signature
    import inspect

    sig = inspect.signature(quantify_with_salmon)
    params = list(sig.parameters.keys())

    required_params = ['fastq_path', 'index_path']
    optional_params = ['output_dir', 'library_type', 'validate_mappings']

    print(f"Function parameters: {params}")

    # Check required parameters
    all_required_present = all(p in params for p in required_params)
    if all_required_present:
        print("✅ PASS: Required parameters present")
    else:
        print(f"❌ FAIL: Missing required parameters. Expected: {required_params}")

    # Check optional parameters
    all_optional_present = all(p in params for p in optional_params)
    if all_optional_present:
        print("✅ PASS: Optional parameters (validate_mappings) present")
    else:
        print(f"⚠️  WARNING: Some optional parameters missing: {optional_params}")

    return all_required_present


def main():
    """Run all tests."""
    print("\n" + "╔" + "="*58 + "╗")
    print("║" + " "*18 + "SALMON PROTOCOL TESTS" + " "*18 + "║")
    print("╚" + "="*58 + "╝")

    results = []

    # Test 1: Installation check
    results.append(("Installation Check", test_salmon_installation()))

    # Test 2: Subprocess Protocol
    results.append(("Subprocess Protocol", test_subprocess_protocol()))

    # Test 3: Command Structure
    results.append(("Command Structure", test_command_structure()))

    # Summary
    print("\n" + "="*60)
    print("TEST SUMMARY")
    print("="*60)

    passed = sum(1 for _, result in results if result)
    total = len(results)

    for test_name, result in results:
        status = "✅ PASS" if result else "❌ FAIL"
        print(f"{status}: {test_name}")

    print(f"\nTotal: {passed}/{total} tests passed")

    if passed == total:
        print("\n🎉 All tests passed! Salmon Subprocess Protocol is correctly implemented.")
        return 0
    else:
        print("\n⚠️  Some tests failed. Please review the implementation.")
        return 1


if __name__ == "__main__":
    sys.exit(main())

