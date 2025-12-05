#!/usr/bin/env python3
"""
Test that all TIER 3 modules can be imported without errors.

This script validates module structure and imports without requiring
all dependencies to be installed or configured.
"""

import sys
from pathlib import Path

# Add parent directory to path
sys.path.insert(0, str(Path(__file__).parent.parent))


def test_import(module_path: str, description: str) -> tuple[bool, str]:
    """Test importing a module."""
    try:
        __import__(module_path)
        return True, f"✅ {description}"
    except ImportError as e:
        return False, f"❌ {description}: {str(e)}"
    except Exception as e:
        return False, f"⚠️  {description}: {str(e)}"


def main() -> None:
    """Test all module imports."""
    print("=" * 80)
    print("TIER 3 MODULE STRUCTURE TEST")
    print("=" * 80)
    print()
    
    tests = [
        # Database layer
        ("amprenta_rag.database.base", "Database base module"),
        ("amprenta_rag.database.models", "Database models"),
        ("amprenta_rag.database", "Database package"),
        
        # Models
        ("amprenta_rag.models.domain", "Domain models"),
        
        # Migration
        ("amprenta_rag.migration.dual_write", "Dual-write migration"),
        
        # RAG integration
        ("amprenta_rag.rag.postgres_builder", "Postgres RAG builder"),
        ("amprenta_rag.rag.postgres_resolver", "Postgres resolver"),
        ("amprenta_rag.rag.hybrid_chunk_collection", "Hybrid chunk collection"),
        
        # Ingestion integration
        ("amprenta_rag.ingestion.postgres_integration", "Postgres ingestion integration"),
        
        # Config
        ("amprenta_rag.config", "Configuration"),
    ]
    
    results = []
    
    print("Testing module imports (lazy initialization)...")
    print("-" * 80)
    
    for module_path, description in tests:
        success, message = test_import(module_path, description)
        print(f"  {message}")
        results.append((description, success))
    
    print()
    
    # Summary
    print("=" * 80)
    print("SUMMARY")
    print("=" * 80)
    
    passed = sum(1 for _, success in results if success)
    total = len(results)
    
    for description, success in results:
        status = "✅" if success else "❌"
        print(f"  {status} {description}")
    
    print()
    print(f"  Total: {passed}/{total} modules import successfully")
    print()
    
    if passed == total:
        print("✅ All modules can be imported!")
        sys.exit(0)
    else:
        print("⚠️  Some modules failed to import (may be expected if dependencies missing)")
        sys.exit(1)


if __name__ == "__main__":
    main()

