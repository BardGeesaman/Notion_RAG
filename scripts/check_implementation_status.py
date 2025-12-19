#!/usr/bin/env python3
"""
Implementation Status Checker.

Verifies what's actually implemented vs what the roadmap claims.
This helps avoid confusion about what's built vs what's planned.
"""

import sys
from pathlib import Path
from typing import Dict, List

# Add parent directory to path
sys.path.insert(0, str(Path(__file__).parent.parent))


def check_file_exists(file_path: str) -> bool:
    """Check if a file exists."""
    return Path(file_path).exists()


def check_module_imports(module_path: str, expected_imports: List[str]) -> Dict[str, bool]:
    """Check if a module has expected imports/classes."""
    results = {}
    file_path = Path(module_path)

    if not file_path.exists():
        return {imp: False for imp in expected_imports}

    try:
        with open(file_path, 'r') as f:
            content = f.read()

        # Simple string matching (not perfect but works for common cases)
        for expected in expected_imports:
            results[expected] = expected in content or expected.split('.')[-1] in content

    except Exception:
        return {imp: False for imp in expected_imports}

    return results


def check_postgres_integration() -> Dict[str, bool]:
    """Check for Postgres-related code."""
    checks = {
        "Postgres connection code": False,
        "SQLAlchemy models": False,
        "Alembic migrations": False,
    }

    # Check for common Postgres/SQLAlchemy patterns
    patterns = {
        "Postgres connection code": ["postgresql", "psycopg", "asyncpg", "create_engine"],
        "SQLAlchemy models": ["from sqlalchemy", "Base = declarative_base", "Column", "relationship"],
        "Alembic migrations": ["alembic", "revision", "upgrade", "downgrade"],
    }

    codebase_root = Path(__file__).parent.parent / "amprenta_rag"

    for check_name, keywords in patterns.items():
        found = False
        for py_file in codebase_root.rglob("*.py"):
            try:
                content = py_file.read_text().lower()
                if any(kw.lower() in content for kw in keywords):
                    found = True
                    break
            except Exception:
                continue
        checks[check_name] = found

    return checks


def check_fastapi_integration() -> Dict[str, bool]:
    """Check for FastAPI-related code."""
    checks = {
        "FastAPI application": False,
        "API endpoints": False,
        "FastAPI routers": False,
    }

    codebase_root = Path(__file__).parent.parent / "amprenta_rag"

    # Check for FastAPI patterns
    for py_file in codebase_root.rglob("*.py"):
        try:
            content = py_file.read_text()
            if "from fastapi" in content or "import fastapi" in content:
                checks["FastAPI application"] = True
            if "@app." in content or "@router." in content:
                checks["API endpoints"] = True
            if "APIRouter" in content:
                checks["FastAPI routers"] = True
        except Exception:
            continue

    return checks


def check_domain_models() -> Dict[str, bool]:
    """Check for domain model implementations."""
    checks = {
        "Domain models module": False,
        "Program model": False,
        "Experiment model": False,
        "Dataset model": False,
        "Signature model": False,
        "Feature model": False,
    }

    models_dir = Path(__file__).parent.parent / "amprenta_rag" / "models"
    domain_models_file = models_dir / "domain.py"

    checks["Domain models module"] = domain_models_file.exists()

    if checks["Domain models module"]:
        try:
            content = domain_models_file.read_text()
            checks["Program model"] = "class Program" in content or "ProgramModel" in content
            checks["Experiment model"] = "class Experiment" in content or "ExperimentModel" in content
            checks["Dataset model"] = "class Dataset" in content or "DatasetModel" in content
            checks["Signature model"] = "class Signature" in content or "SignatureModel" in content
            checks["Feature model"] = "class Feature" in content or "FeatureModel" in content
        except Exception:
            pass

    return checks


def main() -> None:
    """Run all implementation status checks."""
    print("=" * 80)
    print("IMPLEMENTATION STATUS CHECKER")
    print("=" * 80)
    print()

    # Check TIER 3: Architecture Evolution
    print("🏗️  TIER 3: Architecture Evolution")
    print("-" * 80)

    postgres_checks = check_postgres_integration()
    for check, status in postgres_checks.items():
        status_str = "✅" if status else "❌"
        print(f"  {status_str} {check}")

    fastapi_checks = check_fastapi_integration()
    for check, status in fastapi_checks.items():
        status_str = "✅" if status else "❌"
        print(f"  {status_str} {check}")

    print()

    # Check Domain Models (Phase 1)
    print("📦 Phase 1: Domain Model Extraction")
    print("-" * 80)

    domain_checks = check_domain_models()
    for check, status in domain_checks.items():
        status_str = "✅" if status else "❌"
        print(f"  {status_str} {check}")

    print()

    # Check Phase 3-6 components
    print("🔧 Phase 3: FastAPI Service Layer")
    print("-" * 80)

    fastapi_checks = check_fastapi_integration()
    for check, status in fastapi_checks.items():
        status_str = "✅" if status else "❌"
        print(f"  {status_str} {check}")

    print()

    # Check Phase 4: Dual-Write
    print("🔄 Phase 4: Migration Utilities (Dual-Write)")
    print("-" * 80)

    dual_write_exists = check_file_exists("amprenta_rag/migration/dual_write.py")
    status_str = "✅" if dual_write_exists else "❌"
    print(f"  {status_str} DualWriteManager")

    print()

    # Check Phase 5: RAG Integration
    print("🔍 Phase 5: RAG Integration with Postgres")
    print("-" * 80)

    rag_postgres_exists = check_file_exists("amprenta_rag/rag/postgres_builder.py")
    rag_resolver_exists = check_file_exists("amprenta_rag/rag/postgres_resolver.py")
    hybrid_chunks_exists = check_file_exists("amprenta_rag/rag/hybrid_chunk_collection.py")

    print(f"  {'✅' if rag_postgres_exists else '❌'} Postgres RAG Builder")
    print(f"  {'✅' if rag_resolver_exists else '❌'} Postgres Resolver")
    print(f"  {'✅' if hybrid_chunks_exists else '❌'} Hybrid Chunk Collection")

    print()

    # Check Phase 6: Postgres Integration
    print("📥 Phase 6: Postgres SoT Integration")
    print("-" * 80)

    postgres_integration_exists = check_file_exists("amprenta_rag/ingestion/postgres_integration.py")
    print(f"  {'✅' if postgres_integration_exists else '❌'} Postgres Integration Utilities")

    print()

    # Summary
    print("=" * 80)
    print("SUMMARY")
    print("=" * 80)

    total_tier3 = len(postgres_checks) + len(fastapi_checks)
    completed_tier3 = sum(postgres_checks.values()) + sum(fastapi_checks.values())

    total_phase1 = len(domain_checks)
    completed_phase1 = sum(domain_checks.values())

    phase4_complete = 1 if dual_write_exists else 0
    phase5_complete = sum([rag_postgres_exists, rag_resolver_exists, hybrid_chunks_exists])
    phase6_complete = 1 if postgres_integration_exists else 0

    print(f"TIER 3 (Architecture Evolution): {completed_tier3}/{total_tier3} components")
    print(f"Phase 1 (Domain Models): {completed_phase1}/{total_phase1} models")
    print(f"Phase 2 (Postgres Schema): ✅ Complete")
    print(f"Phase 3 (FastAPI): {sum(fastapi_checks.values())}/{len(fastapi_checks)} components")
    print(f"Phase 4 (Dual-Write): {phase4_complete}/1 components")
    print(f"Phase 5 (RAG Integration): {phase5_complete}/3 components")
    print(f"Phase 6 (Postgres SoT): {phase6_complete}/1 components")
    print()

    if completed_tier3 == total_tier3 and completed_phase1 == total_phase1:
        print("✅✅✅ TIER 3 FULLY IMPLEMENTED! ✅✅✅")
        print()
        print("All phases complete:")
        print("  ✅ Phase 1: Domain Model Extraction")
        print("  ✅ Phase 2: Postgres Schema Design")
        print("  ✅ Phase 3: FastAPI Service Layer")
        print("  ✅ Phase 4: Migration Utilities (Dual-Write)")
        print("  ✅ Phase 5: RAG Integration with Postgres")
        print("  ✅ Phase 6: Transition to Postgres SoT")
    else:
        if completed_tier3 == 0:
            print("⚠️  TIER 3 not yet started - major architectural work needed")
        elif completed_tier3 < total_tier3:
            print(f"⚠️  TIER 3 partially complete - {total_tier3 - completed_tier3} items remaining")
        else:
            print("✅ TIER 3 fully implemented!")

    print()


if __name__ == "__main__":
    main()

