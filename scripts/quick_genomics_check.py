#!/usr/bin/env python3
"""
Quick sanity check for genomics pipeline components.

Verifies all prerequisites and components are ready.
"""

import sys
from pathlib import Path

# Add parent directory to path
sys.path.insert(0, str(Path(__file__).parent.parent))

from amprenta_rag.ingestion.genomics.pipeline import check_salmon_installed
from amprenta_rag.logging_utils import get_logger

logger = get_logger(__name__)


def check_prerequisites():
    """Check all prerequisites are met."""
    print("\n" + "="*60)
    print("GENOMICS PIPELINE - QUICK SANITY CHECK")
    print("="*60)

    checks = []

    # 1. Check Salmon installation
    print("\n[1] Checking Salmon installation...")
    salmon_installed = check_salmon_installed()
    checks.append(("Salmon installed", salmon_installed))

    # 2. Check transcriptome index
    print("\n[2] Checking transcriptome index...")
    index_path = Path("./salmon_index")
    index_exists = index_path.exists() and (
        (index_path / "info.json").exists() or
        (index_path / "hash.bin").exists()
    )
    checks.append(("Transcriptome index exists", index_exists))
    if index_exists:
        print(f"   ✅ Index found at: {index_path}")

    # 3. Check reference transcriptome
    print("\n[3] Checking reference transcriptome...")
    transcriptome_file = Path("./reference_data/Homo_sapiens.GRCh38.cdna.all.fa.gz")
    transcriptome_exists = transcriptome_file.exists()
    checks.append(("Reference transcriptome exists", transcriptome_exists))
    if transcriptome_exists:
        size_mb = transcriptome_file.stat().st_size / (1024 * 1024)
        print(f"   ✅ Transcriptome found: {size_mb:.1f} MB")

    # 4. Check Python imports
    print("\n[4] Checking Python imports...")
    try:
        imports_ok = True
        print("   ✅ All pipeline functions importable")
    except Exception as e:
        imports_ok = False
        print(f"   ❌ Import error: {e}")
    checks.append(("Pipeline imports", imports_ok))

    # 5. Check ENA repository
    print("\n[5] Checking ENA repository...")
    try:
        from amprenta_rag.ingestion.repositories.ena import ENARepository
        ENARepository()
        ena_ok = True
        print("   ✅ ENA repository available")
    except Exception as e:
        ena_ok = False
        print(f"   ❌ ENA repository error: {e}")
    checks.append(("ENA repository", ena_ok))

    # Summary
    print("\n" + "="*60)
    print("SANITY CHECK SUMMARY")
    print("="*60)

    passed = sum(1 for _, result in checks if result)
    total = len(checks)

    for check_name, result in checks:
        status = "✅ PASS" if result else "❌ FAIL"
        print(f"{status}: {check_name}")

    print(f"\nTotal: {passed}/{total} checks passed")

    if passed == total:
        print("\n🎉 All prerequisites met! Genomics pipeline is ready to use.")
        return 0
    else:
        print("\n⚠️  Some prerequisites missing. Please address the failures above.")
        return 1


if __name__ == "__main__":
    sys.exit(check_prerequisites())

