#!/usr/bin/env python3
"""
Test Metabolomics Workbench metabolite extraction using the cleaner /data endpoint.

This tests the robust REST API approach for extracting metabolites from MW studies.
"""

import sys
from pathlib import Path

# Add parent directory to path
sys.path.insert(0, str(Path(__file__).parent.parent))

from amprenta_rag.ingestion.repository_feature_extraction import extract_mw_metabolites_from_data_endpoint
from amprenta_rag.logging_utils import get_logger

logger = get_logger(__name__)


def test_mw_metabolite_extraction(study_id: str = "ST000001"):
    """Test extracting metabolites from a Metabolomics Workbench study."""
    print("\n" + "="*70)
    print("TEST: Metabolomics Workbench Metabolite Extraction")
    print("="*70 + "\n")

    print(f"Extracting metabolites from study: {study_id}")
    print(f"Using cleaner REST API /data endpoint\n")

    try:
        metabolite_set = extract_mw_metabolites_from_data_endpoint(
            study_id=study_id,
        )

        if metabolite_set:
            print(f"\n✅ Extraction successful!")
            print(f"   Found {len(metabolite_set)} unique metabolites")
            print(f"\n   Sample metabolites (first 20):")
            for i, metabolite in enumerate(sorted(metabolite_set)[:20], 1):
                print(f"      {i}. {metabolite}")

            if len(metabolite_set) > 20:
                print(f"   ... and {len(metabolite_set) - 20} more metabolites")
        else:
            print(f"\n⚠️  No metabolites extracted")
            print("   This may be expected if:")
            print("   - The study doesn't have metabolite identification data")
            print("   - The study structure is different than expected")

        return metabolite_set

    except Exception as e:
        print(f"\n❌ Extraction failed: {e}")
        import traceback
        traceback.print_exc()
        return set()


def main():
    """Run test."""
    import sys

    print("\n" + "="*70)
    print("METABOLOMICS WORKBENCH METABOLITE EXTRACTION TEST")
    print("="*70)

    # Get study ID from command line argument or use default
    study_id = sys.argv[1] if len(sys.argv) > 1 else "ST000001"

    metabolite_set = test_mw_metabolite_extraction(study_id)

    print("\n" + "="*70)
    print("TEST COMPLETE")
    print("="*70 + "\n")

    if metabolite_set:
        print(f"✅ Successfully extracted {len(metabolite_set)} metabolites using REST API!")
        print(f"\n   Advantages over MetaboLights:")
        print(f"   • Clean JSON response (no file parsing needed)")
        print(f"   • More stable API (NIH-hosted)")
        print(f"   • No brittle endpoints or 500 errors")
    else:
        print("⚠️  No metabolites extracted")
        print("\n   Try a different study ID:")
        print(f"   python scripts/test_mw_metabolite_extraction.py ST000002")


if __name__ == "__main__":
    main()

