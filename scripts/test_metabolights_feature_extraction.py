#!/usr/bin/env python3
"""
Test MetaboLights metabolite extraction.

Tests the extraction of metabolites from MetaboLights ISA-Tab files.
"""

import sys
from pathlib import Path

# Add parent directory to path
sys.path.insert(0, str(Path(__file__).parent.parent))

from amprenta_rag.ingestion.repository_feature_extraction import extract_metabolights_metabolites_from_isa_tab
from amprenta_rag.logging_utils import get_logger

logger = get_logger(__name__)


def test_metabolights_metabolite_extraction(study_id: str = "MTBLS1"):
    """Test extracting metabolites from a MetaboLights study."""
    print("\n" + "="*70)
    print("TEST: MetaboLights Metabolite Extraction")
    print("="*70 + "\n")
    
    print(f"Extracting metabolites from study: {study_id}")
    
    try:
        metabolite_set = extract_metabolights_metabolites_from_isa_tab(
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
            print("   This may be expected if the study doesn't have accessible metabolite data files")
        
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
    print("METABOLIGHTS METABOLITE EXTRACTION TEST")
    print("="*70)
    
    # Get study ID from command line argument or use default
    study_id = sys.argv[1] if len(sys.argv) > 1 else "MTBLS1"
    
    metabolite_set = test_metabolights_metabolite_extraction(study_id)
    
    print("\n" + "="*70)
    print("TEST COMPLETE")
    print("="*70 + "\n")
    
    if metabolite_set:
        print(f"✅ Successfully extracted {len(metabolite_set)} metabolites!")
    else:
        print("⚠️  No metabolites extracted - check if study has accessible data files")


if __name__ == "__main__":
    main()

