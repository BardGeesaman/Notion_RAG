#!/usr/bin/env python3
"""
Test MetaboLights repository implementation with strict protocol.

Tests:
1. Fetch study metadata
2. Fetch study files
"""

import sys
from pathlib import Path

# Add parent directory to path
sys.path.insert(0, str(Path(__file__).parent.parent))

from amprenta_rag.ingestion.repositories.metabolights import MetaboLightsRepository
from amprenta_rag.logging_utils import get_logger

logger = get_logger(__name__)


def test_metabolights_metadata(study_id: str):
    """Test fetching MetaboLights study metadata."""
    print("\n" + "="*70)
    print("TEST 1: Fetch Study Metadata")
    print("="*70 + "\n")
    
    repo = MetaboLightsRepository()
    
    print(f"Fetching metadata for study: {study_id}")
    
    try:
        metadata = repo.fetch_study_metadata(study_id)
        
        if metadata:
            print(f"\n✅ Metadata fetch successful!")
            print(f"   Title: {metadata.title}")
            print(f"   Description: {metadata.summary[:100] if metadata.summary else 'N/A'}...")
            print(f"   Organism: {', '.join(metadata.organism) if metadata.organism else 'N/A'}")
            print(f"   Disease: {', '.join(metadata.disease) if metadata.disease else 'N/A'}")
            print(f"   DOI: {metadata.doi or 'N/A'}")
            print(f"   Publication Date: {metadata.publication_date or 'N/A'}")
            return study_id
        else:
            print(f"\n❌ Metadata fetch returned None")
            return None
            
    except Exception as e:
        print(f"\n❌ Metadata fetch failed: {e}")
        import traceback
        traceback.print_exc()
        return None


def test_metabolights_files(study_id: str):
    """Test fetching MetaboLights study files."""
    print("\n" + "="*70)
    print("TEST 2: Fetch Study Files")
    print("="*70 + "\n")
    
    if not study_id:
        print("⚠️  Skipping - no study ID available")
        return
    
    repo = MetaboLightsRepository()
    
    print(f"Fetching files for study: {study_id}")
    
    try:
        # Get all files first
        all_files = repo.fetch_study_data_files(study_id, file_types=None)
        
        if all_files:
            print(f"\n✅ Files fetch successful!")
            print(f"   Found {len(all_files)} total files:")
            for i, file_info in enumerate(all_files[:10], 1):  # Show first 10
                print(f"   {i}. {file_info.filename}")
                print(f"      Type: {file_info.file_type}, Size: {file_info.size_bytes or 'N/A'} bytes")
            
            if len(all_files) > 10:
                print(f"   ... and {len(all_files) - 10} more files")
        
        # Try to get TSV/CSV files for metabolite extraction
        data_files = repo.fetch_study_data_files(
            study_id,
            file_types=["TSV", "CSV", "TXT"],
        )
        
        print(f"\n   TSV/CSV files found: {len(data_files)}")
        if data_files:
            for i, file_info in enumerate(data_files[:5], 1):
                print(f"   {i}. {file_info.filename}")
        else:
            print("   (No TSV/CSV files - this is OK, may have other formats)")
            
    except Exception as e:
        print(f"\n❌ Files fetch failed: {e}")
        import traceback
        traceback.print_exc()


def main():
    """Run all tests."""
    print("\n" + "="*70)
    print("METABOLIGHTS REPOSITORY TEST - STRICT PROTOCOL")
    print("="*70)
    
    # Use a known public MetaboLights study
    study_id = "MTBLS1"  # Known public study
    
    print(f"\nUsing MetaboLights study: {study_id}")
    
    # Test 1: Metadata
    study_id = test_metabolights_metadata(study_id)
    
    # Test 2: Files
    test_metabolights_files(study_id)
    
    print("\n" + "="*70)
    print("TEST COMPLETE")
    print("="*70 + "\n")


if __name__ == "__main__":
    main()

