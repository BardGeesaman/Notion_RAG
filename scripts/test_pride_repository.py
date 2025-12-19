#!/usr/bin/env python3
"""
Test PRIDE repository implementation with strict protocol.

Tests:
1. Search for projects
2. Fetch project metadata
3. Fetch project files
"""

import sys
from pathlib import Path

# Add parent directory to path
sys.path.insert(0, str(Path(__file__).parent.parent))

from amprenta_rag.ingestion.repositories.pride import PRIDERepository
from amprenta_rag.logging_utils import get_logger

logger = get_logger(__name__)


def test_pride_search():
    """Test PRIDE project search."""
    print("\n" + "="*70)
    print("TEST 1: Search PRIDE Projects")
    print("="*70 + "\n")

    repo = PRIDERepository()

    # Search for a common term
    keywords = ["ALS"]
    print(f"Searching for projects with keywords: {keywords}")

    try:
        project_ids = repo.search_studies(keywords=keywords, max_results=5)

        print("\n✅ Search successful!")
        print(f"   Found {len(project_ids)} projects:")
        for i, pid in enumerate(project_ids, 1):
            print(f"   {i}. {pid}")

        return project_ids[0] if project_ids else None

    except Exception as e:
        print(f"\n❌ Search failed: {e}")
        import traceback
        traceback.print_exc()
        return None


def test_pride_metadata(project_id: str):
    """Test fetching PRIDE project metadata."""
    print("\n" + "="*70)
    print("TEST 2: Fetch Project Metadata")
    print("="*70 + "\n")

    if not project_id:
        print("⚠️  Skipping - no project ID available from search test")
        return None

    repo = PRIDERepository()

    print(f"Fetching metadata for project: {project_id}")

    try:
        metadata = repo.fetch_study_metadata(project_id)

        if metadata:
            print("\n✅ Metadata fetch successful!")
            print(f"   Title: {metadata.title}")
            print(f"   Description: {metadata.summary[:100] if metadata.summary else 'N/A'}...")
            print(f"   Organism: {', '.join(metadata.organism) if metadata.organism else 'N/A'}")
            print(f"   Disease: {', '.join(metadata.disease) if metadata.disease else 'N/A'}")
            print(f"   DOI: {metadata.doi or 'N/A'}")
            print(f"   Publication Date: {metadata.publication_date or 'N/A'}")
            return project_id
        else:
            print("\n❌ Metadata fetch returned None")
            return None

    except Exception as e:
        print(f"\n❌ Metadata fetch failed: {e}")
        import traceback
        traceback.print_exc()
        return None


def test_pride_files(project_id: str):
    """Test fetching PRIDE project files."""
    print("\n" + "="*70)
    print("TEST 3: Fetch Project Files")
    print("="*70 + "\n")

    if not project_id:
        print("⚠️  Skipping - no project ID available")
        return

    repo = PRIDERepository()

    print(f"Fetching files for project: {project_id}")

    try:
        # Try to get TSV/CSV files for protein extraction
        data_files = repo.fetch_study_data_files(
            study_id=project_id,
            file_types=["TSV", "CSV", "TXT"],
        )

        if data_files:
            print("\n✅ Files fetch successful!")
            print(f"   Found {len(data_files)} files:")
            for i, file_info in enumerate(data_files[:5], 1):  # Show first 5
                print(f"   {i}. {file_info.filename}")
                print(f"      Type: {file_info.file_type}, Size: {file_info.size_bytes or 'N/A'} bytes")

            if len(data_files) > 5:
                print(f"   ... and {len(data_files) - 5} more files")
        else:
            print("\n⚠️  No TSV/CSV files found (this is OK - project may only have raw files)")

    except Exception as e:
        print(f"\n❌ Files fetch failed: {e}")
        import traceback
        traceback.print_exc()


def main():
    """Run all tests."""
    print("\n" + "="*70)
    print("PRIDE REPOSITORY TEST - STRICT PROTOCOL")
    print("="*70)

    # Test 1: Search
    project_id = test_pride_search()

    # Test 2: Metadata (use a known project if search fails)
    if not project_id:
        print("\n⚠️  Using known test project: PXD000001")
        project_id = "PXD000001"

    project_id = test_pride_metadata(project_id)

    # Test 3: Files
    test_pride_files(project_id)

    print("\n" + "="*70)
    print("TEST COMPLETE")
    print("="*70 + "\n")


if __name__ == "__main__":
    main()

