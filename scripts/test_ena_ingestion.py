#!/usr/bin/env python3
"""
Test ENA repository ingestion.

Tests:
1. Search for ENA read runs
2. Fetch metadata for a read run
3. Fetch FASTQ file links
"""

import sys
from pathlib import Path

# Add parent directory to path
sys.path.insert(0, str(Path(__file__).parent.parent))

from amprenta_rag.ingestion.repositories import ENARepository
from amprenta_rag.logging_utils import get_logger

logger = get_logger(__name__)


def test_ena_search():
    """Test searching for ENA read runs."""
    print("\n" + "="*60)
    print("TEST 1: ENA Search")
    print("="*60)

    ena = ENARepository()

    # Search for a common organism (e.g., human)
    keywords = ["Homo sapiens"]
    print(f"Searching for: {keywords}")
    print("(This may take a moment due to rate limiting...)\n")

    try:
        study_ids = ena.search_studies(
            keywords=keywords,
            max_results=5,
        )

        if study_ids:
            print(f"✅ Found {len(study_ids)} read runs:")
            for i, study_id in enumerate(study_ids[:5], 1):
                print(f"   {i}. {study_id}")
            return study_ids[0] if study_ids else None
        else:
            print("⚠️  No studies found")
            return None

    except Exception as e:
        print(f"❌ Error during search: {e}")
        import traceback
        traceback.print_exc()
        return None


def test_ena_metadata(study_id: str):
    """Test fetching metadata for an ENA read run."""
    print("\n" + "="*60)
    print("TEST 2: ENA Metadata Fetching")
    print("="*60)

    if not study_id:
        print("⚠️  Skipping metadata test - no study ID available")
        return None

    ena = ENARepository()

    print(f"Fetching metadata for: {study_id}")
    print("(This may take a moment due to rate limiting...)\n")

    try:
        metadata = ena.fetch_study_metadata(study_id)

        if metadata:
            print("✅ Successfully fetched metadata:")
            print(f"   Title: {metadata.title}")
            print(f"   Repository: {metadata.repository}")
            print(f"   Omics Type: {metadata.omics_type}")

            if metadata.organism:
                print(f"   Organism: {', '.join(metadata.organism)}")

            if metadata.platform:
                print(f"   Platform: {metadata.platform}")

            if metadata.raw_metadata:
                run_data = metadata.raw_metadata
                print("\n   Raw Metadata Fields:")
                print(f"   - Run Accession: {run_data.get('run_accession', 'N/A')}")
                print(f"   - Study Accession: {run_data.get('study_accession', 'N/A')}")
                print(f"   - Experiment Accession: {run_data.get('experiment_accession', 'N/A')}")
                print(f"   - Sample Accession: {run_data.get('sample_accession', 'N/A')}")
                print(f"   - Scientific Name: {run_data.get('scientific_name', 'N/A')}")
                print(f"   - Instrument Platform: {run_data.get('instrument_platform', 'N/A')}")
                print(f"   - Library Strategy: {run_data.get('library_strategy', 'N/A')}")

                fastq_ftp = run_data.get('fastq_ftp', '')
                if fastq_ftp:
                    links = [link.strip() for link in fastq_ftp.split(";") if link.strip()]
                    print(f"   - FASTQ Files: {len(links)} file(s) available")

            return metadata
        else:
            print("❌ Failed to fetch metadata")
            return None

    except Exception as e:
        print(f"❌ Error fetching metadata: {e}")
        import traceback
        traceback.print_exc()
        return None


def test_ena_data_files(study_id: str):
    """Test fetching FASTQ file links for an ENA read run."""
    print("\n" + "="*60)
    print("TEST 3: ENA Data Files (FASTQ Links)")
    print("="*60)

    if not study_id:
        print("⚠️  Skipping data files test - no study ID available")
        return None

    ena = ENARepository()

    print(f"Fetching data files for: {study_id}")
    print("(This may take a moment due to rate limiting...)\n")

    try:
        data_files = ena.fetch_study_data_files(study_id)

        if data_files:
            print(f"✅ Found {len(data_files)} data file(s):")
            for i, data_file in enumerate(data_files[:5], 1):
                print(f"\n   File {i}:")
                print(f"   - ID: {data_file.file_id}")
                print(f"   - Filename: {data_file.filename}")
                print(f"   - Type: {data_file.file_type}")
                print(f"   - URL: {data_file.download_url[:80]}..." if len(data_file.download_url) > 80 else f"   - URL: {data_file.download_url}")
                if data_file.description:
                    print(f"   - Description: {data_file.description}")

            if len(data_files) > 5:
                print(f"\n   ... and {len(data_files) - 5} more file(s)")

            print(f"\n✅ Successfully generated {len(data_files)} FASTQ download link(s)")
            print("   Note: Links are generated but files are not downloaded (Master Protocol)")

            return data_files
        else:
            print("⚠️  No data files found")
            return None

    except Exception as e:
        print(f"❌ Error fetching data files: {e}")
        import traceback
        traceback.print_exc()
        return None


def test_ena_with_known_accession():
    """Test with a known ENA accession (if provided)."""
    print("\n" + "="*60)
    print("TEST 4: ENA with Known Accession")
    print("="*60)

    # You can modify this to test with a specific ENA run accession
    # Example: ERR2756784 (a known public ENA run)
    known_accession = None  # Set to a specific accession to test

    if not known_accession:
        print("⚠️  No known accession provided - skipping")
        print("   To test with a specific accession, set 'known_accession' in the script")
        return None

    ena = ENARepository()

    print(f"Testing with known accession: {known_accession}\n")

    try:
        # Fetch metadata
        metadata = ena.fetch_study_metadata(known_accession)
        if metadata:
            print(f"✅ Metadata fetched: {metadata.title}")

        # Fetch data files
        data_files = ena.fetch_study_data_files(known_accession)
        if data_files:
            print(f"✅ Found {len(data_files)} data file(s)")

        return True

    except Exception as e:
        print(f"❌ Error: {e}")
        import traceback
        traceback.print_exc()
        return False


def main():
    """Run all ENA ingestion tests."""
    print("\n" + "="*60)
    print("ENA REPOSITORY INGESTION TEST SUITE")
    print("="*60)
    print("\nTesting ENA repository following Master Protocol:")
    print("  • Uses ENA Browser API")
    print("  • Generates FASTQ FTP links (does not download)")
    print("  • Rate limiting: 1 second between requests")
    print("  • User-Agent headers included")

    results = {}

    # Test 1: Search
    study_id = test_ena_search()
    results['search'] = study_id is not None

    # Test 2: Metadata (if we got a study ID from search)
    metadata = None
    if study_id:
        metadata = test_ena_metadata(study_id)
        results['metadata'] = metadata is not None

    # Test 3: Data Files (if we got a study ID)
    data_files = None
    if study_id:
        data_files = test_ena_data_files(study_id)
        results['data_files'] = data_files is not None

    # Test 4: Known accession (optional)
    known_result = test_ena_with_known_accession()
    if known_result is not None:
        results['known_accession'] = known_result

    # Summary
    print("\n" + "="*60)
    print("TEST SUMMARY")
    print("="*60)

    for test_name, passed in results.items():
        status = "✅ PASS" if passed else "❌ FAIL"
        print(f"{status}: {test_name}")

    total = len(results)
    passed = sum(1 for v in results.values() if v)

    print(f"\n📊 Results: {passed}/{total} tests passed")

    if passed == total:
        print("\n🎉 All ENA ingestion tests passed!")
        return 0
    elif passed > 0:
        print(f"\n⚠️  {total - passed} test(s) failed or skipped")
        return 1
    else:
        print("\n❌ All tests failed")
        return 1


if __name__ == "__main__":
    sys.exit(main())

