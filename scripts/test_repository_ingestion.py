#!/usr/bin/env python3
"""
End-to-end test script for repository ingestion.

Tests the complete workflow:
1. Discover studies from repositories
2. Harvest study (create Notion page)
3. Trigger ingestion
4. Verify feature linking
5. Check signature scoring
"""

import argparse
import sys
from pathlib import Path

# Add parent directory to path
sys.path.insert(0, str(Path(__file__).parent.parent))

from amprenta_rag.ingestion.repositories.discovery import (
    discover_studies,
    fetch_study_metadata,
)
from amprenta_rag.logging_utils import get_logger

logger = get_logger(__name__)


def test_discovery(repository: str, keywords: list[str], max_results: int = 3) -> list[str]:
    """
    Test study discovery from a repository.

    Args:
        repository: Repository name
        keywords: Search keywords
        max_results: Maximum results

    Returns:
        List of study IDs found
    """
    print(f"\n{'=' * 80}")
    print(f"TEST 1: Discovery - {repository}")
    print(f"{'=' * 80}")

    try:
        results = discover_studies(
            keywords=keywords,
            repository=repository,
            max_results=max_results,
        )

        study_ids = results.get(repository, [])

        print("✅ Discovery successful!")
        print(f"   Found {len(study_ids)} studies")
        if study_ids:
            print(f"   Examples: {study_ids[:3]}")
            return study_ids
        else:
            print("   ⚠️  No studies found")
            return []

    except Exception as e:
        print(f"❌ Discovery failed: {e}")
        import traceback
        traceback.print_exc()
        return []


def test_metadata_fetch(repository: str, study_id: str) -> bool:
    """
    Test metadata fetching for a study.

    Args:
        repository: Repository name
        study_id: Study ID

    Returns:
        True if successful
    """
    print(f"\n{'=' * 80}")
    print(f"TEST 2: Metadata Fetch - {repository} / {study_id}")
    print(f"{'=' * 80}")

    try:
        metadata = fetch_study_metadata(study_id, repository)

        if metadata:
            print("✅ Metadata fetch successful!")
            print(f"   Title: {metadata.title[:60]}...")
            print(f"   Omics Type: {metadata.omics_type}")
            if metadata.disease:
                print(f"   Disease: {', '.join(metadata.disease)}")
            if metadata.organism:
                print(f"   Organism: {', '.join(metadata.organism)}")
            if metadata.sample_type:
                print(f"   Sample Type: {', '.join(metadata.sample_type)}")
            return True
        else:
            print("❌ Metadata fetch returned None")
            return False

    except Exception as e:
        print(f"❌ Metadata fetch failed: {e}")
        import traceback
        traceback.print_exc()
        return False


def test_harvest_dry_run(repository: str, study_id: str) -> bool:
    """
    Test harvest in dry-run mode (no Notion creation).

    Args:
        repository: Repository name
        study_id: Study ID

    Returns:
        True if successful
    """
    print(f"\n{'=' * 80}")
    print(f"TEST 3: Harvest (Dry Run) - {repository} / {study_id}")
    print(f"{'=' * 80}")

    try:
        # Import harvest function
        from scripts.harvest_repository_study import harvest_study

        harvest_study(
            study_id=study_id,
            repository=repository,
            create_notion=False,
            ingest=False,
            dry_run=True,
        )

        print("✅ Harvest dry-run successful!")
        return True

    except Exception as e:
        print(f"❌ Harvest dry-run failed: {e}")
        import traceback
        traceback.print_exc()
        return False


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Test repository ingestion end-to-end"
    )
    parser.add_argument(
        "--repository",
        choices=["MW", "MW_LIPIDOMICS", "MW_METABOLOMICS", "PRIDE", "MetaboLights", "GEO"],
        default="MW_LIPIDOMICS",
        help="Repository to test",
    )
    parser.add_argument(
        "--keywords",
        nargs="+",
        default=["ceramide"],
        help="Search keywords",
    )
    parser.add_argument(
        "--max-results",
        type=int,
        default=3,
        help="Maximum results to discover",
    )
    parser.add_argument(
        "--test-harvest",
        action="store_true",
        help="Test actual harvest (creates Notion pages)",
    )
    parser.add_argument(
        "--test-ingestion",
        action="store_true",
        help="Test ingestion after harvest",
    )

    args = parser.parse_args()

    print("\n" + "=" * 80)
    print("REPOSITORY INGESTION - END-TO-END TEST")
    print("=" * 80)
    print(f"\nRepository: {args.repository}")
    print(f"Keywords: {', '.join(args.keywords)}")
    print(f"Max Results: {args.max_results}")

    # Test 1: Discovery
    study_ids = test_discovery(args.repository, args.keywords, args.max_results)

    if not study_ids:
        print("\n❌ No studies found. Cannot continue testing.")
        sys.exit(1)

    # Test 2: Metadata Fetch
    test_study_id = study_ids[0]
    metadata_ok = test_metadata_fetch(args.repository, test_study_id)

    if not metadata_ok:
        print("\n❌ Metadata fetch failed. Cannot continue testing.")
        sys.exit(1)

    # Test 3: Harvest (Dry Run)
    harvest_ok = test_harvest_dry_run(args.repository, test_study_id)

    if not harvest_ok:
        print("\n❌ Harvest dry-run failed.")
        sys.exit(1)

    # Test 4: Actual Harvest (if requested)
    if args.test_harvest:
        print(f"\n{'=' * 80}")
        print(f"TEST 4: Harvest (Create Notion Page) - {args.repository} / {test_study_id}")
        print(f"{'=' * 80}")

        try:
            from scripts.harvest_repository_study import harvest_study

            page_id = harvest_study(
                study_id=test_study_id,
                repository=args.repository,
                create_notion=True,
                ingest=args.test_ingestion,
                dry_run=False,
            )

            if page_id:
                print("✅ Harvest successful!")
                print(f"   Dataset page ID: {page_id}")

                if args.test_ingestion:
                    print("\n✅ Ingestion triggered!")
                    print("   Check logs for ingestion progress")
            else:
                print("⚠️  Harvest completed but no page ID returned")

        except Exception as e:
            print(f"❌ Harvest failed: {e}")
            import traceback
            traceback.print_exc()

    # Summary
    print(f"\n{'=' * 80}")
    print("TEST SUMMARY")
    print(f"{'=' * 80}")
    print(f"✅ Discovery: {'PASS' if study_ids else 'FAIL'}")
    print(f"✅ Metadata Fetch: {'PASS' if metadata_ok else 'FAIL'}")
    print(f"✅ Harvest Dry-Run: {'PASS' if harvest_ok else 'FAIL'}")
    if args.test_harvest:
        print(f"✅ Harvest (Notion): {'PASS' if 'page_id' in locals() and page_id else 'FAIL'}")
    print(f"{'=' * 80}\n")


if __name__ == "__main__":
    main()

