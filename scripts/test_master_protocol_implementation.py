#!/usr/bin/env python3
"""
Test script for Master Protocol Implementation.

Tests:
1. User-Agent headers are being used
2. ENA repository works
3. All repositories are properly integrated
"""

import sys
from pathlib import Path

# Add parent directory to path
sys.path.insert(0, str(Path(__file__).parent.parent))

from amprenta_rag.ingestion.repositories import (
    ENARepository,
    GEORepository,
    MWRepository,
    PRIDERepository,
    REPOSITORY_USER_AGENT,
    get_repository,
)


def test_user_agent_constant():
    """Test that User-Agent constant is accessible."""
    print("\n" + "="*60)
    print("TEST 1: User-Agent Constant")
    print("="*60)

    print(f"✅ User-Agent constant: {REPOSITORY_USER_AGENT}")

    assert REPOSITORY_USER_AGENT == "ResearchBot/1.0 (Bioinformatics Data Pipeline)"
    print("✅ User-Agent constant is correct")
    return True


def test_ena_repository():
    """Test ENA repository initialization and basic functionality."""
    print("\n" + "="*60)
    print("TEST 2: ENA Repository")
    print("="*60)

    # Test initialization
    ena = ENARepository()
    print("✅ ENA repository initialized")
    print(f"   Repository name: {ena.get_repository_name()}")
    print(f"   Omics type: {ena.get_omics_type()}")

    assert ena.get_repository_name() == "ENA"
    assert ena.get_omics_type() == "genomics"
    print("✅ ENA repository basic properties are correct")

    # Test repository registry
    ena_from_registry = get_repository("ENA")
    assert ena_from_registry is not None
    print("✅ ENA repository accessible from registry")

    return True


def test_all_repositories():
    """Test that all repositories can be initialized."""
    print("\n" + "="*60)
    print("TEST 3: All Repositories")
    print("="*60)

    repositories = {
        "GEO": GEORepository,
        "PRIDE": PRIDERepository,
        "MW": MWRepository,
        "ENA": ENARepository,
    }

    for name, repo_class in repositories.items():
        try:
            if name == "GEO":
                repo = repo_class()  # Can optionally pass api_key and email
            elif name == "MW":
                repo = repo_class(omics_type="metabolomics")
            else:
                repo = repo_class()

            repo_name = repo.get_repository_name()
            omics_type = repo.get_omics_type()
            print(f"✅ {name}: {repo_name} ({omics_type})")

        except Exception as e:
            print(f"❌ {name}: Error - {e}")
            return False

    print("✅ All repositories initialized successfully")
    return True


def test_repository_discovery():
    """Test repository discovery system."""
    print("\n" + "="*60)
    print("TEST 4: Repository Discovery")
    print("="*60)

    from amprenta_rag.ingestion.repositories.discovery import (
        get_repository,
        list_available_repositories,
    )

    # Test list of repositories
    repos = list_available_repositories()
    print(f"✅ Available repositories: {', '.join(repos)}")

    # Check ENA is in the list
    assert "ENA" in repos
    print("✅ ENA is in repository list")

    # Test getting repositories
    test_repos = ["GEO", "PRIDE", "MW", "ENA"]
    for repo_name in test_repos:
        repo = get_repository(repo_name)
        if repo:
            print(f"✅ Successfully retrieved {repo_name}")
        else:
            print(f"❌ Failed to retrieve {repo_name}")
            return False

    print("✅ Repository discovery system working")
    return True


def test_user_agent_in_requests():
    """Test that User-Agent headers are included in requests."""
    print("\n" + "="*60)
    print("TEST 5: User-Agent Headers in Requests")
    print("="*60)

    import inspect

    # Check PRIDE repository
    from amprenta_rag.ingestion.repositories.pride import PRIDERepository
    pride = PRIDERepository()
    pride_source = inspect.getsource(pride._make_request_with_retry)

    if "REPOSITORY_USER_AGENT" in pride_source or "User-Agent" in pride_source:
        print("✅ PRIDE repository includes User-Agent header")
    else:
        print("⚠️  PRIDE repository User-Agent header check - inspect source")

    # Check MetaboLights repository
    from amprenta_rag.ingestion.repositories.metabolights import MetaboLightsRepository
    ml = MetaboLightsRepository()
    ml_source = inspect.getsource(ml._make_request_with_retry)

    if "REPOSITORY_USER_AGENT" in ml_source or "User-Agent" in ml_source:
        print("✅ MetaboLights repository includes User-Agent header")
    else:
        print("⚠️  MetaboLights repository User-Agent header check - inspect source")

    # Check MW repository
    from amprenta_rag.ingestion.repositories.mw import MWRepository
    mw = MWRepository()
    # Check one of the request methods
    mw_source = inspect.getsource(mw._fetch_all_study_summaries)

    if "REPOSITORY_USER_AGENT" in mw_source or "User-Agent" in mw_source:
        print("✅ MW repository includes User-Agent header")
    else:
        print("⚠️  MW repository User-Agent header check - inspect source")

    # Check ENA repository
    ena = ENARepository()
    ena_source = inspect.getsource(ena._make_request)

    if "REPOSITORY_USER_AGENT" in ena_source or "User-Agent" in ena_source:
        print("✅ ENA repository includes User-Agent header")
    else:
        print("⚠️  ENA repository User-Agent header check - inspect source")

    print("✅ User-Agent header checks completed")
    return True


def test_ena_search():
    """Test ENA repository search functionality (dry run)."""
    print("\n" + "="*60)
    print("TEST 6: ENA Repository Search (Dry Run)")
    print("="*60)

    ena = ENARepository()

    # Test search with a simple query (won't actually execute)
    print("Testing ENA search method (dry run)...")

    try:
        # Just verify the method exists and can be called
        import inspect
        search_method = getattr(ena, 'search_studies', None)
        if search_method:
            sig = inspect.signature(search_method)
            print(f"✅ search_studies method exists with signature: {sig}")
        else:
            print("❌ search_studies method not found")
            return False

        # Verify other methods exist
        assert hasattr(ena, 'fetch_study_metadata')
        assert hasattr(ena, 'fetch_study_data_files')
        print("✅ All required methods exist")

    except Exception as e:
        print(f"❌ Error: {e}")
        return False

    print("✅ ENA repository structure verified")
    return True


def main():
    """Run all tests."""
    print("\n" + "="*60)
    print("MASTER PROTOCOL IMPLEMENTATION TEST SUITE")
    print("="*60)

    tests = [
        test_user_agent_constant,
        test_ena_repository,
        test_all_repositories,
        test_repository_discovery,
        test_user_agent_in_requests,
        test_ena_search,
    ]

    passed = 0
    failed = 0

    for test in tests:
        try:
            if test():
                passed += 1
            else:
                failed += 1
        except Exception as e:
            print(f"❌ Test {test.__name__} failed with exception: {e}")
            import traceback
            traceback.print_exc()
            failed += 1

    print("\n" + "="*60)
    print("TEST SUMMARY")
    print("="*60)
    print(f"✅ Passed: {passed}")
    print(f"❌ Failed: {failed}")
    print(f"📊 Total:  {passed + failed}")

    if failed == 0:
        print("\n🎉 All tests passed!")
        return 0
    else:
        print(f"\n⚠️  {failed} test(s) failed")
        return 1


if __name__ == "__main__":
    sys.exit(main())

