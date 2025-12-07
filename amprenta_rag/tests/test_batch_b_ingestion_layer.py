"""Batch B: Ingestion Layer Quick Test"""

print("Test 1: Imports...")
try:
    from amprenta_rag.ingestion.signature_matching.matching import find_matching_signatures_for_dataset
    from amprenta_rag.ingestion.screening_ingestion import promote_compounds_to_notion
    print("✓ Imports OK")
except Exception as e:
    print(f"✗ Import failed: {e}")
    raise

print("\nTest 2: promote_compounds_to_notion signature...")
try:
    import inspect
    sig = inspect.signature(promote_compounds_to_notion)
    params = list(sig.parameters.keys())
    # Function should accept campaign_id parameter
    assert 'campaign_id' in params, f"Expected 'campaign_id' parameter, got {params}"
    # Test with non-existent campaign_id - should return empty list
    result = promote_compounds_to_notion("non-existent-campaign-id")
    assert result == [], f"Expected [], got {result}"
    print("✓ Returns empty list (Notion disabled)")
except Exception as e:
    print(f"✗ Test failed: {e}")
    raise

print("\n✅ BATCH B QUICK TESTS PASSED")

