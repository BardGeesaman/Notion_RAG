"""Batch A: Query Layer Tests"""

import logging

logging.basicConfig(level=logging.DEBUG)

# Test 1: Imports
print("Test 1: Imports...")
try:
    from amprenta_rag.query.rag.chunk_collection import collect_chunks
    from amprenta_rag.rag.hybrid_chunk_collection import collect_hybrid_chunks
    from amprenta_rag.query.cross_omics.helpers import fetch_notion_page
    from amprenta_rag.query.cross_omics.context_extraction import extract_aggregated_context
    print("✓ All imports OK")
except Exception as e:
    print(f"✗ Import failed: {e}")
    raise

# Test 2: collect_chunks
print("\nTest 2: collect_chunks...")
try:
    result = collect_chunks([])
    assert result == [], f"Expected [], got {result}"
    print("✓ Returns empty list")
except Exception as e:
    print(f"✗ collect_chunks failed: {e}")
    raise

# Test 3: fetch_notion_page
print("\nTest 3: fetch_notion_page...")
try:
    result = fetch_notion_page('test-id')
    assert isinstance(result, dict), f"Expected dict, got {type(result)}"
    assert result.get('id') == 'test-id', f"Expected id='test-id', got {result.get('id')}"
    assert result.get('properties') == {}, f"Expected empty properties, got {result.get('properties')}"
    print("✓ Returns empty dict")
except Exception as e:
    print(f"✗ fetch_notion_page failed: {e}")
    raise

# Test 4: extract_aggregated_context
print("\nTest 4: extract_aggregated_context...")
try:
    result = extract_aggregated_context(['test-id'])
    assert 'diseases' in result and result['diseases'] == []
    assert 'page_count' in result and result['page_count'] == 0
    print("✓ Returns empty context")
except Exception as e:
    print(f"✗ extract_aggregated_context failed: {e}")
    raise

print("\n✅ BATCH A TESTS PASSED")

