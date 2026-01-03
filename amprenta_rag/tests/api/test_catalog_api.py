"""API tests for Data Catalog endpoints."""

import pytest
from uuid import uuid4
from unittest.mock import patch, MagicMock

from fastapi.testclient import TestClient

from amprenta_rag.api.main import app


@pytest.fixture
def client():
    """Test client fixture."""
    return TestClient(app)


@pytest.fixture
def mock_user():
    """Mock authenticated user."""
    user = MagicMock()
    user.id = uuid4()
    user.email = "test@test.com"
    user.is_admin = True
    return user


@pytest.fixture
def auth_override(mock_user):
    """Override auth dependency."""
    from amprenta_rag.api.dependencies import get_current_user
    app.dependency_overrides[get_current_user] = lambda: mock_user
    yield
    app.dependency_overrides.clear()


# ============================================================================
# CATALOG ENTRY TESTS
# ============================================================================

def test_list_catalog_entries(client):
    """Test GET /catalog/entries returns entries."""
    response = client.get("/api/v1/catalog/entries")
    assert response.status_code == 200
    data = response.json()
    assert isinstance(data, list)


def test_list_catalog_entries_filter_category(client):
    """Test category filter."""
    response = client.get("/api/v1/catalog/entries", params={"category": "Chemistry"})
    assert response.status_code == 200
    data = response.json()
    # Verify all returned entries have Chemistry category (if any)
    for entry in data:
        assert entry.get("category") == "Chemistry"


def test_list_catalog_entries_search(client):
    """Test search filter."""
    response = client.get("/api/v1/catalog/entries", params={"search": "Compound"})
    assert response.status_code == 200
    data = response.json()
    assert isinstance(data, list)


def test_get_catalog_entry(client):
    """Test GET /catalog/entries/{entity_type}."""
    # First get an existing entity type
    list_resp = client.get("/api/v1/catalog/entries")
    entries = list_resp.json()
    if entries:
        entity_type = entries[0]["entity_type"]
        response = client.get(f"/api/v1/catalog/entries/{entity_type}")
        assert response.status_code == 200
        data = response.json()
        assert data["entity_type"] == entity_type
        assert "columns" in data


def test_get_catalog_entry_not_found(client):
    """Test 404 for non-existent entity."""
    response = client.get("/api/v1/catalog/entries/NonExistentEntity")
    assert response.status_code == 404


def test_refresh_catalog(client, auth_override):
    """Test POST /catalog/refresh triggers discovery."""
    response = client.post("/api/v1/catalog/refresh")
    assert response.status_code == 200
    data = response.json()
    assert "entries_updated" in data
    assert "lineage_edges_created" in data
    assert isinstance(data["entries_updated"], int)


# ============================================================================
# COLUMN SEARCH TESTS
# ============================================================================

def test_search_columns(client):
    """Test GET /catalog/columns/search."""
    response = client.get("/api/v1/catalog/columns/search", params={"q": "id"})
    assert response.status_code == 200
    data = response.json()
    assert isinstance(data, list)


def test_search_columns_min_length(client):
    """Test search requires 2+ characters."""
    response = client.get("/api/v1/catalog/columns/search", params={"q": "a"})
    assert response.status_code == 422  # Validation error


def test_get_column_metadata(client):
    """Test GET /catalog/columns/{entity_type}/{column_name}."""
    # Get a column from an existing entity
    list_resp = client.get("/api/v1/catalog/entries")
    entries = list_resp.json()
    if entries:
        entity_type = entries[0]["entity_type"]
        detail_resp = client.get(f"/api/v1/catalog/entries/{entity_type}")
        if detail_resp.ok:
            columns = detail_resp.json().get("columns", [])
            if columns:
                col_name = columns[0]["column_name"]
                response = client.get(f"/api/v1/catalog/columns/{entity_type}/{col_name}")
                assert response.status_code == 200
                data = response.json()
                assert data["column_name"] == col_name


# ============================================================================
# GLOSSARY TESTS
# ============================================================================

def test_list_glossary_terms(client):
    """Test GET /catalog/glossary."""
    response = client.get("/api/v1/catalog/glossary")
    assert response.status_code == 200
    assert isinstance(response.json(), list)


def test_create_glossary_term(client, auth_override):
    """Test POST /catalog/glossary creates term."""
    unique_term = f"TestTerm_{uuid4().hex[:8]}"
    payload = {
        "term": unique_term,
        "definition": "A test definition",
        "category": "Chemistry",
        "synonyms": ["syn1", "syn2"]
    }
    response = client.post("/api/v1/catalog/glossary", json=payload)
    assert response.status_code == 201
    data = response.json()
    assert data["term"] == unique_term
    assert data["definition"] == "A test definition"
    assert data["category"] == "Chemistry"
    
    # Cleanup
    cleanup_resp = client.delete(f"/api/v1/catalog/glossary/{data['id']}")
    assert cleanup_resp.status_code == 204


def test_create_glossary_term_duplicate(client, auth_override):
    """Test duplicate term returns error."""
    unique_term = f"DupeTerm_{uuid4().hex[:8]}"
    payload = {"term": unique_term, "definition": "First"}
    
    # Create first
    resp1 = client.post("/api/v1/catalog/glossary", json=payload)
    assert resp1.status_code == 201
    term_id = resp1.json()["id"]
    
    # Try duplicate
    resp2 = client.post("/api/v1/catalog/glossary", json=payload)
    assert resp2.status_code == 400
    assert "already exists" in resp2.json()["detail"]
    
    # Cleanup
    client.delete(f"/api/v1/catalog/glossary/{term_id}")


def test_update_glossary_term(client, auth_override):
    """Test PUT /catalog/glossary/{term_id}."""
    # Create term
    unique_term = f"UpdateTerm_{uuid4().hex[:8]}"
    create_resp = client.post("/api/v1/catalog/glossary", json={
        "term": unique_term,
        "definition": "Original"
    })
    assert create_resp.status_code == 201
    term_id = create_resp.json()["id"]
    
    # Update
    update_resp = client.put(f"/api/v1/catalog/glossary/{term_id}", json={
        "definition": "Updated definition"
    })
    assert update_resp.status_code == 200
    assert update_resp.json()["definition"] == "Updated definition"
    
    # Cleanup
    client.delete(f"/api/v1/catalog/glossary/{term_id}")


def test_delete_glossary_term(client, auth_override):
    """Test DELETE /catalog/glossary/{term_id}."""
    # Create term
    unique_term = f"DeleteTerm_{uuid4().hex[:8]}"
    create_resp = client.post("/api/v1/catalog/glossary", json={
        "term": unique_term,
        "definition": "To delete"
    })
    assert create_resp.status_code == 201
    term_id = create_resp.json()["id"]
    
    # Delete
    delete_resp = client.delete(f"/api/v1/catalog/glossary/{term_id}")
    assert delete_resp.status_code == 204
    
    # Verify gone
    get_resp = client.get(f"/api/v1/catalog/glossary/{term_id}")
    assert get_resp.status_code == 404


def test_get_glossary_term(client, auth_override):
    """Test GET /catalog/glossary/{term_id}."""
    # Create term
    unique_term = f"GetTerm_{uuid4().hex[:8]}"
    create_resp = client.post("/api/v1/catalog/glossary", json={
        "term": unique_term,
        "definition": "To get"
    })
    assert create_resp.status_code == 201
    term_id = create_resp.json()["id"]
    
    # Get
    get_resp = client.get(f"/api/v1/catalog/glossary/{term_id}")
    assert get_resp.status_code == 200
    data = get_resp.json()
    assert data["term"] == unique_term
    
    # Cleanup
    client.delete(f"/api/v1/catalog/glossary/{term_id}")


# ============================================================================
# LINEAGE TESTS
# ============================================================================

def test_get_lineage_graph(client):
    """Test GET /catalog/lineage/{entity_type}/{entity_id}."""
    # Use a catalog entry ID for testing
    list_resp = client.get("/api/v1/catalog/entries")
    entries = list_resp.json()
    if entries:
        entry = entries[0]
        response = client.get(
            f"/api/v1/catalog/lineage/{entry['entity_type']}/{entry['id']}",
            params={"depth": 2, "direction": "both"}
        )
        assert response.status_code == 200
        data = response.json()
        assert "nodes" in data
        assert "edges" in data
        assert "center_entity" in data


def test_add_lineage_edge(client, auth_override):
    """Test POST /catalog/lineage creates edge."""
    payload = {
        "source_type": "TestSource",
        "source_id": str(uuid4()),
        "target_type": "TestTarget", 
        "target_id": str(uuid4()),
        "relationship_type": "derived_from",
        "transformation": "test_transform"
    }
    response = client.post("/api/v1/catalog/lineage", json=payload)
    # May return 201 or 400 if duplicate
    assert response.status_code in [201, 400]
    if response.status_code == 201:
        data = response.json()
        assert data["created"] is True
        assert "id" in data


# ============================================================================
# INTEGRATION TESTS
# ============================================================================

def test_catalog_endpoints_exist(client):
    """Test that all expected catalog endpoints are registered."""
    # Test that endpoints return proper HTTP status (not 404 for missing routes)
    endpoints_to_test = [
        ("GET", "/api/v1/catalog/entries"),
        ("GET", "/api/v1/catalog/glossary"),
        ("GET", "/api/v1/catalog/columns/search?q=id"),
    ]
    
    for method, path in endpoints_to_test:
        if method == "GET":
            response = client.get(path)
            # Should not be 404 (route not found)
            assert response.status_code != 404


def test_catalog_api_structure(client):
    """Test API returns proper structure."""
    # Test entries structure
    resp = client.get("/api/v1/catalog/entries")
    if resp.status_code == 200:
        entries = resp.json()
        if entries:
            entry = entries[0]
            required_fields = ["id", "entity_type", "table_name", "display_name", "category"]
            for field in required_fields:
                assert field in entry
    
    # Test glossary structure
    resp = client.get("/api/v1/catalog/glossary")
    if resp.status_code == 200:
        terms = resp.json()
        if terms:
            term = terms[0]
            required_fields = ["id", "term", "definition"]
            for field in required_fields:
                assert field in term


def test_api_error_handling(client):
    """Test API error responses."""
    # Test 404 for non-existent resources
    response = client.get("/api/v1/catalog/entries/NonExistentEntity")
    assert response.status_code == 404
    assert "detail" in response.json()
    
    # Test 422 for validation errors
    response = client.get("/api/v1/catalog/columns/search", params={"q": "a"})
    assert response.status_code == 422