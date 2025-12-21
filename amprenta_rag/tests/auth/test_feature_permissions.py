from __future__ import annotations

import sys
from types import ModuleType
import pytest

# Mock scripts.run_dashboard to provide ALL_PAGES/ADMIN_PAGES constants
fake_dash = ModuleType("scripts.run_dashboard")
fake_dash.ALL_PAGES = ["Home", "Settings", "Admin"]
fake_dash.ADMIN_PAGES = ["Admin"]
sys.modules["scripts.run_dashboard"] = fake_dash

from amprenta_rag.auth import feature_permissions as fp
from amprenta_rag.database.models import FeaturePermission


class FakePermission:
    def __init__(self, role, page_name, is_visible=True):
        self.role = role
        self.page_name = page_name
        self.is_visible = is_visible


class FakeQuery:
    def __init__(self, results=None):
        self.results = results or []

    def filter(self, *args, **kwargs):
        # Basic filtering simulation
        filtered = []
        for item in self.results:
            match = True
            # This is a very simplified filter mock; relying on test setup to pre-filter
            # or pass relevant objects.
            # In a real test, we might filter by attributes, but here we trust the setup.
            filtered.append(item)
        return FakeQuery(filtered)

    def first(self):
        return self.results[0] if self.results else None

    def all(self):
        return self.results

    def count(self):
        return len(self.results)


class FakeSession:
    def __init__(self, permissions=None):
        self.permissions = permissions or []
        self.added = []
        self.committed = False
        self.rolled_back = False

    def query(self, model):
        # Return a query object that wraps our permissions list
        # For get_visible_pages, we need to handle specific filters
        return FakeQuery(self.permissions)

    def add(self, obj):
        self.added.append(obj)
        self.permissions.append(obj)

    def commit(self):
        self.committed = True

    def rollback(self):
        self.rolled_back = True


def test_get_visible_pages_no_permissions(monkeypatch):
    # Setup: Empty DB returns no permissions
    # Logic: Should fallback to ALL_PAGES
    session = FakeSession([])
    
    # We need to mock the query behavior specifically for the "has_permissions" check
    # The code does db.query(...).filter(...).first()
    # FakeSession.query returns a FakeQuery with ALL items.
    # We need it to return None for the first check.
    
    def fake_query(model):
        return FakeQuery([]) # Return empty query
        
    session.query = fake_query
    
    pages = fp.get_visible_pages("admin", session)
    assert pages == ["Home", "Settings", "Admin"]


def test_get_visible_pages_with_permissions(monkeypatch):
    # Setup: DB has permissions
    p1 = FakePermission("admin", "Home", True)
    p2 = FakePermission("admin", "Settings", False) # Not visible
    
    session = FakeSession([p1, p2])
    
    # Mock query to handle the two-step logic:
    # 1. Check if permissions exist (first() -> not None)
    # 2. Get visible permissions (all() -> list)
    
    def fake_query(model):
        class SmartQuery:
            def filter(self, *args):
                # If checking visibility (2nd query), return visible only
                # The code checks FeaturePermission.is_visible in filter args
                # Since we can't easily inspect SQLAlchemy binary expressions in a mock,
                # we'll assume the second call is for filtered results.
                
                # Hacky: Assume if we call all(), we want visible ones
                return self
                
            def first(self):
                # Used for "has_permissions" check
                return p1 
                
            def all(self):
                # Return only visible permissions for this role
                return [p1]
                
        return SmartQuery()

    session.query = fake_query
    
    pages = fp.get_visible_pages("admin", session)
    assert pages == ["Home"]


def test_set_page_visibility_update_existing():
    p1 = FakePermission("admin", "Home", True)
    session = FakeSession([p1])
    
    # Mock query to find existing
    def fake_query(model):
        class FoundQuery:
            def filter(self, *args): return self
            def first(self): return p1
        return FoundQuery()
    session.query = fake_query
    
    success = fp.set_page_visibility("admin", "Home", False, session)
    assert success is True
    assert p1.is_visible is False
    assert session.committed is True


def test_set_page_visibility_create_new():
    session = FakeSession([])
    
    # Mock query to NOT find existing
    def fake_query(model):
        class NotFoundQuery:
            def filter(self, *args): return self
            def first(self): return None
        return NotFoundQuery()
    session.query = fake_query
    
    success = fp.set_page_visibility("admin", "Home", True, session)
    assert success is True
    assert len(session.added) == 1
    assert session.added[0].role == "admin"
    assert session.added[0].page_name == "Home"
    assert session.committed is True


def test_set_page_visibility_error_handling():
    session = FakeSession()
    session.commit = lambda: (_ for _ in ()).throw(RuntimeError("Db Error"))
    
    # Mock query to find existing
    def fake_query(model):
        class FoundQuery:
            def filter(self, *args): return self
            def first(self): return FakePermission("r", "p")
        return FoundQuery()
    session.query = fake_query
    
    success = fp.set_page_visibility("admin", "Home", False, session)
    assert success is False
    assert session.rolled_back is True


def test_get_all_permissions():
    p1 = FakePermission("admin", "Home", True)
    p2 = FakePermission("viewer", "Home", False)
    session = FakeSession([p1, p2])
    
    res = fp.get_all_permissions(session)
    assert res["admin"]["Home"] is True
    assert res["viewer"]["Home"] is False


def test_initialize_default_permissions_already_exists():
    p1 = FakePermission("admin", "Home", True)
    session = FakeSession([p1])
    fp.initialize_default_permissions(session)
    # Should perform no additions if count > 0
    assert len(session.added) == 0


def test_initialize_default_permissions_creates_defaults():
    session = FakeSession([]) # count = 0
    fp.initialize_default_permissions(session)
    
    assert len(session.added) > 0
    assert session.committed is True
    
    # Check Admin gets all
    admin_perms = [p for p in session.added if p.role == "admin"]
    assert len(admin_perms) == len(fake_dash.ALL_PAGES)
    
    # Check Researcher gets non-admin
    researcher_perms = [p for p in session.added if p.role == "researcher"]
    expected_researcher = len(fake_dash.ALL_PAGES) - len(fake_dash.ADMIN_PAGES)
    assert len(researcher_perms) == expected_researcher

