from __future__ import annotations

import pytest

from amprenta_rag.utils import optimistic_lock


class FakeEntity:
    def __init__(self, version=1):
        self.version = version
        self.id = "123"
        self.name = "old"


class FakeDB:
    def __init__(self):
        self.committed = False

    def commit(self):
        self.committed = True


def test_check_version_matches():
    ent = FakeEntity(version=2)
    assert optimistic_lock.check_version(ent, 2) is True
    assert optimistic_lock.check_version(ent, 1) is False


def test_update_with_lock_success_increments_and_commits():
    ent = FakeEntity(version=1)
    db = FakeDB()
    updated = optimistic_lock.update_with_lock(ent, {"name": "new"}, expected_version=1, db=db)
    assert updated.name == "new"
    assert updated.version == 2
    assert db.committed is True


def test_update_with_lock_conflict():
    ent = FakeEntity(version=2)
    db = FakeDB()
    with pytest.raises(optimistic_lock.ConflictError):
        optimistic_lock.update_with_lock(ent, {"name": "x"}, expected_version=1, db=db)

