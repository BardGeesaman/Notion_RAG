from __future__ import annotations

from typing import Any, List

from amprenta_rag.utils import health


class FakeQuery:
    def __init__(self, data: List[Any], active_only: bool = False):
        self._data = data
        self._active_only = active_only

    def filter(self, predicate: Any):
        if self._active_only:
            filtered = [u for u in self._data if getattr(u, "is_active", False)]
            return FakeQuery(filtered)
        return self

    def count(self) -> int:
        return len(self._data)


class FakeUser:
    def __init__(self, is_active: bool):
        self.is_active = is_active


class FakeEntity:
    pass


class FakeSession:
    def __init__(self):
        self.users: List[FakeUser] = []
        self.experiments: List[FakeEntity] = []
        self.compounds: List[FakeEntity] = []
        self.signatures: List[FakeEntity] = []
        self.datasets: List[FakeEntity] = []

    def query(self, model: Any) -> FakeQuery:
        if model is health.User:
            return FakeQuery(self.users, active_only=True)
        if model is health.Experiment:
            return FakeQuery(self.experiments)
        if model is health.Compound:
            return FakeQuery(self.compounds)
        if model is health.Signature:
            return FakeQuery(self.signatures)
        if model is health.Dataset:
            return FakeQuery(self.datasets)
        return FakeQuery([])


def test_get_db_stats_counts_active_users_only():
    db = FakeSession()
    db.users.extend([FakeUser(True), FakeUser(False)])
    db.experiments.append(FakeEntity())
    stats = health.get_db_stats(db)
    assert stats["users"] == 1
    assert stats["experiments"] == 1


def test_get_system_info_has_basic_keys():
    info = health.get_system_info()
    assert "python_version" in info and "platform" in info
    assert isinstance(info["python_version"], str)

