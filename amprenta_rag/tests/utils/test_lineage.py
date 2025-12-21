from __future__ import annotations

from uuid import uuid4

import pytest

from amprenta_rag.utils import lineage


class FakeSig:
    def __init__(self, name="sig"):
        self.id = uuid4()
        self.name = name
        self.created_at = None


class FakeDataset:
    def __init__(self, name="ds", signatures=None):
        self.id = uuid4()
        self.name = name
        self.omics_type = "rna"
        self.created_at = None
        self.signatures = signatures or []


class FakeProtocol:
    def __init__(self, name="p", category="cat", version="1"):
        self.id = uuid4()
        self.name = name
        self.category = category
        self.version = version


class FakeProtocolLink:
    def __init__(self, protocol=None):
        self.protocol = protocol


class FakeProject:
    def __init__(self, name="proj"):
        self.id = uuid4()
        self.name = name
        self.created_at = None


class FakeExperiment:
    def __init__(self, name="exp", datasets=None, protocol_links=None, project=None):
        self.id = uuid4()
        self.name = name
        self.datasets = datasets or []
        self.protocol_links = protocol_links or []
        self.project_id = project.id if project else None
        self.project = project
        self.design_type = None
        self.created_at = None


class FakeQuery:
    def __init__(self, obj):
        self._obj = obj

    def filter(self, *args, **kwargs):
        return self

    def first(self):
        return self._obj

    def all(self):
        if isinstance(self._obj, list):
            return self._obj
        return [self._obj] if self._obj else []


class FakeSession:
    def __init__(self, exp=None, project=None, ds=None):
        self.exp = exp
        self.project = project
        self.ds = ds
        self.calls = 0

    def query(self, model):
        self.calls += 1
        if model.__name__ == "Experiment":
            return FakeQuery(self.exp)
        if model.__name__ == "Project":
            return FakeQuery(self.project)
        if model.__name__ == "Dataset":
            return FakeQuery(self.ds)
        return FakeQuery(None)


def test_get_experiment_lineage_builds_graph():
    sig = FakeSig("s1")
    ds = FakeDataset("ds1", signatures=[sig])
    proto = FakeProtocol("prot")
    exp = FakeExperiment("exp1", datasets=[ds], protocol_links=[FakeProtocolLink(proto)])
    project = FakeProject("proj1")
    exp.project_id = project.id
    db = FakeSession(exp=exp, project=project)

    graph = lineage.get_experiment_lineage(exp.id, db)
    assert len(graph["nodes"]) >= 3
    assert any(n["type"] == "dataset" for n in graph["nodes"])
    assert any(e["relationship"] == "has_dataset" for e in graph["edges"])


def test_get_entity_lineage_dispatch_and_unknown():
    exp = FakeExperiment("exp2")
    db = FakeSession(exp=exp)
    res_exp = lineage.get_entity_lineage("experiment", exp.id, db)
    assert res_exp["nodes"]

    res_unknown = lineage.get_entity_lineage("unknown", exp.id, db)
    assert res_unknown == {"nodes": [], "edges": []}

