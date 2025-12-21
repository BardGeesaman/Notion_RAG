from __future__ import annotations

import os

import pytest

from amprenta_rag.client import helpers


class _FakeAnnotator:
    def __init__(self):
        self.calls = []

    def annotate(self, entity_id, text, annotation_type):
        self.calls.append((entity_id, text, annotation_type))
        return {"id": entity_id, "text": text, "type": annotation_type}


class _FakeClient:
    def __init__(self):
        self.experiments = _FakeAnnotator()
        self.datasets = _FakeAnnotator()
        self.signatures = _FakeAnnotator()
        self.compounds = _FakeAnnotator()


def test_get_context_from_url(monkeypatch):
    monkeypatch.setenv("EXPERIMENT_ID", "exp1")
    monkeypatch.setenv("DATASET_ID", "ds1")
    monkeypatch.setenv("COMPOUND_ID", "cmp1")
    ctx = helpers.get_context_from_url()
    assert ctx == {"experiment_id": "exp1", "dataset_id": "ds1", "compound_id": "cmp1"}
    # cleanup
    os.environ.pop("EXPERIMENT_ID")
    os.environ.pop("DATASET_ID")
    os.environ.pop("COMPOUND_ID")


def test_save_annotation_routes_to_entities():
    client = _FakeClient()
    helpers.save_annotation(client, "note", experiment_id="e1")
    helpers.save_annotation(client, "note", dataset_id="d1")
    helpers.save_annotation(client, "note", signature_id="s1")
    helpers.save_annotation(client, "note", compound_id="c1")

    assert client.experiments.calls[0][0] == "e1"
    assert client.datasets.calls[0][0] == "d1"
    assert client.signatures.calls[0][0] == "s1"
    assert client.compounds.calls[0][0] == "c1"


def test_save_annotation_requires_target():
    client = _FakeClient()
    with pytest.raises(ValueError):
        helpers.save_annotation(client, "note")

