from __future__ import annotations

import json
from uuid import uuid4

import pandas as pd
import pytest

from amprenta_rag.utils import data_export


class FakeExperiment:
    def __init__(self, name: str):
        self.id = uuid4()
        self.name = name
        self.type = "type"
        self.description = "desc"
        self.disease = ["covid"]
        self.matrix = None
        self.model_systems = None
        self.targets = None
        self.modality = None
        self.stage = None
        self.biomarker_role = None
        self.treatment_arms = None
        self.design_type = None
        self.design_metadata = {"k": 1}
        self.sample_groups = {"g": []}
        self.created_at = None
        self.updated_at = None


class FakeCompound:
    def __init__(self, compound_id: str):
        self.compound_id = compound_id
        self.smiles = "CCO"
        self.canonical_smiles = "CCO"
        self.inchi_key = "INCHI"
        self.molecular_formula = "C2H6O"
        self.molecular_weight = 46.0
        self.logp = 0.1
        self.hbd_count = 1
        self.hba_count = 1
        self.rotatable_bonds = 0
        self.external_ids = {"foo": "bar"}
        self.created_at = None
        self.updated_at = None


class FakeSignature:
    def __init__(self, name: str):
        self.id = uuid4()
        self.name = name
        self.short_id = "sig1"
        self.description = "sig"
        self.modalities = ["rna"]
        self.biomarker_role = None
        self.phenotype_axes = None
        self.data_ownership = None
        self.created_at = None
        self.updated_at = None


class FakeQuery:
    def __init__(self, data):
        self._data = list(data)

    def filter(self, *args, **kwargs):
        return self

    def all(self):
        return list(self._data)


class FakeSession:
    def __init__(self):
        self.experiments = []
        self.compounds = []
        self.signatures = []

    def query(self, model):
        if model is data_export.Experiment:
            return FakeQuery(self.experiments)
        if model is data_export.Compound:
            return FakeQuery(self.compounds)
        if model is data_export.Signature:
            return FakeQuery(self.signatures)
        return FakeQuery([])


def test_export_to_csv_and_json_roundtrip():
    df = pd.DataFrame([{"a": 1, "b": "x"}])
    csv_bytes = data_export.export_to_csv(df)
    json_bytes = data_export.export_to_json(df)
    assert b"a,b" in csv_bytes
    assert json.loads(json_bytes.decode()) == [{"a": 1, "b": "x"}]


def test_export_experiments_formats(monkeypatch):
    db = FakeSession()
    db.experiments.append(FakeExperiment("exp1"))
    monkeypatch.setattr(data_export, "export_to_excel", lambda df: b"excel-bytes")

    out_json = data_export.export_experiments([], "json", db)
    assert out_json.startswith(b"[")  # empty path

    out_excel = data_export.export_experiments([str(db.experiments[0].id)], "excel", db)
    assert out_excel == b"excel-bytes"

    with pytest.raises(ValueError):
        data_export.export_experiments([], "bogus", db)


def test_export_compounds_and_signatures(monkeypatch):
    db = FakeSession()
    db.compounds.append(FakeCompound("C1"))
    db.signatures.append(FakeSignature("sig"))
    monkeypatch.setattr(data_export, "export_to_excel", lambda df: b"excel")

    comp_bytes = data_export.export_compounds(["C1"], "json", db)
    assert b"C1" in comp_bytes

    sig_bytes = data_export.export_signatures([str(db.signatures[0].id)], "csv", db)
    assert b"sig1" in sig_bytes

