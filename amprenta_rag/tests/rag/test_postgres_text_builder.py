from __future__ import annotations

from datetime import datetime
from uuid import uuid4

from amprenta_rag.rag import postgres_text_builder as ptb


class _FakeQuery:
    def __init__(self, obj):
        self.obj = obj

    def options(self, *args, **kwargs):
        return self

    def filter(self, *args, **kwargs):
        return self

    def first(self):
        return self.obj


class _FakeSession:
    def __init__(self, obj):
        self.obj = obj

    def query(self, model):
        return _FakeQuery(self.obj)


class _FakeObj:
    def __init__(self, **kwargs):
        for k, v in kwargs.items():
            setattr(self, k, v)


def test_build_program_text_happy():
    prog = _FakeObj(
        id=uuid4(),
        name="Prog",
        description="Desc",
        disease=["DiseaseA"],
        experiments=[
            _FakeObj(name="Exp1", type="RNA", matrix=["blood"]),
            _FakeObj(name="Exp2", type=None, matrix=[]),
        ],
        datasets=[
            _FakeObj(name="DS1", description="dataset desc", omics_type="gene"),
            _FakeObj(name="DS2", description=None, omics_type="protein"),
        ],
        signatures=[_FakeObj(name="Sig1", modalities=["gene"])],
        created_at=datetime(2024, 1, 1),
        updated_at=datetime(2024, 1, 2),
    )
    text = ptb.build_program_text(_FakeSession(prog), prog.id)
    assert "# Program: Prog" in text
    assert "Disease Areas" in text
    assert "Experiments (2)" in text
    assert "Datasets (2)" in text
    assert "Signatures (1)" in text
    assert "Program ID" in text


def test_build_program_text_not_found():
    text = ptb.build_program_text(_FakeSession(None), uuid4())
    assert text == ""


def test_build_dataset_text_happy():
    dataset = _FakeObj(
        id=uuid4(),
        name="DS",
        omics_type="metabolite",
        description="desc",
        organism=["human"],
        sample_type=["plasma"],
        disease=["covid"],
        methods="methods",
        summary="summary",
        results="results",
        conclusions="conclusions",
        features=[
            _FakeObj(name="f1", feature_type="lipid"),
            _FakeObj(name="f2", feature_type="lipid"),
        ],
        programs=[_FakeObj(name="P1")],
        experiments=[_FakeObj(name="E1")],
        signatures=[_FakeObj(name="S1")],
        data_origin="external",
        dataset_source_type="repo",
        created_at=datetime(2024, 2, 1),
        updated_at=datetime(2024, 2, 2),
        file_paths=["/tmp/file1"],
    )
    text = ptb.build_dataset_text(_FakeSession(dataset), dataset.id)
    assert "# Dataset: DS" in text
    assert "**Omics Type**: metabolite" in text
    assert "Features (2 total)" in text
    assert "Programs" in text
    assert "Metadata" in text


def test_build_dataset_text_not_found():
    text = ptb.build_dataset_text(_FakeSession(None), uuid4())
    assert text == ""


def test_build_signature_text_happy():
    sig = _FakeObj(
        id=uuid4(),
        name="Sig",
        description="desc",
        modalities=["gene", "protein"],
        short_id="S1",
        biomarker_role=["diagnostic"],
        phenotype_axes=["ax1"],
        data_ownership="internal",
        components=[
            _FakeObj(feature_type="gene", feature_name="G1", direction="up", weight=1.2),
            _FakeObj(feature_type="protein", feature_name="P1", direction=None, weight=None),
        ],
        programs=[_FakeObj(name="P1")],
        datasets=[_FakeObj(name="DS1", omics_type="gene")],
        created_at=datetime(2024, 3, 1),
        updated_at=datetime(2024, 3, 2),
    )
    text = ptb.build_signature_text(_FakeSession(sig), sig.id)
    assert "# Signature: Sig" in text
    assert "Components (2 total)" in text
    assert "Modalities" in text
    assert "Matched Datasets" in text


def test_build_signature_text_not_found():
    text = ptb.build_signature_text(_FakeSession(None), uuid4())
    assert text == ""


def test_build_feature_text_happy():
    feat = _FakeObj(
        id=uuid4(),
        name="Feat",
        feature_type="lipid",
        normalized_name="feat_norm",
        aliases=["fA", "fB"],
        external_ids={"k": "v"},
        datasets=[_FakeObj(name="DS1", omics_type="lipid")],
        signatures=[_FakeObj(name="Sig1")],
        created_at=datetime(2024, 4, 1),
        updated_at=datetime(2024, 4, 2),
    )
    text = ptb.build_feature_text(_FakeSession(feat), feat.id)
    assert "# Feature: Feat" in text
    assert "Normalized Name" in text
    assert "Datasets" in text
    assert "Signatures" in text


def test_build_feature_text_not_found():
    text = ptb.build_feature_text(_FakeSession(None), uuid4())
    assert text == ""

