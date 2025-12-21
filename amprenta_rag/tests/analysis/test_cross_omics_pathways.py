from __future__ import annotations

from types import SimpleNamespace
from uuid import uuid4

from amprenta_rag.analysis import cross_omics_pathways as cop


class FakeFeature:
    def __init__(self, name: str | None, feature_type: str):
        self.name = name
        self.feature_type = feature_type
        self.id = uuid4()


class FakeProgram:
    def __init__(self, datasets):
        self.datasets = datasets


class FakeQuery:
    def __init__(self, data):
        self._data = data

    def join(self, *a, **k):
        return self

    def filter(self, *a, **k):
        return self

    def all(self):
        return list(self._data)

    def first(self):
        return self._data[0] if self._data else None


class FakeSession:
    def __init__(self, program=None, features=None):
        self.program = program
        self.features = features or []

    def query(self, model):
        if model is cop.Program:
            return FakeQuery([self.program] if self.program else [])
        return FakeQuery(self.features)

    def __enter__(self):
        return self

    def __exit__(self, *exc):
        return False


def test_collect_features_by_type(monkeypatch):
    feats = [
        FakeFeature("g1", "gene"),
        FakeFeature("p1", "protein"),
        FakeFeature("m1", "metabolite"),
        FakeFeature("l1", "lipid"),
        FakeFeature(None, "gene"),
    ]
    monkeypatch.setattr(cop, "db_session", lambda: FakeSession(features=feats))
    fs = cop._collect_features_for_datasets([uuid4()])
    assert fs.transcriptomics == ["g1"]
    assert fs.proteomics == ["p1"]
    assert fs.metabolomics == ["m1"]
    assert fs.lipidomics == ["l1"]


def test_analyze_pathway_convergence_maps_matches(monkeypatch):
    features = cop.OmicsFeatureSet(
        transcriptomics=["g1"],
        proteomics=["p1"],
        metabolomics=[],
        lipidomics=[],
    )

    fake_pathway = SimpleNamespace(pathway_id="P1", name="Pathway", source="Reactome")
    fake_res = SimpleNamespace(
        pathway=fake_pathway,
        adjusted_p_value=0.01,
        enrichment_ratio=2.0,
        matched_features=["g1", "p1"],
    )
    monkeypatch.setattr(cop, "perform_pathway_enrichment", lambda **k: [fake_res])

    conv = cop.analyze_pathway_convergence(features)
    assert conv[0].matched_by_omics["transcriptomics"] == ["g1"]
    assert conv[0].matched_by_omics["proteomics"] == ["p1"]


def test_get_cross_omics_enrichment_happy(monkeypatch):
    ds_id = uuid4()
    features = [FakeFeature("g1", "gene")]
    program = FakeProgram(datasets=[type("D", (), {"id": ds_id})()])

    fake_res = SimpleNamespace(
        pathway=SimpleNamespace(pathway_id="P1", name="Pathway", source="KEGG"),
        adjusted_p_value=0.05,
        enrichment_ratio=1.5,
        matched_features=["g1"],
    )

    monkeypatch.setattr(cop, "perform_pathway_enrichment", lambda **k: [fake_res])
    monkeypatch.setattr(cop, "db_session", lambda: FakeSession(program=program, features=features))

    res = cop.get_cross_omics_enrichment(uuid4())
    assert res.pathways[0].pathway_id == "P1"
    assert res.features.transcriptomics == ["g1"]


def test_dataclass_helpers():
    fs = cop.OmicsFeatureSet(["g1"], ["p1"], [], [])
    assert fs.all_features() == {"g1", "p1"}
    assert fs.feature_types() == {"gene", "protein"}

    cp = cop.ConvergentPathway(
        pathway_id="p1",
        name="Path",
        source="KEGG",
        adjusted_p_value=0.01,
        enrichment_ratio=2.0,
        matched_by_omics={"transcriptomics": ["g1"], "proteomics": [], "metabolomics": [], "lipidomics": []},
    )
    assert cp.asdict()["pathway_id"] == "p1"

    enr = cop.CrossOmicsEnrichmentResult(program_id=uuid4(), pathways=[cp], features=fs)
    assert "program_id" in enr.asdict()

