from __future__ import annotations

from types import SimpleNamespace
from uuid import uuid4

from amprenta_rag.analysis import cross_omics_pathways as cop


class FakeFeature:
    def __init__(self, name, feature_type):
        self.name = name
        self.feature_type = feature_type
        self.id = uuid4()


class FakeDataset:
    def __init__(self, id_, features):
        self.id = id_
        self.features = features


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
        self.calls = 0

    def query(self, model):
        # First call for Program in get_cross_omics_enrichment
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
from __future__ import annotations

from contextlib import contextmanager
from types import SimpleNamespace
from uuid import uuid4
from unittest.mock import Mock

import pytest
from fastapi.testclient import TestClient


@pytest.mark.unit
def test_omicsfeatureset_dataclass_helpers():
    from amprenta_rag.analysis.cross_omics_pathways import OmicsFeatureSet

    fs = OmicsFeatureSet(
        transcriptomics=["G1", "G2"],
        proteomics=["P1"],
        metabolomics=["M1"],
        lipidomics=[],
    )
    assert fs.all_features() == {"G1", "G2", "P1", "M1"}
    assert fs.feature_types() == {"gene", "protein", "metabolite"}
    assert fs.asdict() == {
        "transcriptomics": ["G1", "G2"],
        "proteomics": ["P1"],
        "metabolomics": ["M1"],
        "lipidomics": [],
    }


@pytest.mark.unit
def test_convergentpathway_dataclass_asdict():
    from amprenta_rag.analysis.cross_omics_pathways import ConvergentPathway

    cp = ConvergentPathway(
        pathway_id="p1",
        name="Pathway 1",
        source="KEGG",
        adjusted_p_value=0.01,
        enrichment_ratio=2.5,
        matched_by_omics={"transcriptomics": ["G1"], "proteomics": [], "metabolomics": [], "lipidomics": []},
    )
    assert cp.asdict()["pathway_id"] == "p1"
    assert cp.asdict()["matched_by_omics"]["transcriptomics"] == ["G1"]


@pytest.mark.unit
def test_combine_omics_features_groups_features_by_type(monkeypatch):
    from amprenta_rag.analysis import cross_omics_pathways as cop

    # Fake Feature-like objects
    feats = [
        SimpleNamespace(name="G1", feature_type="gene"),
        SimpleNamespace(name="P1", feature_type="protein"),
        SimpleNamespace(name="M1", feature_type="metabolite"),
        SimpleNamespace(name="L1", feature_type="lipid"),
        SimpleNamespace(name=None, feature_type="gene"),  # skipped
        SimpleNamespace(name="X", feature_type="unknown"),  # skipped
    ]

    q = SimpleNamespace(
        join=lambda *args, **kwargs: q,
        filter=lambda *args, **kwargs: q,
        all=lambda: feats,
    )
    db = SimpleNamespace(query=lambda _model: q)

    @contextmanager
    def fake_db_session():
        yield db

    monkeypatch.setattr(cop, "db_session", fake_db_session)

    fs = cop.combine_omics_features([uuid4(), uuid4()])
    assert fs.transcriptomics == ["G1"]
    assert fs.proteomics == ["P1"]
    assert fs.metabolomics == ["M1"]
    assert fs.lipidomics == ["L1"]


@pytest.mark.unit
def test_analyze_pathway_convergence_returns_empty_on_no_features_or_types(monkeypatch):
    from amprenta_rag.analysis.cross_omics_pathways import OmicsFeatureSet, analyze_pathway_convergence

    fs = OmicsFeatureSet(transcriptomics=[], proteomics=[], metabolomics=[], lipidomics=[])
    assert analyze_pathway_convergence(fs) == []


@pytest.mark.unit
def test_analyze_pathway_convergence_builds_matched_by_omics(monkeypatch):
    from amprenta_rag.analysis import cross_omics_pathways as cop

    fs = cop.OmicsFeatureSet(
        transcriptomics=["G1", "Xshared"],
        proteomics=["P1", "Xshared"],
        metabolomics=["M1"],
        lipidomics=[],
    )

    class FakePathway:
        def __init__(self, pathway_id: str, name: str, source: str):
            self.pathway_id = pathway_id
            self.name = name
            self.source = source

    fake_res = SimpleNamespace(
        pathway=FakePathway("PWY-1", "Demo Pathway", "KEGG"),
        adjusted_p_value=0.02,
        enrichment_ratio=1.7,
        matched_features=["G1", "Xshared", "M1", "NOT_IN_INPUT"],
    )

    monkeypatch.setattr(cop, "perform_pathway_enrichment", Mock(return_value=[fake_res]))

    out = cop.analyze_pathway_convergence(fs)
    assert len(out) == 1
    p = out[0]
    assert p.pathway_id == "PWY-1"
    assert p.name == "Demo Pathway"
    assert p.source == "KEGG"
    assert p.adjusted_p_value == 0.02
    assert p.enrichment_ratio == 1.7
    assert p.matched_by_omics["transcriptomics"] == ["G1", "Xshared"]
    assert p.matched_by_omics["proteomics"] == ["Xshared"]
    assert p.matched_by_omics["metabolomics"] == ["M1"]
    assert p.matched_by_omics["lipidomics"] == []


@pytest.fixture()
def api_client():
    from amprenta_rag.api.main import app

    return TestClient(app)


@pytest.mark.api
def test_post_pathways_enrich_200(api_client, monkeypatch):
    from amprenta_rag.api.routers import pathways as pathways_router

    dataset_id = uuid4()

    monkeypatch.setattr(pathways_router, "_load_dataset_features", lambda _id: ({"G1", "G2"}, {"gene"}))

    class FakePathway:
        def __init__(self):
            self.pathway_id = "PWY-1"
            self.name = "Demo"
            self.source = "KEGG"

    fake_res = SimpleNamespace(
        pathway=FakePathway(),
        adjusted_p_value=0.01,
        enrichment_ratio=2.0,
        matched_features=["G1"],
    )
    monkeypatch.setattr(pathways_router, "perform_pathway_enrichment", Mock(return_value=[fake_res]))

    resp = api_client.post(f"/api/v1/pathways/enrich?dataset_id={dataset_id}")
    assert resp.status_code == 200
    assert resp.json() == [
        {
            "pathway_id": "PWY-1",
            "name": "Demo",
            "source": "KEGG",
            "adjusted_p_value": 0.01,
            "enrichment_ratio": 2.0,
            "matched_by_omics": {"all": ["G1"]},
        }
    ]


@pytest.mark.api
def test_post_pathways_enrich_404_when_no_features(api_client, monkeypatch):
    from amprenta_rag.api.routers import pathways as pathways_router

    dataset_id = uuid4()
    monkeypatch.setattr(pathways_router, "_load_dataset_features", lambda _id: (set(), set()))
    resp = api_client.post(f"/api/v1/pathways/enrich?dataset_id={dataset_id}")
    assert resp.status_code == 404
    assert resp.json()["detail"] == "No features for dataset"


@pytest.mark.api
def test_get_program_pathway_analysis_200(api_client, monkeypatch):
    from amprenta_rag.api.routers import pathways as pathways_router
    from amprenta_rag.analysis.cross_omics_pathways import CrossOmicsEnrichmentResult, ConvergentPathway, OmicsFeatureSet

    program_id = uuid4()
    result = CrossOmicsEnrichmentResult(
        program_id=program_id,
        pathways=[
            ConvergentPathway(
                pathway_id="PWY-1",
                name="Demo",
                source="KEGG",
                adjusted_p_value=0.01,
                enrichment_ratio=2.0,
                matched_by_omics={"transcriptomics": ["G1"], "proteomics": [], "metabolomics": [], "lipidomics": []},
            )
        ],
        features=OmicsFeatureSet(transcriptomics=["G1"], proteomics=[], metabolomics=[], lipidomics=[]),
    )

    monkeypatch.setattr(pathways_router, "get_cross_omics_enrichment", Mock(return_value=result))

    resp = api_client.get(f"/api/v1/programs/{program_id}/pathway-analysis")
    assert resp.status_code == 200
    body = resp.json()
    assert body["program_id"] == str(program_id)
    assert body["features"]["transcriptomics"] == ["G1"]
    assert body["pathways"][0]["pathway_id"] == "PWY-1"


@pytest.mark.api
def test_get_pathway_features_returns_combined_features(api_client, monkeypatch):
    from amprenta_rag.api.routers import pathways as pathways_router
    from amprenta_rag.analysis.cross_omics_pathways import OmicsFeatureSet

    ds1 = uuid4()
    ds2 = uuid4()
    monkeypatch.setattr(
        pathways_router,
        "combine_omics_features",
        Mock(return_value=OmicsFeatureSet(transcriptomics=["G1"], proteomics=["P1"], metabolomics=[], lipidomics=[])),
    )

    # Note: FastAPI treats `dataset_ids: List[UUID]` as a body param for GET (complex type),
    # so we must pass the UUID list as JSON body rather than query params.
    resp = api_client.request(
        "GET",
        "/api/v1/pathways/PWY-123/features",
        json=[str(ds1), str(ds2)],
    )
    assert resp.status_code == 200
    body = resp.json()
    assert body["pathway_id"] == "PWY-123"
    assert body["features_by_omics"]["transcriptomics"] == ["G1"]
    assert body["features_by_omics"]["proteomics"] == ["P1"]


