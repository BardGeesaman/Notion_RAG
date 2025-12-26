from __future__ import annotations

from types import SimpleNamespace

import pytest
from fastapi.testclient import TestClient

from amprenta_rag.api.main import app


pytest.importorskip("fastapi")

client = TestClient(app)


def test_structure_endpoint_returns_nodes_edges(monkeypatch):
    from amprenta_rag.api.routers import pathway_maps as pm_router

    fake_struct = SimpleNamespace(
        pathway_id="hsa00010",
        name="Glycolysis",
        organism="hsa",
        nodes=[SimpleNamespace(id="1", name="TP53", type="gene", x=0.1, y=0.2, kegg_ids=["hsa:7157"])],
        edges=[SimpleNamespace(source="1", target="1", type="PPrel", subtype="activation")],
    )
    monkeypatch.setattr(pm_router, "get_pathway_structure", lambda pathway_id: fake_struct)
    monkeypatch.setattr(pm_router, "edge_style", lambda e: {"style": "solid", "color": "#2E8B57"})

    resp = client.get("/api/pathway-maps/structure/hsa00010")
    assert resp.status_code == 200
    data = resp.json()
    assert data["pathway_id"] == "hsa00010"
    assert len(data["nodes"]) == 1
    assert len(data["edges"]) == 1
    assert data["edges"][0]["color"] == "#2E8B57"


def test_overlay_endpoint_applies_colors(monkeypatch):
    from amprenta_rag.api.routers import pathway_maps as pm_router

    fake_struct = SimpleNamespace(pathway_id="hsa00010", name="x", organism="hsa", nodes=[], edges=[])
    monkeypatch.setattr(pm_router, "get_pathway_structure", lambda pathway_id: fake_struct)

    monkeypatch.setattr(
        pm_router,
        "compute_overlay",
        lambda structure, expression_data, colormap="RdBu_r", vmin=-2.0, vmax=2.0: [
            SimpleNamespace(node_id="1", gene_symbol="TP53", value=1.0, color="#ff0000", label="TP53: +1.00")
        ],
    )

    resp = client.post(
        "/api/pathway-maps/overlay",
        json={"pathway_id": "hsa00010", "expression_data": {"TP53": 1.0}, "colormap": "RdBu_r", "vmin": -2, "vmax": 2},
    )
    assert resp.status_code == 200
    data = resp.json()
    assert data["pathway_id"] == "hsa00010"
    assert data["overlays"][0]["color"] == "#ff0000"


def test_search_endpoint_returns_results(monkeypatch):
    from amprenta_rag.api.routers import pathway_maps as pm_router

    class Resp:
        status_code = 200
        text = "path:hsa00010\tGlycolysis / Gluconeogenesis - Homo sapiens (human)\npath:hsa00020\tCitrate cycle\n"

    monkeypatch.setattr(pm_router, "_rate_limit_kegg_search", lambda: None)
    monkeypatch.setattr(pm_router.requests, "get", lambda url, timeout=15: Resp())

    resp = client.get("/api/pathway-maps/search", params={"query": "glycolysis"})
    assert resp.status_code == 200
    data = resp.json()
    assert data[0]["pathway_id"] == "hsa00010"


def test_enriched_endpoint_returns_top_pathways(monkeypatch):
    from amprenta_rag.api.routers import pathway_maps as pm_router

    monkeypatch.setattr(pm_router, "_load_dataset_features", lambda dataset_id: ({"TP53"}, {"gene"}, []))

    fake_res = SimpleNamespace(
        pathway=SimpleNamespace(pathway_id="hsa00010", name="Glycolysis", source="KEGG"),
        input_features=3,
    )
    monkeypatch.setattr(pm_router, "perform_pathway_enrichment", lambda **kwargs: [fake_res])

    resp = client.get(f"/api/pathway-maps/enriched/{'0'*8}-{'0'*4}-{'0'*4}-{'0'*4}-{'0'*12}")
    assert resp.status_code == 200
    data = resp.json()
    assert data[0]["pathway_id"] == "hsa00010"


def test_expression_endpoint_returns_fold_changes(monkeypatch):
    from amprenta_rag.api.routers import pathway_maps as pm_router

    feats = [
        SimpleNamespace(name="TP53", metadata={"log2fc": 1.25}, external_ids=None),
        SimpleNamespace(name="BRCA1", metadata=None, external_ids={"log2fc": -0.5}),
    ]
    monkeypatch.setattr(pm_router, "_load_dataset_features", lambda dataset_id: (set(), set(), feats))

    resp = client.get(f"/api/pathway-maps/expression/{'0'*8}-{'0'*4}-{'0'*4}-{'0'*4}-{'0'*12}")
    assert resp.status_code == 200
    data = resp.json()
    assert data["gene_expression"]["TP53"] == 1.25
    assert data["gene_expression"]["BRCA1"] == -0.5


def test_invalid_pathway_id_handled(monkeypatch):
    from amprenta_rag.api.routers import pathway_maps as pm_router

    fake_struct = SimpleNamespace(pathway_id="bad", name="bad", organism="", nodes=[], edges=[])
    monkeypatch.setattr(pm_router, "get_pathway_structure", lambda pathway_id: fake_struct)

    resp = client.get("/api/pathway-maps/structure/bad")
    assert resp.status_code == 200
    data = resp.json()
    assert data["nodes"] == []


