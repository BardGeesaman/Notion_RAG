from __future__ import annotations

import asyncio
from types import SimpleNamespace
from unittest.mock import MagicMock, patch
from uuid import uuid4

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


class TestAsyncPathwayMapsAPI:
    """Test async execution of pathway maps API endpoints."""

    @pytest.mark.asyncio
    async def test_structure_async(self):
        """Test pathway structure endpoint executes asynchronously."""
        with patch('amprenta_rag.api.routers.pathway_maps._sync_get_pathway_structure') as mock_get_structure:
            # Mock the structure function
            mock_structure = MagicMock()
            mock_structure.pathway_id = "hsa00010"
            mock_structure.name = "Async Test Pathway"
            mock_structure.organism = "hsa"
            # Create proper mock objects with string attributes
            mock_node = MagicMock()
            mock_node.id = "node1"
            mock_node.name = "Test Gene"
            mock_node.type = "gene"
            mock_node.x = 10.0
            mock_node.y = 20.0
            mock_node.kegg_ids = ["hsa:123"]
            mock_structure.nodes = [mock_node]
            
            mock_edge = MagicMock()
            mock_edge.source = "node1"
            mock_edge.target = "node2"
            mock_edge.type = "PPrel"
            mock_edge.subtype = "activation"
            mock_structure.edges = [mock_edge]
            mock_get_structure.return_value = mock_structure
            
            # Mock edge_style function
            with patch('amprenta_rag.api.routers.pathway_maps.edge_style') as mock_edge_style:
                mock_edge_style.return_value = {"style": "solid", "color": "#2E8B57"}
                
                from amprenta_rag.api.routers.pathway_maps import get_structure
                
                result = await get_structure("hsa00010")
                
                # Verify async execution and result
                assert result.pathway_id == "hsa00010"
                assert result.name == "Async Test Pathway"
                assert len(result.nodes) == 1
                assert result.nodes[0].name == "Test Gene"
                assert len(result.edges) == 1
                assert result.edges[0].color == "#2E8B57"
                mock_get_structure.assert_called_once_with("hsa00010")

    @pytest.mark.asyncio
    async def test_search_async(self):
        """Test KEGG search endpoint executes asynchronously."""
        with patch('amprenta_rag.api.routers.pathway_maps._sync_search_kegg_pathways') as mock_search:
            # Mock the search function
            mock_search.return_value = [
                {
                    "pathway_id": "hsa00010",
                    "name": "Async Glycolysis Pathway",
                    "organism": "hsa",
                    "gene_count": 0
                },
                {
                    "pathway_id": "hsa00020",
                    "name": "Async TCA Cycle",
                    "organism": "hsa", 
                    "gene_count": 0
                }
            ]
            
            from amprenta_rag.api.routers.pathway_maps import search_pathways
            
            result = await search_pathways("glycolysis")
            
            # Verify async execution and result
            assert len(result) == 2
            assert result[0].pathway_id == "hsa00010"
            assert result[0].name == "Async Glycolysis Pathway"
            assert result[1].pathway_id == "hsa00020"
            assert result[1].name == "Async TCA Cycle"
            mock_search.assert_called_once_with("glycolysis")

    @pytest.mark.asyncio
    async def test_enrich_async(self):
        """Test pathway enrichment endpoint executes asynchronously."""
        with patch('amprenta_rag.api.routers.pathway_maps._sync_perform_pathway_enrichment') as mock_enrich:
            with patch('amprenta_rag.api.routers.pathway_maps._load_dataset_features') as mock_load:
                # Mock dataset features
                features = {"TP53", "BRCA1", "MYC"}
                types = {"gene"}
                mock_load.return_value = (features, types, [])
                
                # Mock enrichment results
                mock_result = MagicMock()
                mock_result.pathway.pathway_id = "hsa05200"
                mock_result.pathway.name = "Async Cancer Pathway"
                mock_result.input_features = 3
                mock_enrich.return_value = [mock_result]
                
                from amprenta_rag.api.routers.pathway_maps import enriched_pathways
                
                dataset_id = uuid4()
                result = await enriched_pathways(dataset_id)
                
                # Verify async execution and result
                assert len(result) == 1
                assert result[0].pathway_id == "hsa05200"
                assert result[0].name == "Async Cancer Pathway"
                assert result[0].gene_count == 3
                mock_load.assert_called_once_with(dataset_id)
                mock_enrich.assert_called_once_with(features, types)

    @pytest.mark.asyncio
    async def test_concurrent_rate_limited(self):
        """Test concurrent requests execute properly."""
        with patch('amprenta_rag.api.routers.pathway_maps._sync_search_kegg_pathways') as mock_search:
            # Mock function to return different results for different calls
            def mock_search_func(query):
                return [
                    {
                        "pathway_id": f"test_{query.replace(' ', '_')}",
                        "name": f"Pathway for {query}",
                        "organism": "hsa",
                        "gene_count": 0
                    }
                ]
            
            mock_search.side_effect = mock_search_func
            
            from amprenta_rag.api.routers.pathway_maps import search_pathways
            
            async def search_for_query(query: str):
                return await search_pathways(query)
            
            # Make 3 concurrent searches
            tasks = [
                search_for_query("cancer"),
                search_for_query("diabetes"),
                search_for_query("alzheimer"),
            ]
            
            results = await asyncio.gather(*tasks)
            
            # All searches should succeed
            assert len(results) == 3
            assert results[0][0].name == "Pathway for cancer"
            assert results[1][0].name == "Pathway for diabetes"
            assert results[2][0].name == "Pathway for alzheimer"
            
            # Verify all calls were made
            assert mock_search.call_count == 3


