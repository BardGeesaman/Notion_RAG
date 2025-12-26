from __future__ import annotations

import pytest


def test_compute_overlay_matches_genes():
    from amprenta_rag.analysis.pathway.overlay import compute_overlay
    from amprenta_rag.analysis.pathway.structure import PathwayNode, PathwayStructure

    st = PathwayStructure(
        pathway_id="hsa00010",
        name="x",
        organism="hsa",
        nodes=[
            PathwayNode(id="1", name="TP53", type="gene", x=0.1, y=0.1, kegg_ids=["hsa:7157"]),
            PathwayNode(id="2", name="Glucose", type="compound", x=0.2, y=0.2, kegg_ids=["cpd:C00031"]),
        ],
        edges=[],
    )
    out = compute_overlay(st, {"TP53": 1.5})
    assert len(out) == 1
    assert out[0].node_id == "1"
    assert out[0].gene_symbol.lower() == "tp53"


def test_overlay_missing_genes_skipped():
    from amprenta_rag.analysis.pathway.overlay import compute_overlay
    from amprenta_rag.analysis.pathway.structure import PathwayNode, PathwayStructure

    st = PathwayStructure(
        pathway_id="hsa00010",
        name="x",
        organism="hsa",
        nodes=[PathwayNode(id="1", name="TP53", type="gene", x=0.1, y=0.1, kegg_ids=["hsa:7157"])],
        edges=[],
    )
    out = compute_overlay(st, {"NOT_A_GENE": 2.0})
    assert out == []


def test_value_to_color_diverging():
    pytest.importorskip("matplotlib")

    from amprenta_rag.analysis.pathway.overlay import value_to_color

    c1 = value_to_color(-2.0, vmin=-2.0, vmax=2.0, colormap="RdBu_r")
    c2 = value_to_color(2.0, vmin=-2.0, vmax=2.0, colormap="RdBu_r")
    assert isinstance(c1, str) and c1.startswith("#") and len(c1) == 7
    assert isinstance(c2, str) and c2.startswith("#") and len(c2) == 7
    assert c1 != c2


def test_colormap_options():
    pytest.importorskip("matplotlib")

    from amprenta_rag.analysis.pathway.overlay import value_to_color

    for cmap in ["RdBu_r", "coolwarm", "viridis"]:
        c = value_to_color(0.5, vmin=-2.0, vmax=2.0, colormap=cmap)
        assert isinstance(c, str) and c.startswith("#") and len(c) == 7


