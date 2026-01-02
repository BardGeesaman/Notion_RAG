"""Tests for STRING API client."""

from unittest.mock import patch
from amprenta_rag.connectivity.string_client import (
    StringInteraction,
    resolve_gene_symbols,
    get_interactions,
    interactions_to_cytoscape,
)


class TestStringInteraction:
    def test_model_creation(self):
        interaction = StringInteraction(
            protein_a="9606.ENSP00000269305",
            protein_b="9606.ENSP00000253339",
            score=950,
            gene_a="TP53",
            gene_b="MDM2",
        )
        assert interaction.score == 950
        assert interaction.gene_a == "TP53"


class TestResolveGeneSymbols:
    @patch("amprenta_rag.connectivity.string_client._string_request")
    def test_resolve_success(self, mock_request):
        mock_request.return_value = [
            {"queryItem": "TP53", "stringId": "9606.ENSP00000269305"},
            {"queryItem": "BRCA1", "stringId": "9606.ENSP00000312236"},
        ]
        result = resolve_gene_symbols(["TP53", "BRCA1"])
        assert "TP53" in result
        assert result["TP53"] == "9606.ENSP00000269305"

    def test_empty_input(self):
        result = resolve_gene_symbols([])
        assert result == {}


class TestGetInteractions:
    @patch("amprenta_rag.connectivity.string_client._string_request")
    def test_get_interactions_success(self, mock_request):
        mock_request.return_value = [
            {
                "stringId_A": "9606.ENSP00000269305",
                "stringId_B": "9606.ENSP00000253339",
                "preferredName_A": "TP53",
                "preferredName_B": "MDM2",
                "score": 0.95,
            }
        ]
        result = get_interactions(["TP53", "MDM2"])
        assert len(result) == 1
        assert result[0].score == 950
        assert result[0].gene_a == "TP53"

    @patch("amprenta_rag.connectivity.string_client._string_request")
    def test_api_error_returns_empty(self, mock_request):
        mock_request.side_effect = Exception("API Error")
        result = get_interactions(["TP53"])
        assert result == []


class TestInteractionsToCytoscape:
    def test_conversion(self):
        interactions = [
            StringInteraction(
                protein_a="9606.A", protein_b="9606.B",
                score=900, gene_a="GENE1", gene_b="GENE2"
            )
        ]
        nodes, edges = interactions_to_cytoscape(interactions)
        assert len(nodes) == 2
        assert len(edges) == 1
        assert edges[0]["data"]["confidence"] == 0.9
