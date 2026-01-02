"""Tests for expression overlay service."""

from amprenta_rag.services.expression_overlay import (
    compute_node_colors,
    _interpolate_hex_colors,
)


class TestComputeNodeColors:
    def test_diverging_colormap(self):
        expression = {"GENE1": -2.0, "GENE2": 0.0, "GENE3": 2.0}
        colors = compute_node_colors(expression, "diverging")
        assert len(colors) == 3
        # Downregulated should be blue-ish
        assert colors["GENE1"].startswith("#")
        # No change should be white-ish
        assert colors["GENE2"].startswith("#")

    def test_sequential_colormap(self):
        expression = {"GENE1": 0.0, "GENE2": 0.5, "GENE3": 1.0}
        colors = compute_node_colors(expression, "sequential")
        assert len(colors) == 3

    def test_empty_expression(self):
        result = compute_node_colors({})
        assert result == {}


class TestInterpolateHexColors:
    def test_midpoint(self):
        result = _interpolate_hex_colors("#000000", "#ffffff", 0.5)
        # Should be gray-ish (around #7f7f7f)
        assert result.startswith("#")
        assert len(result) == 7

    def test_endpoints(self):
        assert _interpolate_hex_colors("#ff0000", "#0000ff", 0.0) == "#ff0000"
        assert _interpolate_hex_colors("#ff0000", "#0000ff", 1.0) == "#0000ff"
