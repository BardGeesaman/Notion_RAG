"""Expression overlay service for network visualization.

Provides functions to fetch expression data and convert to node colors
for Cytoscape.js network visualization.
"""

from __future__ import annotations

import logging
from typing import Dict, List, Optional, Tuple
from uuid import UUID

from amprenta_rag.database.session import db_session
from amprenta_rag.database.models import Feature

logger = logging.getLogger(__name__)

# Colormap definitions (diverging blue-white-red for fold change)
DIVERGING_COLORMAP = {
    "name": "diverging",
    "colors": [
        (-3.0, "#2166ac"),   # Strong downregulation - blue
        (-2.0, "#4393c3"),
        (-1.0, "#92c5de"),
        (-0.5, "#d1e5f0"),
        (0.0, "#f7f7f7"),    # No change - white
        (0.5, "#fddbc7"),
        (1.0, "#f4a582"),
        (2.0, "#d6604d"),
        (3.0, "#b2182b"),    # Strong upregulation - red
    ]
}

SEQUENTIAL_COLORMAP = {
    "name": "sequential",
    "colors": [
        (0.0, "#f7fcf5"),    # Low - light green
        (0.25, "#c7e9c0"),
        (0.5, "#74c476"),
        (0.75, "#238b45"),
        (1.0, "#00441b"),    # High - dark green
    ]
}

# Dataset-feature association table
try:
    from amprenta_rag.database.models import dataset_feature_assoc
except ImportError:
    dataset_feature_assoc = None


def get_expression_for_genes(
    gene_symbols: List[str],
    dataset_id: UUID,
) -> Dict[str, float]:
    """
    Fetch expression values (log2FC) for a list of genes from a dataset.
    
    Args:
        gene_symbols: List of gene symbols to query
        dataset_id: UUID of the dataset
    
    Returns:
        Dict mapping gene symbol to expression value (log2FC or similar)
    """
    if not gene_symbols:
        return {}
    
    result: Dict[str, float] = {}
    
    try:
        with db_session() as db:
            # Query features associated with the dataset
            if dataset_feature_assoc is not None:
                query = (
                    db.query(Feature)
                    .join(dataset_feature_assoc, Feature.id == dataset_feature_assoc.c.feature_id)
                    .filter(dataset_feature_assoc.c.dataset_id == dataset_id)
                    .filter(Feature.name.in_(gene_symbols))
                )
            else:
                # Fallback: direct query on features with dataset_id
                query = db.query(Feature).filter(
                    Feature.name.in_(gene_symbols),
                    Feature.dataset_id == dataset_id,
                )
            
            features = query.all()
            
            for f in features:
                name = getattr(f, "name", None)
                if not name:
                    continue
                
                # Extract log2FC from feature metadata or value field
                value = _extract_expression_value(f)
                if value is not None:
                    result[str(name)] = float(value)
        
        logger.info(f"Fetched expression for {len(result)}/{len(gene_symbols)} genes from dataset {dataset_id}")
        
    except Exception as e:
        logger.error(f"Failed to fetch expression data: {e}")
    
    return result


def _extract_expression_value(feature: Feature) -> Optional[float]:
    """Extract expression value from a feature object."""
    # Try common attribute names for expression values
    for attr in ["log2fc", "log2_fold_change", "fold_change", "value", "expression"]:
        val = getattr(feature, attr, None)
        if val is not None:
            try:
                return float(val)
            except (ValueError, TypeError):
                continue
    
    # Try metadata dict
    metadata = getattr(feature, "metadata", None) or {}
    if isinstance(metadata, dict):
        for key in ["log2fc", "log2_fold_change", "fold_change", "value"]:
            if key in metadata:
                try:
                    return float(metadata[key])
                except (ValueError, TypeError):
                    continue
    
    return None


def compute_node_colors(
    expression: Dict[str, float],
    colormap: str = "diverging",
    vmin: Optional[float] = None,
    vmax: Optional[float] = None,
) -> Dict[str, str]:
    """
    Map expression values to hex colors.
    
    Args:
        expression: Dict of gene_symbol → expression_value
        colormap: "diverging" (blue-white-red) or "sequential" (green)
        vmin: Minimum value for color scaling (default: min of data or -3)
        vmax: Maximum value for color scaling (default: max of data or 3)
    
    Returns:
        Dict mapping gene_symbol to hex color string
    """
    if not expression:
        return {}
    
    # Select colormap
    cmap = DIVERGING_COLORMAP if colormap == "diverging" else SEQUENTIAL_COLORMAP
    color_stops = cmap["colors"]
    
    # Determine value range
    values = list(expression.values())
    if colormap == "diverging":
        # Symmetric around 0 for fold change
        max_abs = max(abs(v) for v in values) if values else 3.0
        actual_vmin = vmin if vmin is not None else -min(max_abs, 3.0)
        actual_vmax = vmax if vmax is not None else min(max_abs, 3.0)
    else:
        # Sequential: 0 to max
        actual_vmin = vmin if vmin is not None else min(values) if values else 0.0
        actual_vmax = vmax if vmax is not None else max(values) if values else 1.0
    
    # Map values to colors
    result: Dict[str, str] = {}
    for gene, value in expression.items():
        result[gene] = _value_to_color(value, color_stops, actual_vmin, actual_vmax)
    
    return result


def _value_to_color(
    value: float,
    color_stops: List[Tuple[float, str]],
    vmin: float,
    vmax: float,
) -> str:
    """Interpolate a value to a hex color using color stops."""
    # Normalize value to colormap range
    cmap_min = color_stops[0][0]
    cmap_max = color_stops[-1][0]
    
    # Scale value from [vmin, vmax] to [cmap_min, cmap_max]
    if vmax == vmin:
        normalized = (cmap_min + cmap_max) / 2
    else:
        t = (value - vmin) / (vmax - vmin)
        t = max(0.0, min(1.0, t))  # Clamp to [0, 1]
        normalized = cmap_min + t * (cmap_max - cmap_min)
    
    # Find surrounding color stops
    for i in range(len(color_stops) - 1):
        stop_val, stop_color = color_stops[i]
        next_val, next_color = color_stops[i + 1]
        
        if stop_val <= normalized <= next_val:
            # Interpolate between stops
            if next_val == stop_val:
                t = 0.5
            else:
                t = (normalized - stop_val) / (next_val - stop_val)
            return _interpolate_hex_colors(stop_color, next_color, t)
    
    # Return endpoint color if out of range
    if normalized <= color_stops[0][0]:
        return color_stops[0][1]
    return color_stops[-1][1]


def _interpolate_hex_colors(color1: str, color2: str, t: float) -> str:
    """Linearly interpolate between two hex colors."""
    # Parse hex colors
    r1, g1, b1 = int(color1[1:3], 16), int(color1[3:5], 16), int(color1[5:7], 16)
    r2, g2, b2 = int(color2[1:3], 16), int(color2[3:5], 16), int(color2[5:7], 16)
    
    # Interpolate
    r = int(r1 + t * (r2 - r1))
    g = int(g1 + t * (g2 - g1))
    b = int(b1 + t * (b2 - b1))
    
    return f"#{r:02x}{g:02x}{b:02x}"


def get_expression_overlay(
    gene_symbols: List[str],
    dataset_id: UUID,
    colormap: str = "diverging",
) -> Tuple[Dict[str, str], Dict[str, float]]:
    """
    Convenience function to get expression data and colors in one call.
    
    Args:
        gene_symbols: Gene symbols to query
        dataset_id: Dataset UUID
        colormap: "diverging" or "sequential"
    
    Returns:
        Tuple of (node_colors dict, raw_expression dict)
    """
    expression = get_expression_for_genes(gene_symbols, dataset_id)
    colors = compute_node_colors(expression, colormap)
    return colors, expression


__all__ = [
    "get_expression_for_genes",
    "compute_node_colors",
    "get_expression_overlay",
    "DIVERGING_COLORMAP",
    "SEQUENTIAL_COLORMAP",
]
