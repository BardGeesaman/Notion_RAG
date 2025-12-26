"""Overlay computation for pathway structures (e.g., fold-change coloring)."""

from __future__ import annotations

from dataclasses import dataclass
from typing import Dict, List

from amprenta_rag.analysis.pathway.structure import PathwayStructure


@dataclass(frozen=True)
class OverlayData:
    node_id: str
    gene_symbol: str
    value: float
    color: str
    label: str


def value_to_color(value: float, vmin: float, vmax: float, colormap: str) -> str:
    """
    Convert a numeric value into a hex color using a diverging colormap.
    """
    try:
        import matplotlib.cm as cm  # type: ignore
        import matplotlib.colors as mcolors  # type: ignore

        norm = mcolors.TwoSlopeNorm(vmin=float(vmin), vcenter=0.0, vmax=float(vmax))
        cmap = cm.get_cmap(str(colormap))
        rgba = cmap(norm(float(value)))
        return str(mcolors.to_hex(rgba, keep_alpha=False))
    except Exception:
        return "#9E9E9E"


def _candidate_keys(node_name: str, kegg_ids: List[str]) -> List[str]:
    keys: List[str] = []
    if node_name:
        keys.append(node_name)
        keys.append(node_name.upper())
        keys.append(node_name.lower())
    for kid in kegg_ids or []:
        if not kid:
            continue
        keys.append(kid)
        keys.append(kid.upper())
        keys.append(kid.lower())
        if ":" in kid:
            suffix = kid.split(":", 1)[1]
            keys.append(suffix)
            keys.append(suffix.upper())
            keys.append(suffix.lower())
    # Dedup preserving order
    out: List[str] = []
    seen = set()
    for k in keys:
        if k and k not in seen:
            out.append(k)
            seen.add(k)
    return out


def compute_overlay(
    structure: PathwayStructure,
    expression_data: Dict[str, float],
    colormap: str = "RdBu_r",
    vmin: float = -2.0,
    vmax: float = 2.0,
) -> List[OverlayData]:
    """
    Compute overlay data for pathway nodes.

    Best-effort matching:
    - expression_data keys are matched against node.kegg_ids and node.name (case-insensitive).
    """
    if not isinstance(expression_data, dict):
        return []
    # Case-insensitive lookup for gene symbols
    expr_lut: Dict[str, float] = {}
    for k, v in expression_data.items():
        try:
            expr_lut[str(k).lower()] = float(v)
        except Exception:
            continue

    out: List[OverlayData] = []
    for n in structure.nodes or []:
        if (n.type or "").lower() != "gene":
            continue
        matched_key = None
        matched_val = None
        for cand in _candidate_keys(n.name, n.kegg_ids):
            val = expr_lut.get(str(cand).lower())
            if val is not None:
                matched_key = cand
                matched_val = float(val)
                break
        if matched_key is None or matched_val is None:
            continue

        color = value_to_color(matched_val, vmin=vmin, vmax=vmax, colormap=colormap)
        label = f"{matched_key}: {matched_val:+.2f}"
        out.append(OverlayData(node_id=n.id, gene_symbol=str(matched_key), value=float(matched_val), color=color, label=label))
    return out


__all__ = ["OverlayData", "compute_overlay", "value_to_color"]


