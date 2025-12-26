"""Cytoscape.js component for rendering KEGG pathway maps with overlays."""

from __future__ import annotations

import base64
import io
import json
from typing import Any, Dict, List, Optional

import streamlit as st
import streamlit.components.v1 as components


def render_pathway_map(
    nodes: List[Dict[str, Any]],
    edges: List[Dict[str, Any]],
    overlays: Optional[List[Dict[str, Any]]] = None,
    layout: str = "preset",
    *,
    height: int = 700,
    key: str = "pathway_map",
) -> Any:
    """
    Render a KEGG pathway map in Cytoscape.

    Args:
        nodes/edges: API-shaped dicts (PathwayStructureResponse.nodes/edges).
        overlays: Optional overlay dicts with fields {node_id,color,label,...}.
        layout: "preset" (use KGML coords) or "cose".

    Returns:
        Streamlit component value on interactions (selection or PNG export).
    """
    overlay_map: Dict[str, Dict[str, Any]] = {}
    for o in overlays or []:
        if isinstance(o, dict) and o.get("node_id"):
            overlay_map[str(o["node_id"])] = o

    elements: List[Dict[str, Any]] = []
    scale_x = 1000.0
    scale_y = float(height)

    # Nodes
    for n in nodes or []:
        if not isinstance(n, dict):
            continue
        nid = n.get("id")
        if not nid:
            continue
        ntype = (n.get("type") or "").lower()
        shape = "ellipse"
        if ntype == "gene":
            shape = "round-rectangle"
        elif ntype == "compound":
            shape = "ellipse"
        elif ntype in ("sub-pathway", "subpathway", "map"):
            shape = "hexagon"

        overlay = overlay_map.get(str(nid))
        color = overlay.get("color") if isinstance(overlay, dict) else None
        if not color:
            color = "#9E9E9E"

        data = {
            "id": str(nid),
            "label": str(n.get("name") or nid),
            "node_type": ntype or "unknown",
            "shape": shape,
            "color": str(color),
            "kegg_ids": n.get("kegg_ids") or [],
            "overlay_label": overlay.get("label") if isinstance(overlay, dict) else None,
            "overlay_value": overlay.get("value") if isinstance(overlay, dict) else None,
        }
        el: Dict[str, Any] = {"data": data}
        if layout == "preset":
            try:
                x = float(n.get("x") or 0.0) * scale_x
                y = float(n.get("y") or 0.0) * scale_y
                el["position"] = {"x": x, "y": y}
                el["locked"] = True
            except Exception:
                pass
        elements.append(el)

    # Edges
    for e in edges or []:
        if not isinstance(e, dict):
            continue
        src = e.get("source")
        tgt = e.get("target")
        if not src or not tgt:
            continue
        style = (e.get("style") or "solid").lower()
        line_style = "solid" if style != "dashed" else "dashed"
        color = e.get("color") or "#9E9E9E"
        data = {
            "id": e.get("id") or f"{src}->{tgt}",
            "source": str(src),
            "target": str(tgt),
            "edge_type": e.get("type") or "relation",
            "subtype": e.get("subtype"),
            "line_style": line_style,
            "edge_color": str(color),
        }
        elements.append({"data": data})

    cy_elements = json.dumps(elements)
    layout_name = "preset" if layout == "preset" else "cose"
    layout_json = json.dumps({"name": layout_name, "animate": False})

    html = f"""
<div style="display:flex; gap: 8px; align-items:center; margin-bottom:6px;">
  <button id="exportPng" style="padding:6px 10px; border: 1px solid #ccc; border-radius: 6px; background: #fff; cursor:pointer;">Export PNG</button>
  <span style="font-family: system-ui; font-size: 12px; color:#555;">Tip: Use mouse wheel to zoom, drag to pan.</span>
</div>
<div id="cy" style="width: 100%; height: {int(height)}px; border: 1px solid #ddd; border-radius: 6px;"></div>
<div id="hover" style="margin-top:8px; font-family: system-ui; font-size: 13px; color:#333;">Hover: (none)</div>
<div id="selection" style="margin-top:6px; font-family: system-ui; font-size: 13px; color:#333;">Selected: (none)</div>
<link rel="stylesheet" href="https://unpkg.com/cytoscape@3.23.0/dist/cytoscape.min.css" />
<script src="https://unpkg.com/cytoscape@3.23.0/dist/cytoscape.min.js"></script>
<script>
  const elements = {cy_elements};
  const layout = {layout_json};
  const cy = cytoscape({{
    container: document.getElementById('cy'),
    elements: elements,
    style: [
      {{
        selector: 'node',
        style: {{
          'label': 'data(label)',
          'text-valign': 'center',
          'color': '#1b1f23',
          'background-color': 'data(color)',
          'shape': 'data(shape)',
          'width': 44,
          'height': 30,
          'font-size': '10px',
          'text-wrap': 'wrap',
          'text-max-width': 110,
          'border-width': 1,
          'border-color': '#ffffff'
        }}
      }},
      {{
        selector: 'edge',
        style: {{
          'width': 2,
          'line-color': 'data(edge_color)',
          'line-style': 'data(line_style)',
          'curve-style': 'bezier',
          'target-arrow-color': 'data(edge_color)',
          'target-arrow-shape': 'triangle'
        }}
      }},
      {{
        selector: ':selected',
        style: {{
          'border-width': 3,
          'border-color': '#111111',
          'line-color': '#111111',
          'target-arrow-color': '#111111'
        }}
      }}
    ],
    layout: layout,
    wheelSensitivity: 0.2
  }});

  const hoverBox = document.getElementById('hover');
  const selBox = document.getElementById('selection');

  function setValue(payload) {{
    try {{
      window.parent.postMessage({{
        type: 'streamlit:setComponentValue',
        value: payload
      }}, '*');
    }} catch (e) {{}}
  }}

  cy.on('mouseover', 'node', (evt) => {{
    const n = evt.target;
    const d = n.data();
    const extra = d.overlay_label ? ` · ${{d.overlay_label}}` : '';
    hoverBox.innerText = `Hover node: ${{d.label || d.id}}${{extra}}`;
  }});
  cy.on('mouseover', 'edge', (evt) => {{
    const e = evt.target;
    const d = e.data();
    const sub = d.subtype ? ` (${{d.subtype}})` : '';
    hoverBox.innerText = `Hover edge: ${{d.edge_type}}${{sub}}`;
  }});
  cy.on('mouseout', 'node, edge', () => {{
    hoverBox.innerText = 'Hover: (none)';
  }});

  cy.on('tap', 'node', (evt) => {{
    const n = evt.target;
    const d = n.data();
    selBox.innerText = `Selected node: ${{d.id}}`;
    setValue({{kind: 'node', data: d}});
  }});
  cy.on('tap', 'edge', (evt) => {{
    const e = evt.target;
    const d = e.data();
    selBox.innerText = `Selected edge: ${{d.id}} (${{d.source}} → ${{d.target}})`;
    setValue({{kind: 'edge', data: d}});
  }});

  document.getElementById('exportPng').addEventListener('click', () => {{
    try {{
      const dataUrl = cy.png({{ full: true, scale: 2, bg: '#ffffff' }});
      setValue({{kind: 'png', data_url: dataUrl}});
    }} catch (e) {{
      setValue({{kind: 'png', error: 'PNG export failed'}});
    }}
  }});
</script>
"""

    return components.html(html, height=int(height) + 120, key=key)


def render_color_legend(vmin: float, vmax: float, colormap: str) -> None:
    """Render a horizontal diverging color legend."""
    try:
        import matplotlib.cm as cm  # type: ignore
        import matplotlib.colors as mcolors  # type: ignore
        import numpy as np  # type: ignore

        cmap = cm.get_cmap(str(colormap))
        norm = mcolors.TwoSlopeNorm(vmin=float(vmin), vcenter=0.0, vmax=float(vmax))
        grad = np.linspace(float(vmin), float(vmax), 256)
        rgba = cmap(norm(grad))
        img = (rgba[:, :3] * 255).astype("uint8")[None, :, :]

        import matplotlib.pyplot as plt  # type: ignore

        fig, ax = plt.subplots(figsize=(6, 0.4), dpi=150)
        ax.imshow(img, aspect="auto")
        ax.set_axis_off()
        buf = io.BytesIO()
        fig.savefig(buf, format="png", bbox_inches="tight", pad_inches=0)
        plt.close(fig)
        b64 = base64.b64encode(buf.getvalue()).decode("utf-8")
        st.image(f"data:image/png;base64,{b64}", caption=f"{colormap}  (vmin={vmin}, vmax={vmax})")
    except Exception:
        st.caption(f"Legend: {colormap} (vmin={vmin}, vmax={vmax})")


__all__ = ["render_pathway_map", "render_color_legend"]


