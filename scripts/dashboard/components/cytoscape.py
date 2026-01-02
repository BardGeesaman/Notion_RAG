"""Enhanced Cytoscape.js component with expression overlay and cluster highlighting."""

from __future__ import annotations

import json
from typing import Dict, List, Optional

import streamlit.components.v1 as components


# 10-color palette for communities
COMMUNITY_COLORS = [
    "#2E86AB", "#A23B72", "#F18F01", "#C73E1D", "#3B1F2B",
    "#4ECDC4", "#C7F9CC", "#9B5DE5", "#F15BB5", "#00BBF9",
]


def render_cytoscape(
    nodes: List[Dict],
    edges: List[Dict],
    height: int = 600,
    layout: str = "cose",
    node_colors: Optional[Dict[str, str]] = None,
    communities: Optional[Dict[str, int]] = None,
    show_edge_labels: bool = False,
    show_legend: bool = True,
    color_by: str = "default",
) -> None:
    """
    Render an interactive Cytoscape.js network with optional overlays.
    
    Args:
        nodes: List of node dicts with 'id', 'label', optional 'data'
        edges: List of edge dicts with 'source', 'target', optional 'label'
        height: Pixel height
        layout: Cytoscape layout (cose, grid, circle, breadthfirst, concentric)
        node_colors: Dict mapping node_id to hex color (for expression overlay)
        communities: Dict mapping node_id to community_id integer
        show_edge_labels: Whether to display edge labels
        show_legend: Whether to show color legend
        color_by: "default", "expression", or "community"
    """
    # Build elements with color data
    elements = []
    node_colors = node_colors or {}
    communities = communities or {}
    
    for n in nodes:
        nid = n.get("id") or (n.get("data") or {}).get("id")
        if not nid:
            continue
        
        label = n.get("label") or (n.get("data") or {}).get("label") or nid
        
        # Determine node color
        if color_by == "expression" and nid in node_colors:
            color = node_colors[nid]
        elif color_by == "community" and nid in communities:
            cid = communities[nid] % len(COMMUNITY_COLORS)
            color = COMMUNITY_COLORS[cid]
        else:
            color = "#4e79a7"  # Default blue
        
        elements.append({
            "data": {
                "id": nid,
                "label": label,
                "color": color,
                "community": communities.get(nid, 0),
            }
        })
    
    for e in edges:
        src = e.get("source") or (e.get("data") or {}).get("source")
        tgt = e.get("target") or (e.get("data") or {}).get("target")
        if not src or not tgt:
            continue
        
        eid = e.get("id") or f"{src}->{tgt}"
        label = e.get("label") or e.get("data", {}).get("label") or ""
        confidence = e.get("confidence") or e.get("data", {}).get("confidence") or e.get("score", 0)
        
        elements.append({
            "data": {
                "id": eid,
                "source": src,
                "target": tgt,
                "label": label,
                "confidence": confidence,
            }
        })
    
    # Generate legend HTML
    legend_html = ""
    if show_legend and color_by == "community" and communities:
        unique_comms = sorted(set(communities.values()))
        legend_items = "".join([
            f'<span style="display:inline-block;width:12px;height:12px;background:{COMMUNITY_COLORS[c % len(COMMUNITY_COLORS)]};margin-right:4px;"></span>Cluster {c}&nbsp;&nbsp;'
            for c in unique_comms[:10]
        ])
        legend_html = f'<div style="font-size:12px;margin-bottom:8px;">{legend_items}</div>'
    elif show_legend and color_by == "expression" and node_colors:
        legend_html = '''
        <div style="font-size:12px;margin-bottom:8px;">
            <span style="background:linear-gradient(to right,#2166ac,#f7f7f7,#b2182b);display:inline-block;width:100px;height:12px;"></span>
            <span style="margin-left:8px;">Low → High Expression</span>
        </div>
        '''
    
    # Build HTML
    cy_elements = json.dumps(elements)
    layout_json = json.dumps({"name": layout})
    edge_label_style = "'label': 'data(label)'," if show_edge_labels else ""
    
    html = f"""
    {legend_html}
    <div style="margin-bottom:8px;">
        <button id="zoom-in" style="padding:4px 10px;margin-right:4px;cursor:pointer;">+</button>
        <button id="zoom-out" style="padding:4px 10px;margin-right:4px;cursor:pointer;">−</button>
        <button id="fit" style="padding:4px 10px;cursor:pointer;">Fit</button>
    </div>
    <div id="cy" style="width:100%;height:{height}px;border:1px solid #ddd;border-radius:6px;"></div>
    <div id="selection" style="margin-top:8px;font-family:system-ui;font-size:14px;color:#333;">Selected: (none)</div>
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
              'width': '32px',
              'height': '32px',
              'font-size': '11px',
              'text-wrap': 'wrap',
              'text-max-width': '80px',
              'border-width': 1,
              'border-color': '#ffffff'
            }}
          }},
          {{
            selector: 'edge',
            style: {{
              'width': 2,
              'line-color': '#9aa0a6',
              'curve-style': 'bezier',
              'target-arrow-color': '#9aa0a6',
              'target-arrow-shape': 'triangle',
              {edge_label_style}
              'font-size': '10px',
              'text-background-color': '#fff',
              'text-background-opacity': 0.7
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

      // Zoom controls
      document.getElementById('zoom-in').onclick = () => cy.zoom(cy.zoom() * 1.2);
      document.getElementById('zoom-out').onclick = () => cy.zoom(cy.zoom() / 1.2);
      document.getElementById('fit').onclick = () => cy.fit();

      // Selection display
      const selBox = document.getElementById('selection');
      cy.on('tap', 'node', (evt) => {{
        const n = evt.target;
        const d = n.data();
        selBox.innerText = `Selected node: ${{d.id}} (label=${{d.label}}, community=${{d.community}})`;
      }});
      cy.on('tap', 'edge', (evt) => {{
        const e = evt.target;
        const d = e.data();
        selBox.innerText = `Selected edge: ${{d.source}} → ${{d.target}} (confidence=${{d.confidence}})`;
      }});
    </script>
    """
    components.html(html, height=height + 100)


__all__ = ["render_cytoscape", "COMMUNITY_COLORS"]