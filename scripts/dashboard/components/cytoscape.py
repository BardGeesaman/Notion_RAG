from __future__ import annotations

import json
from typing import List, Dict

import streamlit as st
import streamlit.components.v1 as components


def render_cytoscape(nodes: List[Dict], edges: List[Dict], height: int = 600, layout: str = "cose") -> None:
    """
    Render an interactive Cytoscape.js network inside Streamlit.

    Args:
        nodes: List of node dicts; each should have at least an 'id' and optional 'label'.
        edges: List of edge dicts; each should have 'source' and 'target' ids.
        height: Pixel height of the visualization.
        layout: Cytoscape layout name (e.g., "cose", "grid", "circle").
    """
    elements = []
    for n in nodes:
        nid = n.get("id")
        if not nid:
            continue
        elements.append({"data": {"id": nid, "label": n.get("label", nid)}})
    for e in edges:
        src = e.get("source")
        tgt = e.get("target")
        if not src or not tgt:
            continue
        eid = e.get("id") or f"{src}->{tgt}"
        elements.append({"data": {"id": eid, "source": src, "target": tgt}})

    cy_elements = json.dumps(elements)
    layout_json = json.dumps({"name": layout})

    html = f"""
    <div id="cy" style="width: 100%; height: {height}px; border: 1px solid #ddd; border-radius: 6px;"></div>
    <div id="selection" style="margin-top:8px; font-family: system-ui; font-size: 14px; color:#333;">Selected: (none)</div>
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
              'background-color': '#4e79a7',
              'width': '32px',
              'height': '32px',
              'font-size': '12px',
              'text-wrap': 'wrap'
            }}
          }},
          {{
            selector: 'edge',
            style: {{
              'width': 2,
              'line-color': '#9aa0a6',
              'curve-style': 'bezier',
              'target-arrow-color': '#9aa0a6',
              'target-arrow-shape': 'triangle'
            }}
          }},
          {{
            selector: ':selected',
            style: {{
              'background-color': '#e15759',
              'line-color': '#e15759',
              'target-arrow-color': '#e15759'
            }}
          }}
        ],
        layout: layout,
        wheelSensitivity: 0.2
      }});

      const selBox = document.getElementById('selection');
      cy.on('tap', 'node', (evt) => {{
        const n = evt.target;
        selBox.innerText = `Selected node: ${'{'}n.id(){'}'} (label=${'{'}n.data('label'){'}'})`;
      }});
      cy.on('tap', 'edge', (evt) => {{
        const e = evt.target;
        selBox.innerText = `Selected edge: ${'{'}e.id(){'}'} (${ '{'}e.data('source'){'}'} → ${'{'}e.data('target'){'}'})`;
      }});
    </script>
    """
    components.html(html, height=height + 60)

