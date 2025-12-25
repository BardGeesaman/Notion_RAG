"""Enhanced Cytoscape.js component for compound-target networks."""

from __future__ import annotations

import json
from typing import Any, Dict, List, Optional

import streamlit as st
import streamlit.components.v1 as components


DEFAULT_LAYOUTS = ["cose", "concentric", "circle"]


def render_compound_target_network(
    nodes: List[Dict[str, Any]],
    edges: List[Dict[str, Any]],
    layout: str = "cose",
    options: Optional[Dict[str, Any]] = None,
    *,
    height: int = 650,
    key: str = "compound_target_network",
) -> Any:
    """
    Render a Cytoscape network with drug-discovery styling and selection output.

    Notes:
    - Max ~500 nodes recommended; caller should pre-trim.
    - Selection is emitted via Streamlit component value (postMessage) and also displayed.
    """
    opts = options or {}
    lay = layout if layout in DEFAULT_LAYOUTS else "cose"

    # Flatten to Cytoscape elements list
    elements: List[Dict[str, Any]] = []
    for n in nodes or []:
        if not isinstance(n, dict):
            continue
        data = n.get("data") if isinstance(n.get("data"), dict) else n
        nid = data.get("id")
        if not nid:
            continue
        elements.append({"data": data})

    for e in edges or []:
        if not isinstance(e, dict):
            continue
        data = e.get("data") if isinstance(e.get("data"), dict) else e
        if not data.get("source") or not data.get("target"):
            continue
        if not data.get("id"):
            data["id"] = f"{data['source']}->{data['target']}"
        elements.append({"data": data})

    if len([el for el in elements if "source" not in el["data"]]) > 500:
        st.warning("Network is large; showing first 500 nodes to prevent browser overload.")
        node_elems = [el for el in elements if "source" not in el["data"]][:500]
        edge_elems = [el for el in elements if "source" in el["data"]]
        keep = {el["data"]["id"] for el in node_elems}
        edge_elems = [el for el in edge_elems if el["data"].get("source") in keep and el["data"].get("target") in keep]
        elements = node_elems + edge_elems

    cy_elements = json.dumps(elements)
    layout_json = json.dumps({"name": lay})
    opt_json = json.dumps(opts)

    html = f"""
<div id="cy" style="width: 100%; height: {int(height)}px; border: 1px solid #ddd; border-radius: 6px;"></div>
<div id="hover" style="margin-top:8px; font-family: system-ui; font-size: 13px; color:#333;">Hover: (none)</div>
<div id="selection" style="margin-top:6px; font-family: system-ui; font-size: 13px; color:#333;">Selected: (none)</div>
<link rel="stylesheet" href="https://unpkg.com/cytoscape@3.23.0/dist/cytoscape.min.css" />
<script src="https://unpkg.com/cytoscape@3.23.0/dist/cytoscape.min.js"></script>
<script>
  const elements = {cy_elements};
  const layout = {layout_json};
  const options = {opt_json};
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
          'shape': 'data(shape)',
          'width': 34,
          'height': 34,
          'font-size': '11px',
          'text-wrap': 'wrap',
          'text-max-width': 90,
          'border-width': 1,
          'border-color': '#ffffff'
        }}
      }},
      {{
        selector: 'node[node_type = "target"]',
        style: {{
          'background-color': '#59a14f'
        }}
      }},
      {{
        selector: 'edge',
        style: {{
          'width': 'data(pic50_width)',
          'line-color': 'data(edge_color)',
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
    wheelSensitivity: 0.2,
    ...options
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
    const label = d.label || d.id;
    const extra = d.compound_id ? `compound_id=${{d.compound_id}}` : (d.target_name ? `target=${{d.target_name}}` : '');
    hoverBox.innerText = `Hover node: ${{label}} ${{extra ? '(' + extra + ')' : ''}}`;
  }});
  cy.on('mouseover', 'edge', (evt) => {{
    const e = evt.target;
    const d = e.data();
    const best = d.best_ic50_nm ? `${{d.best_ic50_nm}} nM` : 'n/a';
    const cnt = d.assays_count || d.provenance?.assays_count || 'n/a';
    hoverBox.innerText = `Hover edge: IC50(best)=${{best}} · assays=${{cnt}}`;
  }});
  cy.on('mouseout', 'node, edge', () => {{
    hoverBox.innerText = 'Hover: (none)';
  }});

  cy.on('tap', 'node', (evt) => {{
    const n = evt.target;
    const d = n.data();
    selBox.innerText = `Selected node: ${{d.id}} (label=${{d.label || ''}})`;
    setValue({{kind: 'node', data: d}});
  }});
  cy.on('tap', 'edge', (evt) => {{
    const e = evt.target;
    const d = e.data();
    selBox.innerText = `Selected edge: ${{d.id}} (${{d.source}} → ${{d.target}})`;
    setValue({{kind: 'edge', data: d}});
  }});
</script>
"""

    return components.html(html, height=int(height) + 80, key=key)


__all__ = ["render_compound_target_network", "DEFAULT_LAYOUTS"]


