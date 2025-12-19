from __future__ import annotations

import streamlit as st

from scripts.dashboard.components.cytoscape import render_cytoscape
from scripts.dashboard.db_session import db_session
from amprenta_rag.database.models import Program


def _load_hierarchy(program_id: str | None):
    with db_session() as db:
        query = db.query(Program)
        if program_id:
            query = query.filter(Program.id == program_id)
        programs = query.order_by(Program.name).all()
        data = []
        for p in programs:
            data.append(
                {
                    "type": "program",
                    "id": str(p.id),
                    "label": p.name or "Program",
                    "ref": p,
                }
            )
            for exp in p.experiments:
                data.append(
                    {
                        "type": "experiment",
                        "id": str(exp.id),
                        "label": exp.name or "Experiment",
                        "ref": exp,
                        "parent": str(p.id),
                    }
                )
                for ds in exp.datasets:
                    data.append(
                        {
                            "type": "dataset",
                            "id": str(ds.id),
                            "label": ds.name or "Dataset",
                            "ref": ds,
                            "parent": str(exp.id),
                        }
                    )
        # signatures linked to datasets via dataset_signature_assoc relationship
        for p in programs:
            for ds in p.datasets:
                for sig in ds.signatures:
                    data.append(
                        {
                            "type": "signature",
                            "id": str(sig.id),
                            "label": sig.name or "Signature",
                            "ref": sig,
                            "parent": str(ds.id),
                        }
                    )
        return data


def _build_cytoscape_elements(data):
    colors = {
        "program": "#4e79a7",
        "experiment": "#f28e2b",
        "dataset": "#59a14f",
        "signature": "#e15759",
    }
    nodes = []
    edges = []
    seen = set()
    refs = {}
    for item in data:
        nid = item["id"]
        if nid in seen:
            continue
        seen.add(nid)
        refs[nid] = item["ref"]
        nodes.append(
            {
                "id": nid,
                "label": item.get("label", ""),
                "bg": colors.get(item["type"], "#6c757d"),
            }
        )
        parent = item.get("parent")
        if parent:
            edges.append({"source": parent, "target": nid})
    return nodes, edges, refs


def render() -> None:
    st.header("Entity Relationship Explorer")
    st.caption("Program → Experiment → Dataset → Signature relationships")

    with db_session() as db:
        programs = db.query(Program).order_by(Program.name).all()
    prog_opts = {p.name or str(p.id): str(p.id) for p in programs}
    selected_prog = st.selectbox("Filter by Program", ["All"] + list(prog_opts.keys()))
    prog_id = prog_opts.get(selected_prog) if selected_prog != "All" else None

    data = _load_hierarchy(prog_id)
    if not data:
        st.info("No entities found.")
        return

    nodes, edges, refs = _build_cytoscape_elements(data)

    # Use Cytoscape component
    render_cytoscape(nodes, edges, height=700, layout="cose")

    st.sidebar.markdown("### Entity Details")
    st.sidebar.info("Click a node to view details. (Selection display is in the component.)")

