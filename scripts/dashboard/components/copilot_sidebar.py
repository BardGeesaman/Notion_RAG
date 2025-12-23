"""Copilot sidebar component for notebook/code generation.

This component provides quick actions + a prompt box to generate notebook-ready
Python code using the LLM copilot utilities in `amprenta_rag.notebook`.

Intended usage:
    from components.copilot_sidebar import render_copilot_sidebar
"""

from __future__ import annotations

from typing import Dict, List, Optional, Tuple

import streamlit as st

from amprenta_rag.notebook.context import AnalysisContext
from amprenta_rag.notebook.copilot import explain_cell, synthesize_cell
from scripts.dashboard.db_session import db_session


def _copy_code(code: str) -> None:
    """Copy code to clipboard if possible; otherwise render a JS clipboard helper."""
    try:
        import pyperclip  # type: ignore

        pyperclip.copy(code)
        st.success("Copied to clipboard.")
        return
    except Exception:
        # Fallback: browser clipboard via HTML (may be blocked depending on browser perms).
        import streamlit.components.v1 as components

        escaped = (
            code.replace("\\", "\\\\")
            .replace("`", "\\`")
            .replace("${", "\\${")
            .replace("\n", "\\n")
        )
        components.html(
            f"""
            <div style="display:flex; gap:8px; align-items:center;">
              <button id="copilot-copy-btn" style="padding:6px 10px;">Copy</button>
              <span id="copilot-copy-status" style="font-size:12px; color:#666;"></span>
            </div>
            <script>
              const txt = `{escaped}`;
              const btn = document.getElementById("copilot-copy-btn");
              const status = document.getElementById("copilot-copy-status");
              btn.onclick = async () => {{
                try {{
                  await navigator.clipboard.writeText(txt.replace(/\\n/g, "\\n"));
                  status.textContent = "Copied";
                }} catch (e) {{
                  status.textContent = "Copy failed (clipboard permissions)";
                }}
              }};
            </script>
            """,
            height=42,
        )


def _load_entity_options() -> Dict[str, List[Tuple[str, str]]]:
    """Load entity options for dropdown selectors (name label, id)."""
    # Use dashboard import pattern: query via db_session with models imported from db_models module.
    import amprenta_rag.database.models as db_models

    Dataset = db_models.Dataset
    Experiment = db_models.Experiment
    HTSCampaign = db_models.HTSCampaign
    Compound = db_models.Compound

    opts: Dict[str, List[Tuple[str, str]]] = {
        "dataset": [],
        "experiment": [],
        "campaign": [],
        "compound": [],
    }

    with db_session() as db:
        for ds in db.query(Dataset).order_by(Dataset.created_at.desc()).limit(200).all():
            label = f"{ds.name} ({str(ds.id)[:8]})"
            opts["dataset"].append((label, str(ds.id)))

        for exp in db.query(Experiment).order_by(Experiment.created_at.desc()).limit(200).all():
            label = f"{exp.name} ({str(exp.id)[:8]})"
            opts["experiment"].append((label, str(exp.id)))

        for camp in db.query(HTSCampaign).order_by(HTSCampaign.created_at.desc()).limit(200).all():
            label = f"{camp.name} ({str(camp.id)[:8]})"
            opts["campaign"].append((label, str(camp.id)))

        for cpd in db.query(Compound).order_by(Compound.created_at.desc()).limit(200).all():
            name = cpd.name or cpd.corporate_id or "Compound"
            label = f"{name} ({str(cpd.id)[:8]})"
            opts["compound"].append((label, str(cpd.id)))

    return opts


def _default_intent(action: str, ctx: AnalysisContext) -> str:
    action_norm = action.strip().lower().replace("-", "_")
    if action_norm == "load_dataset":
        return "Load dataset and display summary + features table."
    if action_norm == "hts_qc":
        return f"Run HTS QC for campaign_id={ctx.campaign_id or ctx.entity_id} and summarize key metrics."
    if action_norm == "dose_response":
        return f"Fit dose-response for compound_id={ctx.compound_id or ctx.entity_id} and summarize results."
    if action_norm == "publish":
        return "Publish the current analysis results to RAG with tags and a title."
    return f"Generate a notebook cell for action '{action_norm}'."


def render_copilot_sidebar() -> Optional[str]:
    """Render the Copilot sidebar; returns generated code (or None)."""
    st.subheader("Copilot")
    st.caption("Generate notebook-ready code from context + intent.")

    if "copilot_action" not in st.session_state:
        st.session_state["copilot_action"] = "load_dataset"
    if "copilot_generated_code" not in st.session_state:
        st.session_state["copilot_generated_code"] = ""
    if "copilot_code_explanation" not in st.session_state:
        st.session_state["copilot_code_explanation"] = ""

    # Quick actions
    c1, c2 = st.columns(2)
    with c1:
        if st.button("Load Dataset", use_container_width=True):
            st.session_state["copilot_action"] = "load_dataset"
    with c2:
        if st.button("HTS QC", use_container_width=True):
            st.session_state["copilot_action"] = "hts_qc"

    c3, c4 = st.columns(2)
    with c3:
        if st.button("Dose-response", use_container_width=True):
            st.session_state["copilot_action"] = "dose_response"
    with c4:
        if st.button("Publish", use_container_width=True):
            st.session_state["copilot_action"] = "publish"

    st.markdown("---")

    action = st.selectbox(
        "Action",
        ["load_dataset", "hts_qc", "dose_response", "publish"],
        index=["load_dataset", "hts_qc", "dose_response", "publish"].index(st.session_state["copilot_action"]),
        key="copilot_action_select",
    )
    st.session_state["copilot_action"] = action

    # Entity selectors
    with st.spinner("Loading entities..."):
        opts = _load_entity_options()

    entity_type: str
    entity_id: str
    campaign_id: Optional[str] = None
    compound_id: Optional[str] = None

    if action == "load_dataset":
        entity_type = "dataset"
        labels = [x[0] for x in opts["dataset"]] or ["(no datasets found)"]
        sel = st.selectbox("Dataset", labels, key="copilot_dataset_select")
        entity_id = next((id_ for (lbl, id_) in opts["dataset"] if lbl == sel), "")
    elif action == "hts_qc":
        entity_type = "campaign"
        labels = [x[0] for x in opts["campaign"]] or ["(no campaigns found)"]
        sel = st.selectbox("Campaign", labels, key="copilot_campaign_select")
        entity_id = next((id_ for (lbl, id_) in opts["campaign"] if lbl == sel), "")
        campaign_id = entity_id
    elif action == "dose_response":
        entity_type = "compound"
        labels = [x[0] for x in opts["compound"]] or ["(no compounds found)"]
        sel = st.selectbox("Compound", labels, key="copilot_compound_select")
        entity_id = next((id_ for (lbl, id_) in opts["compound"] if lbl == sel), "")
        compound_id = entity_id
    else:
        entity_type = st.selectbox(
            "Entity Type",
            ["experiment", "dataset", "campaign", "compound"],
            index=0,
            key="copilot_publish_entity_type",
        )
        labels = [x[0] for x in opts.get(entity_type, [])] or [f"(no {entity_type}s found)"]
        sel = st.selectbox("Entity", labels, key="copilot_publish_entity_select")
        entity_id = next((id_ for (lbl, id_) in opts.get(entity_type, []) if lbl == sel), "")
        campaign_id = entity_id if entity_type == "campaign" else None
        compound_id = entity_id if entity_type == "compound" else None

    prompt = st.text_area(
        "Prompt (optional)",
        value="",
        placeholder="Describe what you want the cell to do...",
        height=120,
        key="copilot_prompt",
    )

    tags = st.text_input("Tags (publish)", value="analysis", key="copilot_tags")
    title = st.text_input("Title (publish)", value="Analysis Results", key="copilot_title")

    ctx = AnalysisContext(
        entity_type=entity_type,
        entity_id=entity_id or "",
        campaign_id=campaign_id,
        compound_id=compound_id,
        metadata={},
    )

    if st.button("Generate code", type="primary", use_container_width=True):
        if not ctx.entity_id:
            st.error("No entity selected.")
        else:
            intent = prompt.strip() or _default_intent(action, ctx)
            if action == "publish":
                ctx.metadata = {
                    "tags": [t.strip() for t in tags.split(",") if t.strip()],
                    "title": title,
                }
            with st.spinner("Generating..."):
                try:
                    st.session_state["copilot_generated_code"] = synthesize_cell(intent=intent, context=ctx)
                    st.session_state["copilot_code_explanation"] = ""
                except Exception as e:
                    st.error(f"Copilot failed: {e}")

    code = st.session_state.get("copilot_generated_code", "")
    if code:
        st.markdown("### Generated code")
        st.code(code, language="python")
        if st.button("Copy", use_container_width=True, key="copilot_copy_btn"):
            _copy_code(code)

        if st.button("Explain", use_container_width=True, key="copilot_explain_btn"):
            with st.spinner("Explaining..."):
                try:
                    st.session_state["copilot_code_explanation"] = explain_cell(code)
                except Exception as e:
                    st.error(f"Explain failed: {e}")

        explanation = st.session_state.get("copilot_code_explanation", "")
        if explanation:
            with st.expander("Explanation", expanded=True):
                st.write(explanation)

    return code or None


