"""Notebook planning for Query→Notebook generation.

Given a user query, an AnalysisContext, and retrieved RAG chunks, this module
creates a lightweight "plan" (cell intents) for a full notebook.
"""

from __future__ import annotations

import json
import os
import re
from textwrap import dedent
from typing import Any, Dict, List, Sequence

from amprenta_rag.llm.model_registry import call_model
from amprenta_rag.notebook.context import AnalysisContext


DEFAULT_PLANNER_MODEL = os.getenv("AMPRENTA_NOTEBOOK_PLANNER_MODEL", "gpt-4o-mini")


def _safe_json_loads(raw: str) -> Any:
    try:
        return json.loads(raw)
    except Exception:
        m = re.search(r"\[[\s\S]*\]", raw)
        if m:
            return json.loads(m.group(0))
        raise


def plan_notebook(query: str, context: AnalysisContext, rag_chunks: Sequence[str]) -> List[Dict[str, str]]:
    """Plan a notebook as a list of cell intents.

    Args:
        query: User question
        context: AnalysisContext (primary entity resolved)
        rag_chunks: Retrieved chunk texts/snippets (strings)

    Returns:
        List of dicts: [{type: "code"|"markdown", intent: str, section: str}, ...]
    """
    # Keep chunks short; this is planning, not synthesis.
    chunk_snips = "\n\n".join((c[:600] for c in rag_chunks[:6]))

    system = dedent(
        """\
        You are a notebook planner. Create a concise outline for a Jupyter notebook.
        Return ONLY valid JSON (no markdown), as a list of objects:
        [{"type":"markdown|code","section":"...","intent":"..."}, ...]

        Rules:
        - Include sections: Context, Data Loading, Analysis, Visualization, Summary
        - Keep intents action-oriented and specific.
        - Do NOT include code; only intents.
        - Keep list length between 6 and 12 cells.
        """
    )

    user = dedent(
        f"""\
        User query:
        {query}

        AnalysisContext (JSON):
        {context.to_json()}

        Retrieved evidence snippets:
        {chunk_snips}
        """
    )

    raw = call_model(
        DEFAULT_PLANNER_MODEL,
        messages=[{"role": "system", "content": system}, {"role": "user", "content": user}],
        temperature=0.2,
    )
    plan = _safe_json_loads(raw)

    out: List[Dict[str, str]] = []
    if isinstance(plan, list):
        for item in plan:
            if not isinstance(item, dict):
                continue
            t = str(item.get("type") or "code").strip().lower()
            if t not in {"code", "markdown"}:
                t = "code"
            section = str(item.get("section") or "").strip() or "Analysis"
            intent = str(item.get("intent") or "").strip()
            if not intent:
                continue
            out.append({"type": t, "section": section, "intent": intent})

    if out:
        return out

    # Fallback plan
    return [
        {"type": "markdown", "section": "Context", "intent": "Describe the question and analysis context."},
        {"type": "code", "section": "Data Loading", "intent": "Load the primary entity using amprenta_rag.notebook helpers and display a summary."},
        {"type": "code", "section": "Analysis", "intent": "Run exploratory analysis relevant to the query and compute key metrics."},
        {"type": "code", "section": "Visualization", "intent": "Create a simple plot/table that supports answering the query."},
        {"type": "markdown", "section": "Summary", "intent": "Summarize findings, caveats, and recommended follow-ups."},
    ]


__all__ = ["plan_notebook"]


