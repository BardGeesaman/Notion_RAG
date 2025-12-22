"""Query→Notebook generator.

Orchestrates:
1) entity extraction
2) RAG retrieval for evidence snippets
3) notebook planning (cell intents)
4) cell synthesis (code) and assembly into an nbformat notebook
"""

from __future__ import annotations

import datetime as dt
import re
from typing import Any, Dict, List, Optional, Sequence

import nbformat

from amprenta_rag.notebook.context import AnalysisContext, generate_context_cell
from amprenta_rag.notebook.copilot import synthesize_cell
from amprenta_rag.notebook.entity_extractor import extract_entities
from amprenta_rag.notebook.notebook_planner import plan_notebook
from amprenta_rag.query.rag.query import query_rag


def _slugify(s: str, max_len: int = 60) -> str:
    s2 = re.sub(r"[^a-zA-Z0-9]+", "-", s.strip()).strip("-").lower()
    if not s2:
        return "notebook"
    return s2[:max_len].strip("-")


def generate_notebook(query: str) -> nbformat.NotebookNode:
    """Generate a full notebook from a natural language query."""
    extracted = extract_entities(query)
    ctx = AnalysisContext.from_dict(extracted.get("context") or {"entity_type": "dataset", "entity_id": "unknown"})

    # Retrieve evidence but avoid synthesis (no OpenAI) at this stage.
    rag_result = query_rag(
        query,
        top_k=10,
        generate_answer=False,
        use_cache=False,
    )
    rag_chunks: List[str] = []
    try:
        rag_chunks = list(getattr(rag_result, "context_chunks", []) or [])
    except Exception:
        rag_chunks = []

    plan = plan_notebook(query=query, context=ctx, rag_chunks=rag_chunks)

    nb = nbformat.v4.new_notebook()
    nb["metadata"] = {
        "kernelspec": {"name": "python3", "display_name": "Python 3"},
        "language_info": {"name": "python", "pygments_lexer": "python"},
        "amprenta": {
            "generated_by": "query_to_notebook",
            "generated_at": dt.datetime.utcnow().isoformat(),
        },
    }

    cells: List[nbformat.NotebookNode] = []

    # Always include the context cell first.
    cells.append(nbformat.v4.new_code_cell(generate_context_cell(ctx)))

    # Add a title cell
    cells.append(
        nbformat.v4.new_markdown_cell(
            f"## Query Notebook\n\n**Question:** {query}\n"
        )
    )

    last_section: Optional[str] = None
    for item in plan:
        cell_type = (item.get("type") or "code").strip().lower()
        section = (item.get("section") or "Analysis").strip()
        intent = (item.get("intent") or "").strip()
        if not intent:
            continue

        if section and section != last_section:
            cells.append(nbformat.v4.new_markdown_cell(f"### {section}\n"))
            last_section = section

        if cell_type == "markdown":
            # Keep markdown deterministic; avoid an extra LLM call.
            cells.append(nbformat.v4.new_markdown_cell(intent))
        else:
            code = synthesize_cell(intent=intent, context=ctx)
            cells.append(nbformat.v4.new_code_cell(code))

    nb["cells"] = cells
    return nb


def default_filename(query: str) -> str:
    """Suggested filename for a generated notebook."""
    ts = dt.datetime.utcnow().strftime("%Y%m%d_%H%M%S")
    return f"{_slugify(query, max_len=50)}_{ts}.ipynb"


__all__ = ["generate_notebook", "default_filename"]


