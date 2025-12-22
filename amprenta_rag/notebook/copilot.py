"""Notebook copilot functions powered by LLM synthesis.

This module provides thin wrappers around `amprenta_rag.llm.model_registry.call_model`
to generate, fix, and explain notebook cells using an `AnalysisContext`.
"""

from __future__ import annotations

import json
import os
from textwrap import dedent
from typing import Any, Dict, List, Optional, Sequence, Union

from amprenta_rag.llm.model_registry import call_model
from amprenta_rag.notebook.context import AnalysisContext


DEFAULT_NOTEBOOK_MODEL = os.getenv("AMPRENTA_NOTEBOOK_MODEL", "gpt-4o-mini")


def synthesize_cell(intent: str, context: AnalysisContext) -> str:
    """Synthesize a single notebook cell (Python code) for the given intent and context."""
    system = dedent(
        """\
        You are a notebook copilot for the Amprenta codebase.
        Generate ONE Python notebook cell as plain code (no markdown, no backticks).
        Prefer using vetted helpers from `amprenta_rag.notebook` when possible.
        Keep imports minimal and at the top of the cell.
        Do not access the network. Do not read local files. Do not use secrets.
        """
    )

    user = dedent(
        f"""\
        Intent: {intent}

        AnalysisContext (JSON):
        {context.to_json()}

        Requirements:
        - Return ONLY Python code for a single cell.
        - If you need the context, assume `CONTEXT` exists as an AnalysisContext instance.
        - Use helper functions like: load_dataset/load_experiment/load_campaign/load_compound,
          run_hts_qc, fit_dose_response, publish_to_rag.
        """
    )

    return call_model(
        DEFAULT_NOTEBOOK_MODEL,
        messages=[
            {"role": "system", "content": system},
            {"role": "user", "content": user},
        ],
        temperature=0.2,
    )


def fix_cell(code: str, error: str) -> str:
    """Fix a notebook cell given the failing code and error output."""
    system = dedent(
        """\
        You are a senior Python engineer fixing a single Jupyter cell.
        Return ONLY the corrected Python code for the cell (no markdown, no backticks).
        Keep the change minimal, and preserve the intent of the original cell.
        Prefer using `amprenta_rag.notebook` helpers where appropriate.
        """
    )

    user = dedent(
        f"""\
        The following cell failed:

        {code}

        Error:
        {error}

        Return only the corrected cell code.
        """
    )

    return call_model(
        DEFAULT_NOTEBOOK_MODEL,
        messages=[
            {"role": "system", "content": system},
            {"role": "user", "content": user},
        ],
        temperature=0.1,
    )


def explain_cell(code: str) -> str:
    """Explain what a notebook cell does."""
    system = dedent(
        """\
        You explain notebook cells concisely for scientists.
        Return a short plain-text explanation (no markdown code fences).
        """
    )

    user = dedent(
        f"""\
        Explain what this cell does:

        {code}
        """
    )

    return call_model(
        DEFAULT_NOTEBOOK_MODEL,
        messages=[
            {"role": "system", "content": system},
            {"role": "user", "content": user},
        ],
        temperature=0.2,
    )


def _coerce_cell_source(src: Any) -> str:
    if src is None:
        return ""
    if isinstance(src, str):
        return src
    if isinstance(src, list):
        return "".join(str(x) for x in src)
    return str(src)


def _try_parse_json_dict(text: str) -> Optional[Dict[str, Any]]:
    """Parse a JSON object from model output; tolerate extra text."""
    s = (text or "").strip()
    if not s:
        return None
    try:
        obj = json.loads(s)
        return obj if isinstance(obj, dict) else None
    except Exception:
        pass

    # Try to salvage by extracting the first {...} block.
    try:
        start = s.find("{")
        end = s.rfind("}")
        if start != -1 and end != -1 and end > start:
            obj2 = json.loads(s[start : end + 1])
            return obj2 if isinstance(obj2, dict) else None
    except Exception:
        return None
    return None


def summarize_notebook(notebook: Union[Dict[str, Any], Any]) -> Dict[str, Any]:
    """Summarize an nbformat notebook (dict or NotebookNode) into structured fields."""
    nb: Dict[str, Any] = notebook if isinstance(notebook, dict) else dict(notebook)  # type: ignore[arg-type]
    cells: Sequence[Dict[str, Any]] = nb.get("cells", []) or []
    cell_count = len(cells)

    extracted: List[str] = []
    for i, cell in enumerate(cells, start=1):
        ctype = str(cell.get("cell_type") or "")
        if ctype not in {"markdown", "code"}:
            continue
        src = _coerce_cell_source(cell.get("source"))
        if not src.strip():
            continue
        # Keep context bounded.
        src = src.strip()
        if len(src) > 2000:
            src = src[:2000] + "\n... (truncated)"
        extracted.append(f"[{i}:{ctype}]\n{src}")

    # If there's nothing to summarize, return a deterministic fallback without calling the LLM.
    if cell_count == 0 or not extracted:
        return {
            "title": "Notebook Summary",
            "entity_summary": "",
            "methods": "",
            "key_findings": "",
            "cell_count": cell_count,
        }

    context = "\n\n".join(extracted)
    # Guard against overly large prompts.
    if len(context) > 12000:
        context = context[:12000] + "\n... (truncated)"

    system = dedent(
        """\
        You summarize scientific Jupyter notebooks concisely for scientists.
        Return ONLY a JSON object (no markdown, no code fences).
        Required keys: title, entity_summary, methods, key_findings.
        Values should be plain strings (use bullet-like newlines if helpful).
        """
    )

    user = dedent(
        f"""\
        Summarize this notebook.

        Notebook cell count: {cell_count}

        Extracted cells (markdown + code):
        {context}

        Return JSON with:
        - title: short descriptive title
        - entity_summary: what entity/entities this notebook analyzes
        - methods: key methods/steps used
        - key_findings: main results/conclusions
        """
    )

    raw = call_model(
        DEFAULT_NOTEBOOK_MODEL,
        messages=[
            {"role": "system", "content": system},
            {"role": "user", "content": user},
        ],
        temperature=0.3,
    )

    parsed = _try_parse_json_dict(raw)
    if not parsed:
        return {
            "title": "Notebook Summary",
            "entity_summary": "",
            "methods": "",
            "key_findings": (raw or "").strip(),
            "cell_count": cell_count,
        }

    return {
        "title": str(parsed.get("title") or "Notebook Summary"),
        "entity_summary": str(parsed.get("entity_summary") or ""),
        "methods": str(parsed.get("methods") or ""),
        "key_findings": str(parsed.get("key_findings") or ""),
        "cell_count": cell_count,
    }


__all__ = [
    "AnalysisContext",
    "synthesize_cell",
    "fix_cell",
    "explain_cell",
    "summarize_notebook",
]

