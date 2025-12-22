"""Notebook copilot functions powered by LLM synthesis.

This module provides thin wrappers around `amprenta_rag.llm.model_registry.call_model`
to generate, fix, and explain notebook cells using an `AnalysisContext`.
"""

from __future__ import annotations

import os
from textwrap import dedent

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


