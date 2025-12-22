"""Notebook copilot API endpoints.

Provides LLM-backed endpoints to generate, fix, and explain notebook cells,
plus a simple templates listing for clients.
"""

from __future__ import annotations

from typing import Any, Dict, List, Optional

from fastapi import APIRouter, HTTPException
from pydantic import BaseModel, Field

from amprenta_rag.notebook.context import AnalysisContext
from amprenta_rag.notebook.copilot import explain_cell, fix_cell, summarize_notebook, synthesize_cell
from amprenta_rag.notebook.query_to_notebook import default_filename, generate_notebook

router = APIRouter()


class CopilotGenerateRequest(BaseModel):
    intent: str = Field(..., min_length=1)
    context: Dict[str, Any]


class CopilotGenerateResponse(BaseModel):
    code: str


class CopilotFixRequest(BaseModel):
    code: str = Field(..., min_length=1)
    error: str = Field(..., min_length=1)


class CopilotFixResponse(BaseModel):
    code: str


class CopilotExplainRequest(BaseModel):
    code: str = Field(..., min_length=1)


class CopilotExplainResponse(BaseModel):
    explanation: str


class NotebookTemplate(BaseModel):
    id: str
    title: str
    description: str
    required_fields: List[str]
    optional_fields: List[str] = []


@router.post("/copilot/generate", response_model=CopilotGenerateResponse)
def copilot_generate(payload: CopilotGenerateRequest) -> CopilotGenerateResponse:
    """Generate a notebook cell for a given intent and analysis context."""
    try:
        ctx = AnalysisContext.from_dict(payload.context)
    except Exception as e:
        raise HTTPException(status_code=422, detail=f"Invalid context: {e}")

    try:
        code = synthesize_cell(intent=payload.intent, context=ctx)
    except Exception as e:
        raise HTTPException(status_code=500, detail=f"Copilot generation failed: {e}")

    return CopilotGenerateResponse(code=code)


@router.post("/copilot/fix", response_model=CopilotFixResponse)
def copilot_fix(payload: CopilotFixRequest) -> CopilotFixResponse:
    """Fix a failing notebook cell based on the error output."""
    try:
        code = fix_cell(code=payload.code, error=payload.error)
    except Exception as e:
        raise HTTPException(status_code=500, detail=f"Copilot fix failed: {e}")
    return CopilotFixResponse(code=code)


@router.post("/copilot/explain", response_model=CopilotExplainResponse)
def copilot_explain(payload: CopilotExplainRequest) -> CopilotExplainResponse:
    """Explain what a notebook cell does."""
    try:
        explanation = explain_cell(code=payload.code)
    except Exception as e:
        raise HTTPException(status_code=500, detail=f"Copilot explain failed: {e}")
    return CopilotExplainResponse(explanation=explanation)


@router.get("/templates", response_model=List[NotebookTemplate])
def list_templates(entity_type: Optional[str] = None) -> List[NotebookTemplate]:
    """List available cell templates for clients."""
    templates: List[NotebookTemplate] = [
        NotebookTemplate(
            id="context",
            title="Analysis Context",
            description="Sets the AnalysisContext for this notebook session.",
            required_fields=["entity_type", "entity_id"],
            optional_fields=["campaign_id", "plate_id", "compound_id", "metadata"],
        ),
        NotebookTemplate(
            id="hts_qc",
            title="HTS QC",
            description="Runs quality control metrics for an HTS campaign.",
            required_fields=["campaign_id"],
        ),
        NotebookTemplate(
            id="dose_response",
            title="Dose Response",
            description="Fits a dose-response curve for a compound (optionally filtered by campaign).",
            required_fields=["compound_id"],
            optional_fields=["campaign_id"],
        ),
        NotebookTemplate(
            id="publish",
            title="Publish to RAG",
            description="Publishes analysis results to the RAG index (MVP logs publish).",
            required_fields=["data", "tags"],
            optional_fields=["title", "entity_type", "entity_id"],
        ),
    ]

    if entity_type:
        et = entity_type.strip().lower()
        if et == "campaign":
            return [t for t in templates if t.id in {"context", "hts_qc"}]
        if et == "compound":
            return [t for t in templates if t.id in {"context", "dose_response"}]
        if et in {"dataset", "experiment"}:
            return [t for t in templates if t.id in {"context", "publish"}]

    return templates


class GenerateNotebookFromQueryRequest(BaseModel):
    query: str = Field(..., min_length=1)


class GenerateNotebookFromQueryResponse(BaseModel):
    notebook_json: Dict[str, Any]
    filename: str


@router.post("/notebook/generate-from-query", response_model=GenerateNotebookFromQueryResponse)
def generate_from_query(payload: GenerateNotebookFromQueryRequest) -> GenerateNotebookFromQueryResponse:
    """Generate a full Jupyter notebook from a natural language query."""
    try:
        nb = generate_notebook(payload.query)
        # NotebookNode is dict-like but not always plain-JSON; normalize via nbformat serialization.
        import nbformat as _nbformat
        import json as _json

        nb_json = _json.loads(_nbformat.writes(nb))
        fname = default_filename(payload.query)
        return GenerateNotebookFromQueryResponse(notebook_json=nb_json, filename=fname)
    except HTTPException:
        raise
    except Exception as e:
        raise HTTPException(status_code=500, detail=f"Notebook generation failed: {e}")


class SummarizeNotebookRequest(BaseModel):
    notebook_json: Dict[str, Any]


class SummarizeNotebookResponse(BaseModel):
    title: str
    entity_summary: str
    methods: str
    key_findings: str
    cell_count: int


@router.post("/notebook/summarize", response_model=SummarizeNotebookResponse)
def summarize(payload: SummarizeNotebookRequest) -> SummarizeNotebookResponse:
    """Summarize a Jupyter notebook JSON payload into a structured summary."""
    try:
        import nbformat as _nbformat

        nb = _nbformat.from_dict(payload.notebook_json)
        _nbformat.validate(nb)
    except Exception as e:
        raise HTTPException(status_code=422, detail=f"Invalid notebook: {e}")

    try:
        summary = summarize_notebook(nb)
        return SummarizeNotebookResponse(**summary)
    except HTTPException:
        raise
    except Exception as e:
        raise HTTPException(status_code=500, detail=f"Notebook summarization failed: {e}")


