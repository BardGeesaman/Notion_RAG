"""Structured LLM entity extractor for parsed documents."""

from __future__ import annotations

import json
import logging
import os
import re
from typing import Any, Dict, List

from pydantic import BaseModel, Field, ValidationError

from amprenta_rag.clients.openai_client import get_default_models, get_openai_client

logger = logging.getLogger(__name__)


class ExtractedCompound(BaseModel):
    name: str
    identifiers: List[str] = Field(default_factory=list)
    properties: Dict[str, Any] = Field(default_factory=dict)


class ExtractedGene(BaseModel):
    symbol: str
    aliases: List[str] = Field(default_factory=list)
    organism: str | None = None


class ExtractedDisease(BaseModel):
    name: str
    identifiers: List[str] = Field(default_factory=list)


class ExtractedEntity(BaseModel):
    entity_type: str
    name: str
    metadata: Dict[str, Any] = Field(default_factory=dict)


class ExtractionResult(BaseModel):
    compounds: List[ExtractedCompound] = Field(default_factory=list)
    genes: List[ExtractedGene] = Field(default_factory=list)
    diseases: List[ExtractedDisease] = Field(default_factory=list)
    other: List[ExtractedEntity] = Field(default_factory=list)


def _empty_result() -> ExtractionResult:
    return ExtractionResult(compounds=[], genes=[], diseases=[], other=[])


def _extract_json_object(raw: str) -> Dict[str, Any]:
    try:
        return json.loads(raw)
    except Exception:
        # Best-effort: extract first JSON object substring
        m = re.search(r"\{[\s\S]*\}", raw)
        if not m:
            raise
        return json.loads(m.group(0))


def extract_entities(text: str, doc_type: str = "generic") -> ExtractionResult:
    """
    Extract entities from text using an LLM and validate against ExtractionResult.

    Retries up to 3x on ValidationError with error feedback.
    On failure, returns an empty ExtractionResult (and logs the error).
    """
    if not text or not text.strip():
        return _empty_result()

    cfg_model, _ = get_default_models()
    model = os.getenv("AMPRENTA_STRUCTURED_EXTRACTOR_MODEL", cfg_model or "gpt-4o-mini")

    schema = ExtractionResult.model_json_schema()
    system = (
        "You extract structured biomedical entities from documents.\n"
        "Return ONLY valid JSON (no markdown) that conforms to the provided JSON schema.\n"
        "Rules:\n"
        "- Be conservative: only include entities explicitly mentioned.\n"
        "- Prefer normalized names/symbols.\n"
        "- If unsure, leave lists empty.\n"
    )

    base_user = {
        "doc_type": doc_type,
        "schema": schema,
        "text": text,
    }

    messages: List[Dict[str, str]] = [
        {"role": "system", "content": system},
        {
            "role": "user",
            "content": (
                "Extract entities from the following document.\n"
                "Respond with JSON only.\n\n"
                f"{json.dumps(base_user, ensure_ascii=False)}"
            ),
        },
    ]

    client = get_openai_client()
    last_err: Exception | None = None

    for attempt in range(1, 4):
        try:
            resp = client.chat.completions.create(
                model=model,
                messages=messages,
                temperature=0.0,
            )
            raw = (resp.choices[0].message.content or "").strip()
            data = _extract_json_object(raw)
            return ExtractionResult.model_validate(data)
        except ValidationError as ve:
            last_err = ve
            messages.append(
                {
                    "role": "user",
                    "content": (
                        "Your previous JSON did not validate against the schema.\n"
                        f"Validation error:\n{ve}\n\n"
                        "Return corrected JSON ONLY."
                    ),
                }
            )
            continue
        except Exception as e:  # noqa: BLE001
            last_err = e
            messages.append(
                {
                    "role": "user",
                    "content": (
                        "The previous response could not be parsed as JSON.\n"
                        f"Error:\n{e}\n\n"
                        "Return JSON ONLY that conforms to the schema."
                    ),
                }
            )
            continue

    logger.exception("Structured extraction failed after retries (doc_type=%s): %r", doc_type, last_err)
    return _empty_result()


def _flatten_parsed(parsed_doc: dict) -> str:
    parts: List[str] = []

    txt = parsed_doc.get("text")
    if isinstance(txt, str) and txt.strip():
        parts.append(txt.strip())

    # Tables from docx parser
    tables = parsed_doc.get("tables")
    if isinstance(tables, list):
        for t in tables:
            if not isinstance(t, dict):
                continue
            rows = t.get("rows")
            if not isinstance(rows, list):
                continue
            for row in rows:
                if isinstance(row, list):
                    parts.append("\t".join(str(c) for c in row))

    # Slides from pptx parser
    slides = parsed_doc.get("slides")
    if isinstance(slides, list):
        for s in slides:
            if not isinstance(s, dict):
                continue
            for k in ("title", "content", "notes"):
                v = s.get(k)
                if isinstance(v, str) and v.strip():
                    parts.append(v.strip())

    # Sheets from excel parser
    sheets = parsed_doc.get("sheets")
    if isinstance(sheets, list):
        for sh in sheets:
            if not isinstance(sh, dict):
                continue
            name = sh.get("name")
            if isinstance(name, str) and name.strip():
                parts.append(f"[Sheet] {name.strip()}")
            cols = sh.get("columns")
            if isinstance(cols, list) and cols:
                parts.append("Columns: " + ", ".join(str(c) for c in cols))
            sample_rows = sh.get("sample_rows")
            if isinstance(sample_rows, list) and sample_rows:
                parts.append("Sample rows: " + json.dumps(sample_rows[:10], ensure_ascii=False))

    # CSV parser output
    cols = parsed_doc.get("columns")
    if isinstance(cols, list) and cols:
        parts.append("Columns: " + ", ".join(str(c) for c in cols))
    sample = parsed_doc.get("sample")
    if isinstance(sample, list) and sample:
        parts.append("Sample rows: " + json.dumps(sample[:10], ensure_ascii=False))

    return "\n\n".join(parts).strip()


def extract_from_parsed(parsed_doc: dict, doc_type: str) -> ExtractionResult:
    """
    Combine extracted content from a parsed document dict and run structured extraction.
    """
    text = _flatten_parsed(parsed_doc or {})
    return extract_entities(text, doc_type=doc_type)


__all__ = [
    "ExtractedCompound",
    "ExtractedGene",
    "ExtractedDisease",
    "ExtractedEntity",
    "ExtractionResult",
    "extract_entities",
    "extract_from_parsed",
]


