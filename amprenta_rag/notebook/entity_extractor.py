"""Entity extraction for Query→Notebook generation.

Uses an LLM to extract entity mentions from a natural language query and then
resolves those mentions against the database to produce concrete UUIDs.
"""

from __future__ import annotations

import json
import os
import re
import uuid
from dataclasses import dataclass
from textwrap import dedent
from typing import Any, Dict, List, Optional, Sequence, Tuple

from amprenta_rag.database.session import db_session
from amprenta_rag.llm.model_registry import call_model
from amprenta_rag.notebook.context import AnalysisContext


DEFAULT_ENTITY_MODEL = os.getenv("AMPRENTA_ENTITY_EXTRACTOR_MODEL", "gpt-4o-mini")


def _looks_like_uuid(value: str) -> bool:
    try:
        uuid.UUID(str(value))
        return True
    except Exception:
        return False


def _normalize_list(x: Any) -> List[str]:
    if x is None:
        return []
    if isinstance(x, list):
        return [str(v).strip() for v in x if str(v).strip()]
    if isinstance(x, str):
        return [x.strip()] if x.strip() else []
    return [str(x).strip()]


def _llm_extract_candidates(query: str) -> Dict[str, Any]:
    """LLM-only: return candidate names/IDs and primary entity hint."""
    system = dedent(
        """\
        You extract entity references from user questions for a scientific platform.
        Return ONLY valid JSON (no markdown) with this schema:
        {
          "datasets": [<uuid-or-name>],
          "experiments": [<uuid-or-name>],
          "compounds": [<uuid-or-name>],
          "campaigns": [<uuid-or-name>],
          "primary": {"type": "dataset|experiment|compound|campaign|null", "value": <uuid-or-name|null>}
        }
        Rules:
        - Prefer UUIDs if user provides them.
        - If user provides names, include the names as strings.
        - If unsure, return empty lists and primary.type=null.
        """
    )
    user = f"Query: {query}"
    raw = call_model(
        DEFAULT_ENTITY_MODEL,
        messages=[{"role": "system", "content": system}, {"role": "user", "content": user}],
        temperature=0.0,
    )
    try:
        return json.loads(raw)
    except Exception:
        # best-effort: attempt to extract a JSON object substring
        m = re.search(r"\{[\s\S]*\}", raw)
        if m:
            return json.loads(m.group(0))
        raise


def _resolve_by_name(db: Any, model: Any, names: Sequence[str]) -> List[str]:
    """Resolve entity names to UUIDs via a best-effort DB lookup."""
    if not names:
        return []
    out: List[str] = []
    for name in names:
        if not name:
            continue
        # UUID passthrough
        if _looks_like_uuid(name):
            out.append(str(uuid.UUID(name)))
            continue
        # Best-effort match by name (case-insensitive)
        q = db.query(model).filter(model.name.ilike(name))  # exact-ish
        row = q.first()
        if row is None:
            q2 = db.query(model).filter(model.name.ilike(f"%{name}%"))
            row = q2.first()
        if row is not None:
            out.append(str(getattr(row, "id")))
    # Deduplicate, preserve order
    seen: set[str] = set()
    deduped: List[str] = []
    for x in out:
        if x not in seen:
            seen.add(x)
            deduped.append(x)
    return deduped


def _pick_primary(resolved: Dict[str, List[str]], primary_hint: Dict[str, Any] | None) -> Tuple[str, str]:
    """Return (entity_type, entity_id). Falls back to first available resolved ID."""
    if primary_hint and isinstance(primary_hint, dict):
        t = (primary_hint.get("type") or "").strip().lower()
        v = primary_hint.get("value")
        if t in {"dataset", "experiment", "compound", "campaign"} and isinstance(v, str) and v.strip():
            # If UUID or name is provided, we may only have UUID in resolved lists.
            # Prefer resolved IDs where possible.
            if _looks_like_uuid(v):
                return t, str(uuid.UUID(v))
            # Try to find the resolved ID that likely corresponds to the named entity
            pool_key = {
                "dataset": "dataset_ids",
                "experiment": "experiment_ids",
                "compound": "compound_ids",
                "campaign": "campaign_ids",
            }[t]
            if resolved.get(pool_key):
                return t, resolved[pool_key][0]

    for t, k in [
        ("dataset", "dataset_ids"),
        ("experiment", "experiment_ids"),
        ("compound", "compound_ids"),
        ("campaign", "campaign_ids"),
    ]:
        if resolved.get(k):
            return t, resolved[k][0]

    # Default: unknown context (still valid for notebook generation)
    return "dataset", "unknown"


def extract_entities(query: str) -> Dict[str, Any]:
    """Extract entity IDs from a query and return an AnalysisContext for notebook generation.

    Returns:
        Dict with keys:
        - dataset_ids, experiment_ids, compound_ids, campaign_ids: lists of UUID strings
        - context: AnalysisContext serialized dict (includes primary entity)
    """
    candidates = _llm_extract_candidates(query)
    ds = _normalize_list(candidates.get("datasets"))
    exps = _normalize_list(candidates.get("experiments"))
    cmps = _normalize_list(candidates.get("compounds"))
    camps = _normalize_list(candidates.get("campaigns"))
    primary_hint = candidates.get("primary") if isinstance(candidates, dict) else None

    resolved: Dict[str, List[str]] = {
        "dataset_ids": [],
        "experiment_ids": [],
        "compound_ids": [],
        "campaign_ids": [],
    }

    # DB resolution (best-effort). Separated so tests can monkeypatch _resolve_by_name.
    from amprenta_rag.database.models import Compound, Dataset, Experiment, HTSCampaign

    with db_session() as db:
        resolved["dataset_ids"] = _resolve_by_name(db, Dataset, ds)
        resolved["experiment_ids"] = _resolve_by_name(db, Experiment, exps)
        resolved["compound_ids"] = _resolve_by_name(db, Compound, cmps)
        resolved["campaign_ids"] = _resolve_by_name(db, HTSCampaign, camps)

    primary_type, primary_id = _pick_primary(resolved, primary_hint if isinstance(primary_hint, dict) else None)

    ctx = AnalysisContext(
        entity_type=primary_type,
        entity_id=primary_id,
        campaign_id=(resolved["campaign_ids"][0] if resolved["campaign_ids"] else None),
        compound_id=(resolved["compound_ids"][0] if resolved["compound_ids"] else None),
        version=1,
        metadata={
            "query": query,
            "candidates": {
                "datasets": ds,
                "experiments": exps,
                "compounds": cmps,
                "campaigns": camps,
                "primary": primary_hint,
            },
            "resolved": resolved,
        },
    )

    return {**resolved, "context": ctx.to_dict()}


__all__ = ["extract_entities"]


