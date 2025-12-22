"""IPython magic for the notebook copilot (%%copilot).

This module is intended to be loaded as an IPython extension:

    %load_ext amprenta_rag.notebook.copilot_magic

Example:

    %%copilot
    action: dose_response
    compound_id: abc-123

The cell body is parsed as a YAML-like mapping (simple key/value pairs).
If PyYAML is installed, `yaml.safe_load` is used. Otherwise, a minimal parser
supports the common "key: value" format used in this project.
"""

from __future__ import annotations

from dataclasses import asdict
from typing import Any, Dict, Optional

from amprenta_rag.notebook.context import AnalysisContext
from amprenta_rag.notebook.copilot import synthesize_cell


def _parse_simple_yaml_mapping(text: str) -> Dict[str, Any]:
    """Parse a minimal YAML mapping with `key: value` pairs.

    Supports:
    - Comments starting with '#'
    - Blank lines
    - Unquoted strings, or quoted strings with single/double quotes
    """
    out: Dict[str, Any] = {}
    for raw in text.splitlines():
        line = raw.strip()
        if not line or line.startswith("#"):
            continue
        if ":" not in line:
            raise ValueError(f"Invalid YAML line (expected 'key: value'): {raw!r}")
        key, value = line.split(":", 1)
        key = key.strip()
        value = value.strip()

        if not key:
            raise ValueError(f"Invalid YAML line (empty key): {raw!r}")

        if value == "":
            out[key] = None
            continue

        # Strip quotes if present
        if (value.startswith('"') and value.endswith('"')) or (value.startswith("'") and value.endswith("'")):
            out[key] = value[1:-1]
        else:
            out[key] = value

    return out


def _parse_yaml(cell_text: str) -> Dict[str, Any]:
    """Parse YAML using PyYAML if available, otherwise a minimal parser."""
    try:
        import yaml  # type: ignore

        parsed = yaml.safe_load(cell_text)  # type: ignore[attr-defined]
        if parsed is None:
            return {}
        if not isinstance(parsed, dict):
            raise ValueError("%%copilot YAML must be a mapping (key/value pairs).")
        return dict(parsed)
    except ImportError:
        return _parse_simple_yaml_mapping(cell_text)


def _infer_entity(mapping: Dict[str, Any]) -> tuple[str, str]:
    """Infer (entity_type, entity_id) from common id fields."""
    if mapping.get("entity_type") and mapping.get("entity_id"):
        return str(mapping["entity_type"]), str(mapping["entity_id"])

    if mapping.get("experiment_id"):
        return "experiment", str(mapping["experiment_id"])
    if mapping.get("dataset_id"):
        return "dataset", str(mapping["dataset_id"])
    if mapping.get("campaign_id"):
        return "campaign", str(mapping["campaign_id"])
    if mapping.get("compound_id"):
        return "compound", str(mapping["compound_id"])

    raise ValueError("%%copilot requires an entity id (e.g., experiment_id/dataset_id/campaign_id/compound_id).")


def _intent_from_action(action: str, mapping: Dict[str, Any]) -> str:
    action_norm = action.strip().lower().replace("-", "_")
    if action_norm == "dose_response":
        return f"Fit dose-response for compound_id={mapping.get('compound_id')} and summarize results."
    if action_norm in {"hts_qc", "qc"}:
        return f"Run HTS QC for campaign_id={mapping.get('campaign_id')} and summarize key metrics."
    if action_norm in {"publish", "publish_to_rag"}:
        return "Publish the current analysis results to RAG with tags and a title."
    return f"Generate a notebook cell for action '{action}'."


def load_ipython_extension(ipython) -> None:
    """IPython extension entrypoint."""
    from IPython.core.magic import register_cell_magic

    @register_cell_magic
    def copilot(line: str, cell: str) -> str:  # noqa: ANN001
        """%%copilot cell magic: parse YAML config and synthesize a cell."""
        mapping = _parse_yaml(cell)

        action = str(mapping.get("action") or mapping.get("intent") or "").strip()
        if not action:
            raise ValueError("%%copilot requires an 'action' (or 'intent') field.")

        entity_type, entity_id = _infer_entity(mapping)

        context = AnalysisContext(
            entity_type=entity_type,
            entity_id=entity_id,
            campaign_id=(str(mapping["campaign_id"]) if mapping.get("campaign_id") else None),
            plate_id=(str(mapping["plate_id"]) if mapping.get("plate_id") else None),
            compound_id=(str(mapping["compound_id"]) if mapping.get("compound_id") else None),
            metadata={k: v for k, v in mapping.items() if k not in {"action", "intent"}},
        )

        intent = _intent_from_action(action, mapping)
        code = synthesize_cell(intent=intent, context=context)

        # Save for debugging / reuse
        ipython.user_ns["_copilot_last_request"] = {"intent": intent, "context": asdict(context)}
        ipython.user_ns["_copilot_last_code"] = code

        print(code)
        return code


