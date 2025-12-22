"""Analysis context schema for notebook copilot.

Provides a standardized context object that captures the current analysis
environment (entity being analyzed, related IDs, version info) for use
in notebook cell generation and copilot interactions.
"""

from __future__ import annotations

from dataclasses import asdict, dataclass, field
from datetime import datetime
from typing import Any, Dict, Optional
import json


@dataclass
class AnalysisContext:
    """Context for notebook analysis sessions.

    This context is typically set in the first cell of a notebook and
    provides information about what entity is being analyzed.

    Attributes:
        entity_type: Type of entity ("experiment", "campaign", "dataset", "compound")
        entity_id: UUID of the primary entity being analyzed
        campaign_id: Optional HTS campaign ID for screening workflows
        plate_id: Optional plate ID for plate-level analysis
        compound_id: Optional compound ID for chemistry workflows
        version: Context version for compatibility
        timestamp: ISO format timestamp when context was created
        metadata: Additional key-value pairs for custom context
    """

    entity_type: str
    entity_id: str
    campaign_id: Optional[str] = None
    plate_id: Optional[str] = None
    compound_id: Optional[str] = None
    version: int = 1
    timestamp: str = field(default_factory=lambda: datetime.utcnow().isoformat())
    metadata: Dict[str, Any] = field(default_factory=dict)

    def to_dict(self) -> Dict[str, Any]:
        """Convert context to dictionary."""
        return asdict(self)

    def to_json(self) -> str:
        """Convert context to JSON string."""
        return json.dumps(self.to_dict(), indent=2)

    @classmethod
    def from_dict(cls, data: Dict[str, Any]) -> "AnalysisContext":
        """Create context from dictionary."""
        return cls(
            entity_type=data["entity_type"],
            entity_id=data["entity_id"],
            campaign_id=data.get("campaign_id"),
            plate_id=data.get("plate_id"),
            compound_id=data.get("compound_id"),
            version=data.get("version", 1),
            timestamp=data.get("timestamp", datetime.utcnow().isoformat()),
            metadata=data.get("metadata", {}),
        )

    @classmethod
    def from_json(cls, json_str: str) -> "AnalysisContext":
        """Create context from JSON string."""
        return cls.from_dict(json.loads(json_str))


def generate_context_cell(context: AnalysisContext) -> str:
    """Generate Python code for a notebook cell that sets the analysis context.

    This is typically the first cell in an analysis notebook.

    Args:
        context: The AnalysisContext to embed in the cell

    Returns:
        Python code string for the cell
    """
    return f'''# Analysis Context - DO NOT MODIFY
# This cell defines the context for this analysis session

from amprenta_rag.notebook.context import AnalysisContext

CONTEXT = AnalysisContext.from_dict({context.to_dict()!r})

print(f"Analysis: {{CONTEXT.entity_type}} {{CONTEXT.entity_id}}")
'''


