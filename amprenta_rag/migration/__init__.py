"""
Migration utilities for TIER 3 architecture evolution.

Provides dual-write capability to write to both Notion and Postgres
during the transition period. New data should be ingested into both systems
simultaneously until Postgres becomes the source of truth.
"""

from amprenta_rag.migration.dual_write import (
    DualWriteManager,
)

__all__ = [
    "DualWriteManager",
]
