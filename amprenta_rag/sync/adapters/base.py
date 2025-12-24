from __future__ import annotations

from abc import ABC, abstractmethod
from datetime import datetime
from typing import AsyncIterator
from uuid import UUID


class BaseSyncAdapter(ABC):
    """Base adapter for synchronizing records from an external data source."""

    source: str  # "chembl", "pubchem", etc.

    @abstractmethod
    async def fetch_records(self, since: datetime | None) -> AsyncIterator[dict]:
        """Fetch records from external source, optionally since last sync."""

    @abstractmethod
    def compute_checksum(self, record: dict) -> str:
        """Compute MD5 hash for change detection."""

    @abstractmethod
    def map_to_entity(self, record: dict, db_session) -> tuple[str, UUID | None]:
        """Map external record to local entity. Returns (entity_type, entity_id)."""


