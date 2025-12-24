"""External source synchronization framework (MVP).

This package provides a `SyncManager` orchestrator and adapter interfaces for
integrating external sources (e.g., ChEMBL, PubChem) into the local database.
"""

from amprenta_rag.sync.manager import SyncManager
from amprenta_rag.sync.adapters.base import BaseSyncAdapter

__all__ = ["SyncManager", "BaseSyncAdapter"]


