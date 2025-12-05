"""
FastAPI service layer for the multi-omics platform.

This package provides REST API endpoints for accessing the multi-omics
platform data through Postgres (TIER 3 architecture evolution).
"""

from amprenta_rag.api.main import app

__all__ = ["app"]

