"""
FastAPI service layer for the multi-omics platform.

This package provides REST API endpoints for accessing the multi-omics
platform data through Postgres (TIER 3 architecture evolution).
"""

try:
    # Importing the full FastAPI app can require optional extras (e.g. python-multipart).
    # Keep the package import-safe so modules like `amprenta_rag.api.schemas` can be imported in minimal envs.
    from amprenta_rag.api.main import app  # type: ignore
except Exception:  # noqa: BLE001
    app = None  # type: ignore[assignment]

__all__ = ["app"]

