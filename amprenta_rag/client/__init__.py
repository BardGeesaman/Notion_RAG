"""
amprenta_rag.client

Lightweight Python client for the Amprenta Multi-Omics Platform API.
"""

from __future__ import annotations

from typing import Optional

from amprenta_rag.client.base import APIError, BaseHTTPClient
from amprenta_rag.client.datasets import DatasetsClient
from amprenta_rag.client.experiments import ExperimentsClient


class RAGClient:
    """
    Top-level API client.

    Example:
        client = RAGClient(api_url="http://localhost:8000", api_key="...")
        datasets = client.datasets.list(limit=10)
    """

    def __init__(self, api_url: str, api_key: Optional[str] = None):
        self.http = BaseHTTPClient(api_url=api_url, api_key=api_key)
        self.datasets = DatasetsClient(self.http)
        self.experiments = ExperimentsClient(self.http)

    def close(self) -> None:
        self.http.close()

    def __enter__(self) -> "RAGClient":
        return self

    def __exit__(self, exc_type, exc, tb) -> None:
        self.close()


__all__ = [
    "APIError",
    "RAGClient",
    "DatasetsClient",
    "ExperimentsClient",
]


