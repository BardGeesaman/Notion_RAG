"""
amprenta_rag.client

Lightweight Python client for the Amprenta Multi-Omics Platform API.
"""

from __future__ import annotations

from typing import Optional

from amprenta_rag.client.base import APIError, BaseHTTPClient
from amprenta_rag.client.compounds import CompoundsClient
from amprenta_rag.client.datasets import DatasetsClient
from amprenta_rag.client.experiments import ExperimentsClient
from amprenta_rag.client.features import FeaturesClient
from amprenta_rag.client.programs import ProgramsClient
from amprenta_rag.client.screening import ScreeningClient
from amprenta_rag.client.signatures import SignaturesClient


class RAGClient:
    """
    Top-level API client.

    Example:
        client = RAGClient(api_url="http://localhost:8000", api_key="...")
        datasets = client.datasets.list(limit=10)
    """

    def __init__(self, api_url: str, api_key: Optional[str] = None):
        self.http = BaseHTTPClient(api_url=api_url, api_key=api_key)
        self.programs = ProgramsClient(self.http)
        self.datasets = DatasetsClient(self.http)
        self.experiments = ExperimentsClient(self.http)
        self.features = FeaturesClient(self.http)
        self.signatures = SignaturesClient(self.http)
        self.compounds = CompoundsClient(self.http)
        self.screening = ScreeningClient(self.http)

    def close(self) -> None:
        self.http.close()

    def __enter__(self) -> "RAGClient":
        return self

    def __exit__(self, exc_type, exc, tb) -> None:
        self.close()


__all__ = [
    "APIError",
    "RAGClient",
    "CompoundsClient",
    "DatasetsClient",
    "ExperimentsClient",
    "FeaturesClient",
    "ProgramsClient",
    "ScreeningClient",
    "SignaturesClient",
]


