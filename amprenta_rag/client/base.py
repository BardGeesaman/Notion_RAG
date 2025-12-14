"""
Shared HTTP client utilities for the Amprenta API client.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Any, Dict, Optional

import httpx


@dataclass(frozen=True)
class APIError(Exception):
    """HTTP error from the API."""

    status_code: int
    message: str
    response_text: str = ""

    def __str__(self) -> str:  # pragma: no cover
        return f"APIError(status_code={self.status_code}, message={self.message})"


class BaseHTTPClient:
    """
    Thin wrapper around httpx.Client with consistent auth + error handling.
    """

    def __init__(
        self,
        api_url: str,
        api_key: Optional[str] = None,
        *,
        timeout: float = 30.0,
        headers: Optional[Dict[str, str]] = None,
    ):
        self.api_url = api_url.rstrip("/")
        self.api_key = api_key

        final_headers: Dict[str, str] = {
            "Accept": "application/json",
        }
        if headers:
            final_headers.update(headers)

        # Optional auth headers (server may ignore if auth is disabled)
        if api_key:
            final_headers.setdefault("Authorization", f"Bearer {api_key}")
            final_headers.setdefault("X-API-Key", api_key)

        self._client = httpx.Client(base_url=self.api_url, timeout=timeout, headers=final_headers)

    def close(self) -> None:
        self._client.close()

    def __enter__(self) -> "BaseHTTPClient":
        return self

    def __exit__(self, exc_type, exc, tb) -> None:
        self.close()

    def _request(
        self,
        method: str,
        path: str,
        *,
        params: Optional[Dict[str, Any]] = None,
        json: Optional[Dict[str, Any]] = None,
    ) -> Any:
        resp = self._client.request(method, path, params=params, json=json)
        if 200 <= resp.status_code < 300:
            # 204 No Content
            if resp.status_code == 204:
                return None
            if not resp.content:
                return None
            return resp.json()

        # Try to extract FastAPI-style {"detail": "..."} message
        message = f"Request failed: {method} {path}"
        response_text = ""
        try:
            data = resp.json()
            if isinstance(data, dict) and "detail" in data:
                message = str(data["detail"])
            else:
                message = str(data)
        except Exception:
            response_text = resp.text or ""

        raise APIError(status_code=resp.status_code, message=message, response_text=response_text)


