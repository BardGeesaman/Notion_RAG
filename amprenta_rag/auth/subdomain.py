"""Subdomain parsing utilities for tenant routing."""

from __future__ import annotations

from typing import Optional

from fastapi import Request


def extract_subdomain(request: Request) -> Optional[str]:
    """Extract subdomain from Host header.

    Examples:
      - "acme.amprenta.app" -> "acme"
      - "localhost:8000" -> None
      - "127.0.0.1:8000" -> None
    """
    host = request.headers.get("host") or ""
    host = host.split(":", 1)[0].strip().lower()
    if not host or host in ("localhost",):
        return None
    if all(part.isdigit() for part in host.split(".")):
        return None
    parts = host.split(".")
    if len(parts) < 3:
        return None
    sub = parts[0].strip()
    return sub or None


__all__ = ["extract_subdomain"]


