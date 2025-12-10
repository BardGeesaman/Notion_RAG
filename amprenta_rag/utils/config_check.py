"""Lightweight configuration validation helpers."""
from __future__ import annotations

import os

REQUIRED_VARS = [
    "DATABASE_URL",
    "OPENAI_API_KEY",
    "PINECONE_API_KEY",
]


def validate_config() -> bool:
    """Validate presence of required environment variables.

    Prints a status line for each variable and returns True if all are set.
    """
    all_ok = True
    for var in REQUIRED_VARS:
        val = os.getenv(var)
        if val:
            print(f"✅ {var} is set")
        else:
            print(f"❌ {var} is NOT set")
            all_ok = False
    return all_ok


__all__ = ["validate_config", "REQUIRED_VARS"]
