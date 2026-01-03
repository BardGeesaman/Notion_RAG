"""Configuration module."""
from amprenta_rag.config_secrets.secrets import get_secret, clear_secret_cache

__all__ = ["get_secret", "clear_secret_cache"]
