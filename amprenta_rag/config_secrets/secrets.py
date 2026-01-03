"""Unified secrets management with priority: AWS -> env vars -> dev defaults."""

import os
import logging
from functools import lru_cache
from typing import Optional

logger = logging.getLogger(__name__)

# Dev-only defaults (NEVER use in production)
_DEV_DEFAULTS = {
    "SIGNATURE_SECRET_KEY": "dev-signature-key-not-for-production",
    "JWT_SECRET_KEY": "dev-jwt-key-not-for-production",
}


def _get_from_aws(name: str) -> Optional[str]:
    """Attempt to get secret from AWS Secrets Manager."""
    try:
        import boto3
        from botocore.exceptions import ClientError
        
        client = boto3.client("secretsmanager")
        env = os.getenv("ENVIRONMENT", "dev")
        aws_name = f"amprenta/{env}/{name.lower()}"
        
        response = client.get_secret_value(SecretId=aws_name)
        return response.get("SecretString")
    except Exception as e:
        logger.debug(f"AWS Secrets Manager not available for {name}: {e}")
        return None


@lru_cache(maxsize=128)
def get_secret(name: str, default: Optional[str] = None) -> Optional[str]:
    """Get secret with priority: AWS -> env var -> dev default -> default.
    
    Priority order:
    1. AWS Secrets Manager (production)
    2. Environment variable (fallback)
    3. Dev defaults when DISABLE_AUTH=true (last resort, dev/test only)
    4. Provided default parameter
    
    Args:
        name: Secret name (e.g., "JWT_SECRET_KEY")
        default: Default value if not found anywhere
        
    Returns:
        Secret value or default
    """
    # Priority 1: AWS Secrets Manager (production)
    value = _get_from_aws(name)
    if value:
        return value
    
    # Priority 2: Environment variable
    value = os.getenv(name)
    if value:
        return value
    
    # Priority 3: Dev defaults (ONLY when DISABLE_AUTH=true)
    if os.getenv("DISABLE_AUTH", "").lower() in ("1", "true"):
        if name in _DEV_DEFAULTS:
            logger.warning(f"Using dev default for {name} - NOT FOR PRODUCTION")
            return _DEV_DEFAULTS[name]
    
    # Priority 4: Provided default
    return default


def clear_secret_cache() -> None:
    """Clear the secret cache. Useful for testing."""
    get_secret.cache_clear()
