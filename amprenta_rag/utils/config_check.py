"""Configuration and secrets validation utilities."""

import os
import logging
from typing import List, Tuple

logger = logging.getLogger(__name__)


def validate_required_secrets() -> Tuple[bool, List[str]]:
    """Validate all required secrets are present at startup.
    
    Returns:
        Tuple of (all_valid, list_of_missing_secrets)
    """
    required = [
        ("SIGNATURE_SECRET_KEY", "Electronic signature HMAC key"),
        ("JWT_SECRET_KEY", "JWT authentication signing key"),
    ]
    
    # Database password required unless using full URL
    if not os.getenv("DATABASE_URL") and not os.getenv("POSTGRES_URL"):
        required.append(("POSTGRES_PASSWORD", "Database password"))
    
    missing = []
    for var_name, description in required:
        if not os.getenv(var_name):
            missing.append(f"{var_name} ({description})")
            logger.error(f"Missing required secret: {var_name}")
    
    # Warnings for recommended secrets
    recommended = [
        ("OPENAI_API_KEY", "LLM features disabled"),
        ("NCBI_EMAIL", "NCBI/GEO queries may fail"),
    ]
    
    for var_name, consequence in recommended:
        if not os.getenv(var_name):
            logger.warning(f"Recommended secret missing: {var_name} - {consequence}")
    
    return len(missing) == 0, missing


def validate_config():
    """Legacy function for backward compatibility."""
    pass