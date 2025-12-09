"""
Configuration validation helpers.

Provides utilities to validate configuration settings and provide
helpful error messages for common configuration issues.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import List, Optional

from amprenta_rag.config import get_config
from amprenta_rag.logging_utils import get_logger

logger = get_logger(__name__)


class ConfigValidationError(Exception):
    """Raised when configuration validation fails."""
    pass


# Legacy API expected by utils/__init__.py
class ValidationError(ConfigValidationError):
    """Alias for backward compatibility."""
    pass


@dataclass
class ValidationResult:
    success: bool
    errors: List[str]


class ConfigValidator:
    """Lightweight validator wrapper for backward compatibility."""

    def __init__(self):
        self.errors: List[str] = []

    def validate(self) -> ValidationResult:
        issues = validate_all_config()
        flat_errors: List[str] = []
        for category, errs in issues.items():
            for err in errs:
                flat_errors.append(f"{category}: {err}")
        self.errors = flat_errors
        return ValidationResult(success=not flat_errors, errors=flat_errors)


def validate_configuration() -> ValidationResult:
    """Run all validators and return ValidationResult."""
    validator = ConfigValidator()
    return validator.validate()


def validate_postgres_config() -> List[str]:
    """
    Validate Postgres configuration and return list of issues.
    
    Returns:
        List of validation error messages (empty if valid)
    """
    issues: List[str] = []
    cfg = get_config()
    
    if cfg.pipeline.use_postgres_as_sot:
        # Check if Postgres connection details are provided
        if not cfg.database.postgres_host:
            issues.append("POSTGRES_HOST is not set but USE_POSTGRES_AS_SOT=true")
        
        if not cfg.database.postgres_db:
            issues.append("POSTGRES_DB is not set but USE_POSTGRES_AS_SOT=true")
        
        if not cfg.database.postgres_user:
            issues.append("POSTGRES_USER is not set but USE_POSTGRES_AS_SOT=true")
    
    return issues


def validate_feature_linking_config() -> List[str]:
    """
    Validate feature linking configuration and return list of issues.
    
    Returns:
        List of validation error messages (empty if valid)
    """
    issues: List[str] = []
    cfg = get_config()
    
    if cfg.pipeline.enable_feature_linking:
        if cfg.pipeline.feature_linking_max_workers < 1:
            issues.append("FEATURE_LINKING_MAX_WORKERS must be >= 1")
        
        if cfg.pipeline.feature_linking_max_workers > 100:
            issues.append("FEATURE_LINKING_MAX_WORKERS is very high (>100), may cause performance issues")
    
    return issues


def validate_notion_sync_config() -> List[str]:
    """
    Validate Notion sync configuration and return list of issues.
    
    Returns:
        List of validation error messages (empty if valid)
    """
    issues: List[str] = []
    cfg = get_config()
    
    if cfg.pipeline.enable_notion_sync or cfg.pipeline.enable_dual_write:
        # Check if Notion API key is set
        from amprenta_rag.config import NOTION_API_KEY
        if not NOTION_API_KEY:
            issues.append("NOTION_API_KEY is not set but Notion sync is enabled")
    
    return issues


def validate_all_config() -> dict[str, List[str]]:
    """
    Validate all configuration sections and return issues by category.
    
    Returns:
        Dictionary mapping category -> list of issues
    """
    return {
        "postgres": validate_postgres_config(),
        "feature_linking": validate_feature_linking_config(),
        "notion_sync": validate_notion_sync_config(),
    }


def print_config_summary() -> None:
    """
    Print a summary of current configuration.
    """
    cfg = get_config()
    
    logger.info("=" * 60)
    logger.info("Configuration Summary")
    logger.info("=" * 60)
    
    logger.info(f"Postgres as SoT: {cfg.pipeline.use_postgres_as_sot}")
    logger.info(f"Notion Sync: {cfg.pipeline.enable_notion_sync}")
    logger.info(f"Dual Write: {cfg.pipeline.enable_dual_write}")
    logger.info(f"Feature Linking: {cfg.pipeline.enable_feature_linking}")
    if cfg.pipeline.enable_feature_linking:
        logger.info(f"  Max Workers: {cfg.pipeline.feature_linking_max_workers}")
    
    logger.info("=" * 60)
    
    # Validate and show issues
    all_issues = validate_all_config()
    has_issues = any(all_issues.values())
    
    if has_issues:
        logger.warning("Configuration Issues Found:")
        for category, issues in all_issues.items():
            if issues:
                logger.warning(f"  {category.upper()}:")
                for issue in issues:
                    logger.warning(f"    - {issue}")
    else:
        logger.info("✓ Configuration is valid")
