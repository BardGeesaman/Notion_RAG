"""Shared error utilities for CLI and application messaging."""
from __future__ import annotations

from typing import Dict, List, TypedDict

from amprenta_rag.logging_utils import get_logger

logger = get_logger(__name__)

class ErrorEntry(TypedDict):
    message: str
    suggestions: List[str]


ERROR_MESSAGES: Dict[str, ErrorEntry] = {
    "db_connection": {
        "message": "Database connection failed: {details}",
        "suggestions": [
            "Verify database URL and credentials",
            "Check network access to the database host",
            "Ensure migrations are up to date",
        ],
    },
    "missing_api_key": {
        "message": "API key is missing or invalid for {service}",
        "suggestions": [
            "Set the required API key in environment variables",
            "Double-check the key value and permissions",
            "Restart the service after setting credentials",
        ],
    },
    "file_not_found": {
        "message": "File not found: {path}",
        "suggestions": [
            "Confirm the path exists and is readable",
            "Check working directory and relative paths",
            "Regenerate or download the missing file",
        ],
    },
    "parse_error": {
        "message": "Failed to parse input: {details}",
        "suggestions": [
            "Validate input format and encoding",
            "Check for truncated or malformed content",
            "Retry with a known-good sample",
        ],
    },
    "validation_error": {
        "message": "Validation failed: {details}",
        "suggestions": [
            "Review required fields and data types",
            "Fix invalid values and re-run",
            "Run schema or contract validation tools",
        ],
    },
    "network_error": {
        "message": "Network request failed: {details}",
        "suggestions": [
            "Check internet connectivity and proxies",
            "Retry the request; transient failures are common",
            "Increase timeout or add retries with backoff",
        ],
    },
}


def format_error(category: str, **kwargs: str) -> str:
    """Format an error message with suggestions for the given category."""
    details = kwargs.get("details", "unspecified error")
    service = kwargs.get("service", "service")
    path = kwargs.get("path", "")

    entry = ERROR_MESSAGES.get(category)
    if not entry:
        return f"[{category}] {details}"

    message = entry["message"].format(details=details, service=service, path=path)
    suggestions = entry.get("suggestions", [])
    if not suggestions:
        return message

    tips = "\n".join(f"- {tip}" for tip in suggestions)
    return f"{message}\nSuggestions:\n{tips}"


def render_cli_error(category: str, **kwargs: str) -> None:
    """Print a colored CLI error with suggestions."""
    text = format_error(category, **kwargs)
    red = "\033[91m"
    reset = "\033[0m"
    logger.error(f"{red}{text}{reset}")


__all__ = ["ERROR_MESSAGES", "format_error", "render_cli_error"]
