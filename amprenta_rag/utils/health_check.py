"""
System health check utilities.

Provides comprehensive health checks for all system components,
useful for monitoring and diagnostics.
"""

from __future__ import annotations

import time
from dataclasses import dataclass
from typing import Any, Dict, List, Optional

from amprenta_rag.config import get_config
from amprenta_rag.logging_utils import get_logger

logger = get_logger(__name__)


@dataclass
class HealthCheckResult:
    """Result of a health check."""
    component: str
    status: str  # "healthy", "degraded", "unhealthy"
    message: str
    latency_ms: Optional[float] = None
    details: Optional[Dict[str, Any]] = None


class HealthChecker:
    """Perform health checks on system components."""

    def __init__(self):
        self.results: List[HealthCheckResult] = []

    def check_all(self) -> Dict[str, Any]:
        """Run all health checks."""
        logger.info("[HEALTH-CHECK] Running comprehensive health checks...")

        self.results = []

        # Check external services
        self._check_notion()
        # Pinecone deprecated - using pgvector instead
        self._check_openai()

        # Check internal systems
        self._check_configuration()

        # Aggregate results
        healthy = sum(1 for r in self.results if r.status == "healthy")
        degraded = sum(1 for r in self.results if r.status == "degraded")
        unhealthy = sum(1 for r in self.results if r.status == "unhealthy")

        overall_status = "healthy"
        if unhealthy > 0:
            overall_status = "unhealthy"
        elif degraded > 0:
            overall_status = "degraded"

        return {
            "overall_status": overall_status,
            "components": {
                "healthy": healthy,
                "degraded": degraded,
                "unhealthy": unhealthy,
                "total": len(self.results),
            },
            "checks": [self._result_to_dict(r) for r in self.results],
            "timestamp": time.time(),
        }

    def _check_notion(self):
        """Check Notion API health."""
        start = time.time()
        try:
            import requests

            cfg = get_config()
            url = f"{cfg.notion.base_url}/users/me"
            headers = {
                "Authorization": f"Bearer {cfg.notion.api_key}",
                "Notion-Version": cfg.notion.version,
            }

            resp = requests.get(url, headers=headers, timeout=5)
            latency = (time.time() - start) * 1000

            if resp.status_code == 200:
                user_data = resp.json()
                self.results.append(
                    HealthCheckResult(
                        component="Notion API",
                        status="healthy",
                        message="Connected successfully",
                        latency_ms=latency,
                        details={"user_id": user_data.get("id", "unknown")},
                    )
                )
            elif resp.status_code == 401:
                self.results.append(
                    HealthCheckResult(
                        component="Notion API",
                        status="unhealthy",
                        message="Authentication failed (invalid API key)",
                        latency_ms=latency,
                    )
                )
            else:
                self.results.append(
                    HealthCheckResult(
                        component="Notion API",
                        status="degraded",
                        message=f"HTTP {resp.status_code}",
                        latency_ms=latency,
                    )
                )
        except requests.exceptions.Timeout:
            latency = (time.time() - start) * 1000
            self.results.append(
                HealthCheckResult(
                    component="Notion API",
                    status="unhealthy",
                    message="Connection timeout",
                    latency_ms=latency,
                )
            )
        except Exception as e:
            latency = (time.time() - start) * 1000
            self.results.append(
                HealthCheckResult(
                    component="Notion API",
                    status="unhealthy",
                    message=f"Error: {str(e)[:100]}",
                    latency_ms=latency,
                )
            )

    # Pinecone support removed - using pgvector backend
    pass

    def _check_openai(self):
        """Check OpenAI API health."""
        start = time.time()
        try:
            import openai

            cfg = get_config()
            client = openai.OpenAI(api_key=cfg.openai.api_key)

            # Simple test call (list models is lightweight)
            models = client.models.list()
            latency = (time.time() - start) * 1000

            model_count = len(list(models))
            self.results.append(
                HealthCheckResult(
                    component="OpenAI API",
                    status="healthy",
                    message="Connected successfully",
                    latency_ms=latency,
                    details={"available_models": model_count},
                )
            )
        except openai.AuthenticationError:
            latency = (time.time() - start) * 1000
            self.results.append(
                HealthCheckResult(
                    component="OpenAI API",
                    status="unhealthy",
                    message="Authentication failed (invalid API key)",
                    latency_ms=latency,
                )
            )
        except Exception as e:
            latency = (time.time() - start) * 1000
            self.results.append(
                HealthCheckResult(
                    component="OpenAI API",
                    status="degraded",
                    message=f"Error: {str(e)[:100]}",
                    latency_ms=latency,
                )
            )

    def _check_configuration(self):
        """Check configuration health."""
        try:
            cfg = get_config()

            # Check critical config values
            issues = []

            if not cfg.openai.api_key:
                issues.append("OpenAI API key missing")
            # Pinecone deprecated - using pgvector
            if not cfg.notion.api_key:
                issues.append("Notion API key missing")

            if issues:
                self.results.append(
                    HealthCheckResult(
                        component="Configuration",
                        status="unhealthy",
                        message="; ".join(issues),
                    )
                )
            else:
                self.results.append(
                    HealthCheckResult(
                        component="Configuration",
                        status="healthy",
                        message="All critical configuration present",
                    )
                )
        except Exception as e:
            self.results.append(
                HealthCheckResult(
                    component="Configuration",
                    status="unhealthy",
                    message=f"Configuration error: {str(e)[:100]}",
                )
            )

    def _result_to_dict(self, result: HealthCheckResult) -> Dict[str, Any]:
        """Convert HealthCheckResult to dictionary."""
        return {
            "component": result.component,
            "status": result.status,
            "message": result.message,
            "latency_ms": result.latency_ms,
            "details": result.details,
        }

    def get_summary(self) -> str:
        """Get a formatted health check summary."""
        lines = ["System Health Check Summary", "=" * 60, ""]

        for result in self.results:
            status_icon = {
                "healthy": "✅",
                "degraded": "⚠️",
                "unhealthy": "❌",
            }.get(result.status, "❓")

            lines.append(f"{status_icon} {result.component}: {result.status.upper()}")
            lines.append(f"   {result.message}")
            if result.latency_ms:
                lines.append(f"   Latency: {result.latency_ms:.1f}ms")
            if result.details:
                for key, value in result.details.items():
                    lines.append(f"   {key}: {value}")
            lines.append("")

        healthy = sum(1 for r in self.results if r.status == "healthy")
        total = len(self.results)

        lines.append(f"Overall: {healthy}/{total} components healthy")

        return "\n".join(lines)


def check_system_health() -> Dict[str, Any]:
    """
    Perform comprehensive system health check.

    Returns:
        Dictionary with health check results
    """
    checker = HealthChecker()
    results = checker.check_all()

    logger.info("\n%s", checker.get_summary())

    return results

