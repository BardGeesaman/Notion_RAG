#!/usr/bin/env python3
"""
System health check script.

Checks the health of all system components and external services.
Useful for monitoring and diagnostics.
"""

import sys
from pathlib import Path

# Add parent directory to path
sys.path.insert(0, str(Path(__file__).parent.parent))

from amprenta_rag.logging_utils import get_logger
from amprenta_rag.utils.health_check import check_system_health

logger = get_logger(__name__)


def main():
    """Run system health checks."""
    print("\n" + "=" * 60)
    print("System Health Check")
    print("=" * 60 + "\n")

    try:
        results = check_system_health()

        overall_status = results["overall_status"]

        if overall_status == "healthy":
            print("\n✅ System is healthy!")
            sys.exit(0)
        elif overall_status == "degraded":
            print("\n⚠️  System is degraded (some components have issues).")
            sys.exit(1)
        else:
            print("\n❌ System is unhealthy (critical components failing).")
            sys.exit(2)

    except Exception as e:
        logger.error("Health check failed: %r", e)
        print(f"\n❌ Health check failed: {e}\n")
        sys.exit(3)


if __name__ == "__main__":
    main()

