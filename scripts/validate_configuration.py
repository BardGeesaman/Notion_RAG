#!/usr/bin/env python3
"""
Configuration validation script.

Validates all configuration and system readiness before starting
the application.
"""

import sys
from pathlib import Path

# Add parent directory to path
sys.path.insert(0, str(Path(__file__).parent.parent))

from amprenta_rag.logging_utils import get_logger
from amprenta_rag.utils.config_validation import validate_configuration

logger = get_logger(__name__)


def main():
    """Run configuration validation."""
    print("\n" + "=" * 60)
    print("Configuration Validation")
    print("=" * 60 + "\n")

    try:
        passed = validate_configuration()

        if passed:
            print("\n✅ All configuration checks passed!")
            print("System is ready to use.\n")
            sys.exit(0)
        else:
            print("\n⚠️  Some configuration checks failed.")
            print("Please review the warnings above.\n")
            sys.exit(1)

    except Exception as e:
        logger.error("Configuration validation failed: %r", e)
        print(f"\n❌ Configuration validation failed: {e}\n")
        sys.exit(1)


if __name__ == "__main__":
    main()

