#!/usr/bin/env python3
"""
Configuration validation script.

Validates all configuration settings and provides helpful error messages
for common configuration issues.

Usage:
    python scripts/validate_config.py
"""

import sys
from pathlib import Path

# Add project root to path
project_root = Path(__file__).parent.parent
sys.path.insert(0, str(project_root))

from amprenta_rag.utils.config_validation import (
    print_config_summary,
    validate_all_config,
)
from amprenta_rag.logging_utils import get_logger

logger = get_logger(__name__)


def main():
    """Main validation function."""
    print("\n" + "=" * 60)
    print("Configuration Validation")
    print("=" * 60 + "\n")
    
    try:
        # Print configuration summary
        print_config_summary()
        
        # Validate all config sections
        all_issues = validate_all_config()
        has_issues = any(all_issues.values())
        
        print("\n" + "=" * 60)
        if has_issues:
            print("❌ Configuration validation found issues")
            print("=" * 60 + "\n")
            
            for category, issues in all_issues.items():
                if issues:
                    print(f"\n{category.upper()} Issues:")
                    for issue in issues:
                        print(f"  ⚠️  {issue}")
            
            print("\n" + "=" * 60)
            print("Please fix the issues above and try again.")
            print("=" * 60 + "\n")
            sys.exit(1)
        else:
            print("✅ Configuration is valid!")
            print("=" * 60 + "\n")
            sys.exit(0)
            
    except Exception as e:
        logger.error("Configuration validation failed: %r", e)
        import traceback
        traceback.print_exc()
        sys.exit(1)


if __name__ == "__main__":
    main()

