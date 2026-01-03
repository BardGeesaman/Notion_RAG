"""Test configuration defaults.

Import this module in conftest.py to set up test environment.
These values are ONLY for testing - never used in production.
"""

import os

# Test-only secrets (not real values)
TEST_DEFAULTS = {
    "SIGNATURE_SECRET_KEY": "test-only-signature-key-do-not-use-in-prod",
    "JWT_SECRET_KEY": "test-only-jwt-key-do-not-use-in-prod",
    "NOTEBOOK_REVIEW_SECRET": "test-only-review-key-do-not-use-in-prod",
}


def configure_test_environment():
    """Set test defaults for missing env vars."""
    for key, value in TEST_DEFAULTS.items():
        if key not in os.environ:
            os.environ[key] = value
