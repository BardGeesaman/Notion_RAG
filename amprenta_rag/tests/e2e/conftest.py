"""Re-export dashboard E2E fixtures for tests/e2e/.

We keep the underlying server-launch logic in `tests/dashboard/conftest.py`.
"""

from __future__ import annotations

# Re-export fixtures so they are discoverable under this test package.
from amprenta_rag.tests.dashboard.conftest import (  # noqa: F401
    fastapi_server,
    streamlit_server,
    browser_type_launch_args,
)


