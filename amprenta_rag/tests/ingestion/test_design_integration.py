from __future__ import annotations

from amprenta_rag.ingestion import design_integration as di


def test_module_imports():
    # Module currently only contains thin helpers; smoke-test import
    assert hasattr(di, "__doc__")

