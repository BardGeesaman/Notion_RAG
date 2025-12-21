from __future__ import annotations

import importlib


def test_automation_init_exports():
    mod = importlib.import_module("amprenta_rag.automation")
    assert set(mod.__all__) == {"fire_trigger", "register_action"}
    assert hasattr(mod, "fire_trigger")
    assert hasattr(mod, "register_action")


def test_ingestion_genomics_init_imports():
    mod = importlib.import_module("amprenta_rag.ingestion.genomics")
    assert hasattr(mod, "__all__")
    assert mod.__all__ == []


def test_metadata_init_exports():
    mod = importlib.import_module("amprenta_rag.metadata")
    assert "classify_and_update_all_literature" in mod.__all__
    assert hasattr(mod, "classify_and_update_all_literature")


def test_migration_init_exports():
    mod = importlib.import_module("amprenta_rag.migration")
    assert "DualWriteManager" in mod.__all__
    assert hasattr(mod, "DualWriteManager")


def test_export_init_exports():
    import pytest
    pytest.importorskip("pptx")
    mod = importlib.import_module("amprenta_rag.export")
    assert set(mod.__all__) == {"generate_experiment_slides", "generate_dataset_slides"}
    assert hasattr(mod, "generate_experiment_slides")
    assert hasattr(mod, "generate_dataset_slides")

