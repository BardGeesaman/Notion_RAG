"""Tests for export engine."""

from __future__ import annotations

import json
from unittest.mock import MagicMock
from uuid import uuid4

import pytest

from amprenta_rag.export.export_engine import (
    create_export_package,
    export_compounds,
    export_dataset,
    export_experiment,
    generate_manifest,
)


class TestExportDataset:
    """Tests for export_dataset."""

    def test_export_dataset_csv(self):
        """Test dataset export to CSV."""
        mock_feat = MagicMock()
        mock_feat.name = "gene1"
        mock_feat.external_ids = {"value": 2.5}
        
        mock_dataset = MagicMock()
        mock_dataset.id = uuid4()
        mock_dataset.name = "Test Dataset"
        mock_dataset.omics_type = "transcriptomics"
        mock_dataset.features = [mock_feat]
        
        mock_db = MagicMock()
        mock_db.query.return_value.filter.return_value.first.return_value = mock_dataset
        
        result = export_dataset(mock_db, str(mock_dataset.id), format="csv")
        
        assert b"feature_name" in result
        assert b"gene1" in result

    def test_export_dataset_not_found(self):
        """Test export with non-existent dataset."""
        mock_db = MagicMock()
        mock_db.query.return_value.filter.return_value.first.return_value = None
        
        with pytest.raises(ValueError, match="not found"):
            export_dataset(mock_db, str(uuid4()), format="csv")


class TestExportExperiment:
    """Tests for export_experiment."""

    def test_export_experiment_json(self):
        """Test experiment export to JSON."""
        mock_exp = MagicMock()
        mock_exp.id = uuid4()
        mock_exp.name = "Test Experiment"
        mock_exp.type = "RNA-seq"
        mock_exp.datasets = []
        
        mock_db = MagicMock()
        mock_db.query.return_value.filter.return_value.first.return_value = mock_exp
        
        result = export_experiment(mock_db, str(mock_exp.id), format="json")
        
        data = json.loads(result)
        assert data["name"] == "Test Experiment"


class TestExportCompounds:
    """Tests for export_compounds."""

    def test_export_compounds_csv(self):
        """Test compounds export to CSV."""
        mock_comp = MagicMock()
        mock_comp.id = uuid4()
        mock_comp.compound_id = "CMPD001"
        mock_comp.canonical_smiles = "CCO"
        mock_comp.smiles = "CCO"
        mock_comp.name = "Ethanol"
        
        mock_db = MagicMock()
        mock_db.query.return_value.filter.return_value.all.return_value = [mock_comp]
        
        result = export_compounds(mock_db, [str(mock_comp.id)], format="csv")
        
        assert b"compound_id" in result
        assert b"CMPD001" in result


class TestGenerateManifest:
    """Tests for generate_manifest."""

    def test_generate_manifest(self):
        """Test manifest generation."""
        items = ["dataset1.csv", "dataset2.csv"]
        checksums = {
            "dataset1.csv": "abc123",
            "dataset2.csv": "def456",
        }
        
        manifest = generate_manifest(items, checksums)
        
        assert "export_timestamp" in manifest
        assert "version" in manifest
        assert manifest["items"] == items
        assert manifest["checksums"] == checksums


class TestCreateExportPackage:
    """Tests for create_export_package."""

    def test_create_zip_package(self):
        """Test ZIP package creation."""
        exports = {
            "dataset1.csv": b"name,value\ngene1,2.5",
            "dataset2.csv": b"name,value\ngene2,1.8",
        }
        
        result = create_export_package(exports, format="zip")
        
        # Should be valid ZIP
        assert result.startswith(b"PK")  # ZIP magic number
        assert len(result) > 0

    def test_create_package_unsupported_format(self):
        """Test package creation with unsupported format."""
        exports = {"test.csv": b"data"}
        
        with pytest.raises(ValueError, match="Unsupported package format"):
            create_export_package(exports, format="tar")

