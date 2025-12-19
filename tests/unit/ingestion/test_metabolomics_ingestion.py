"""
Unit tests for metabolomics ingestion.

This module tests the metabolomics ingestion pipeline in isolation,
using mocks for external services (Postgres, Notion, Pinecone).
"""

import pytest
from unittest.mock import Mock, patch
from uuid import uuid4

# Import the function we're testing
from amprenta_rag.ingestion.metabolomics.ingestion import ingest_metabolomics_file


class TestMetabolomicsIngestion:
    """Test suite for metabolomics ingestion."""

    @pytest.fixture
    def sample_metabolomics_file(self, tmp_path):
        """Create a sample metabolomics CSV file for testing."""
        file_path = tmp_path / "metabolomics_sample.csv"
        file_path.write_text("""
Metabolite,Value,StdDev
Glucose,1.5,0.1
Lactate,2.3,0.2
Pyruvate,0.8,0.05
""".strip())
        return str(file_path)

    @pytest.fixture
    def mock_postgres_dataset(self):
        """Mock Postgres dataset object."""
        dataset = Mock()
        dataset.id = uuid4()
        dataset.name = "Test Dataset"
        dataset.omics_type = "METABOLOMICS"
        return dataset

    def test_ingest_metabolomics_file_missing_file(self):
        """Test that ingestion fails gracefully when file doesn't exist."""
        with pytest.raises(FileNotFoundError, match="Metabolomics file not found"):
            ingest_metabolomics_file("nonexistent_file.csv")

    @patch("amprenta_rag.ingestion.metabolomics.ingestion.get_config")
    @patch("amprenta_rag.ingestion.metabolomics.ingestion.create_or_update_dataset_in_postgres")
    def test_ingest_metabolomics_file_postgres_only(
        self,
        mock_create_dataset,
        mock_get_config,
        sample_metabolomics_file,
        mock_postgres_dataset,
    ):
        """Test metabolomics ingestion with Postgres-only mode."""
        # Setup mocks
        mock_config = Mock()
        mock_config.pipeline.use_postgres_as_sot = True
        mock_config.pipeline.enable_notion_sync = False
        mock_config.pipeline.enable_dual_write = False
        mock_config.pipeline.enable_feature_linking = False
        mock_get_config.return_value = mock_config

        mock_create_dataset.return_value = mock_postgres_dataset

        # Mock embedding
        with patch("amprenta_rag.ingestion.metabolomics.ingestion.embed_dataset_with_postgres_metadata"):
            # Run ingestion
            result = ingest_metabolomics_file(
                file_path=sample_metabolomics_file,
                create_page=False,
            )

            # Verify Postgres dataset was created
            mock_create_dataset.assert_called_once()
            assert result is not None

    @patch("amprenta_rag.ingestion.metabolomics.ingestion.get_config")
    def test_ingest_metabolomics_file_empty_file(self, mock_get_config, tmp_path):
        """Test that ingestion fails when file has no metabolites."""
        # Create empty file
        empty_file = tmp_path / "empty.csv"
        empty_file.write_text("Metabolite,Value\n")

        mock_config = Mock()
        mock_config.pipeline.use_postgres_as_sot = True
        mock_get_config.return_value = mock_config

        with pytest.raises(ValueError, match="No metabolites extracted"):
            ingest_metabolomics_file(str(empty_file))

    # TODO: Add more tests:
    # - Test with Notion sync enabled
    # - Test feature linking
    # - Test program/experiment linking
    # - Test error handling for Postgres failures
    # - Test error handling for Notion failures

