"""
Pytest configuration and shared fixtures.

This file is automatically discovered by pytest and provides:
- Shared fixtures for all tests
- Test configuration
- Mock setup helpers
"""

import pytest
import os
from pathlib import Path
from unittest.mock import Mock, patch
from typing import Generator

# Set test environment variables before imports
os.environ.setdefault("USE_POSTGRES_AS_SOT", "true")
os.environ.setdefault("ENABLE_NOTION_SYNC", "false")
os.environ.setdefault("ENABLE_DUAL_WRITE", "false")


@pytest.fixture(scope="session")
def test_data_dir() -> Path:
    """Return path to test data directory."""
    return Path(__file__).parent / "fixtures" / "sample_data"


@pytest.fixture
def sample_metabolomics_data(tmp_path) -> Path:
    """Create sample metabolomics CSV file."""
    file_path = tmp_path / "metabolomics_sample.csv"
    file_path.write_text("""
Metabolite,Value,StdDev
Glucose,1.5,0.1
Lactate,2.3,0.2
Pyruvate,0.8,0.05
Citrate,0.4,0.02
""".strip())
    return file_path


@pytest.fixture
def sample_proteomics_data(tmp_path) -> Path:
    """Create sample proteomics CSV file."""
    file_path = tmp_path / "proteomics_sample.csv"
    file_path.write_text("""
Protein,Log2FC,PValue
P12345,1.2,0.001
Q67890,-0.8,0.05
R11111,0.5,0.01
""".strip())
    return file_path


@pytest.fixture
def sample_transcriptomics_data(tmp_path) -> Path:
    """Create sample transcriptomics CSV file."""
    file_path = tmp_path / "transcriptomics_sample.csv"
    file_path.write_text("""
Gene,Log2FC,PValue,AdjPValue
GENE1,2.5,0.001,0.005
GENE2,-1.8,0.01,0.02
GENE3,0.9,0.05,0.08
""".strip())
    return file_path


@pytest.fixture
def sample_lipidomics_data(tmp_path) -> Path:
    """Create sample lipidomics CSV file."""
    file_path = tmp_path / "lipidomics_sample.csv"
    file_path.write_text("""
Lipid,Concentration,Units
Cer(d18:1/16:0),1.2,ng/mL
Cer(d18:1/18:0),0.8,ng/mL
SM(d18:1/16:0),2.5,ng/mL
""".strip())
    return file_path


@pytest.fixture
def mock_pinecone_client():
    """Mock Pinecone client."""
    with patch("amprenta_rag.clients.pinecone_client") as mock:
        yield mock


@pytest.fixture
def mock_postgres_session():
    """Mock Postgres database session."""
    with patch("amprenta_rag.database.base.get_db") as mock_get_db:
        mock_session = Mock()
        mock_get_db.return_value = iter([mock_session])
        yield mock_session


@pytest.fixture
def mock_config_postgres_only():
    """Mock config with Postgres-only mode."""
    with patch("amprenta_rag.config.get_config") as mock_get_config:
        mock_config = Mock()
        mock_config.pipeline.use_postgres_as_sot = True
        mock_config.pipeline.enable_notion_sync = False
        mock_config.pipeline.enable_dual_write = False
        mock_config.pipeline.enable_feature_linking = False
        mock_get_config.return_value = mock_config
        yield mock_config


@pytest.fixture
def mock_config_notion_sync():
    """Deprecated: legacy sync is disabled in tests; kept only for compatibility."""
    yield None


# Markers for test organization (Notion-related marker removed)
def pytest_configure(config):
    """Configure custom pytest markers."""
    config.addinivalue_line(
        "markers", "unit: Unit tests (fast, isolated)"
    )
    config.addinivalue_line(
        "markers", "integration: Integration tests (require external services)"
    )
    config.addinivalue_line(
        "markers", "slow: Slow tests that take >1 second"
    )
    config.addinivalue_line(
        "markers", "requires_postgres: Tests requiring Postgres database"
    )
    config.addinivalue_line(
        "markers", "requires_pinecone: Tests requiring Pinecone API"
    )

