"""
Tests for Flow Cytometry FCS parsing and transformation pipeline.

Tests cover FCS file parsing, metadata extraction, data transformations,
Parquet I/O, and validation functions.
"""

import tempfile
import pytest
import numpy as np
import pandas as pd
from pathlib import Path
from unittest.mock import patch, MagicMock
from datetime import datetime

from amprenta_rag.flow_cytometry.fcs_parser import (
    load_fcs,
    extract_metadata,
    validate_fcs,
    save_events_parquet,
    load_events_parquet,
    FlowMetadata,
)

from amprenta_rag.flow_cytometry.transforms import (
    logicle_transform,
    arcsinh_transform,
    apply_compensation,
    auto_logicle_params,
    subsample_events,
)


@pytest.fixture
def mock_fcs_metadata():
    """Create mock FCS metadata dictionary."""
    return {
        '$FIL': 'FCS3.1',
        '$TOT': '10000',
        '$PAR': '4',
        '$P1N': 'FSC-A',
        '$P1S': 'FSC-A',
        '$P1R': '262144',
        '$P2N': 'SSC-A',
        '$P2S': 'SSC-A', 
        '$P2R': '262144',
        '$P3N': 'CD3-FITC',
        '$P3S': 'CD3',
        '$P3R': '262144',
        '$P4N': 'CD4-PE',
        '$P4S': 'CD4',
        '$P4R': '262144',
        '$DATE': '01-JAN-2024',
        '$BTIM': '10:30:00',
        '$CYT': 'FACSCanto II',
        '$CYTSN': 'A12345',
        '$SRC': 'FACSDiva',
        '$SMNO': 'Sample001',
        '$VOL': '100.0',
        '$SPILLOVER': '4,FSC-A,SSC-A,CD3-FITC,CD4-PE,1.0,0.1,0.05,0.02,0.1,1.0,0.15,0.03,0.05,0.15,1.0,0.2,0.02,0.03,0.2,1.0',
    }


@pytest.fixture
def mock_fcs_events():
    """Create mock FCS event data."""
    np.random.seed(42)
    n_events = 1000
    n_params = 4
    
    # Generate realistic flow cytometry data
    events = np.random.lognormal(mean=8, sigma=1, size=(n_events, n_params)).astype(np.float32)
    
    # Add some negative values for testing
    events[:50, 2] = -np.random.exponential(100, 50)
    
    return events


def test_flow_metadata_creation():
    """Test FlowMetadata dataclass creation and defaults."""
    metadata = FlowMetadata(
        filename="test.fcs",
        file_size_bytes=1024,
        fcs_version="FCS3.1",
        n_events=1000,
        n_parameters=4
    )
    
    assert metadata.filename == "test.fcs"
    assert metadata.file_size_bytes == 1024
    assert metadata.fcs_version == "FCS3.1"
    assert metadata.n_events == 1000
    assert metadata.n_parameters == 4
    assert metadata.parameter_names == []
    assert metadata.parameter_ranges == {}
    assert metadata.keywords == {}


@patch('amprenta_rag.flow_cytometry.fcs_parser.parse')
def test_load_fcs_success(mock_parse, mock_fcs_metadata, mock_fcs_events):
    """Test successful FCS file loading."""
    mock_parse.return_value = (mock_fcs_metadata, pd.DataFrame(mock_fcs_events))
    
    with tempfile.NamedTemporaryFile(suffix='.fcs', delete=False) as tmp:
        tmp_path = Path(tmp.name)
    
    try:
        meta, events = load_fcs(tmp_path)
        
        assert meta == mock_fcs_metadata
        assert events.shape == mock_fcs_events.shape
        assert events.dtype == np.float32
        mock_parse.assert_called_once_with(str(tmp_path), reformat_meta=True)
        
    finally:
        tmp_path.unlink()


def test_load_fcs_file_not_found():
    """Test FCS loading with non-existent file."""
    with pytest.raises(FileNotFoundError, match="FCS file not found"):
        load_fcs("nonexistent.fcs")


@patch('amprenta_rag.flow_cytometry.fcs_parser.parse')
def test_load_fcs_parse_error(mock_parse):
    """Test FCS loading with parsing error."""
    mock_parse.side_effect = Exception("Parse failed")
    
    with tempfile.NamedTemporaryFile(suffix='.fcs', delete=False) as tmp:
        tmp_path = Path(tmp.name)
    
    try:
        with pytest.raises(RuntimeError, match="Failed to parse FCS file"):
            load_fcs(tmp_path)
    finally:
        tmp_path.unlink()


def test_extract_metadata_complete(mock_fcs_metadata):
    """Test metadata extraction with complete FCS metadata."""
    metadata = extract_metadata(mock_fcs_metadata, "test.fcs", 1024)
    
    assert metadata.filename == "test.fcs"
    assert metadata.file_size_bytes == 1024
    assert metadata.fcs_version == "FCS3.1"
    assert metadata.n_events == 10000
    assert metadata.n_parameters == 4
    assert metadata.parameter_names == ['FSC-A', 'SSC-A', 'CD3', 'CD4']
    assert metadata.cytometer_model == "FACSCanto II"
    assert metadata.cytometer_serial == "A12345"
    assert metadata.sample_id == "Sample001"
    assert metadata.sample_volume == 100.0
    assert metadata.acquisition_date == datetime(2024, 1, 1, 10, 30, 0)
    assert metadata.spillover_matrix is not None
    assert metadata.spillover_matrix.shape == (4, 4)
    # Check that diagonal elements are 1.0 (self-compensation)
    assert np.allclose(np.diag(metadata.spillover_matrix), 1.0)


def test_extract_metadata_minimal():
    """Test metadata extraction with minimal FCS metadata."""
    minimal_meta = {
        '$FIL': 'FCS2.0',
        '$TOT': '5000',
        '$PAR': '2',
        '$P1N': 'FSC',
        '$P2N': 'SSC',
    }
    
    metadata = extract_metadata(minimal_meta)
    
    assert metadata.fcs_version == "FCS2.0"
    assert metadata.n_events == 5000
    assert metadata.n_parameters == 2
    assert metadata.parameter_names == ['FSC', 'SSC']
    assert metadata.acquisition_date is None
    assert metadata.cytometer_model is None
    assert metadata.spillover_matrix is None


@patch('amprenta_rag.flow_cytometry.fcs_parser.load_fcs')
def test_validate_fcs_success(mock_load_fcs, mock_fcs_metadata, mock_fcs_events):
    """Test FCS validation with valid file."""
    # Make sure metadata matches the mock events
    mock_fcs_metadata['$TOT'] = str(mock_fcs_events.shape[0])  # Match event count
    mock_fcs_metadata['$PAR'] = str(mock_fcs_events.shape[1])  # Match parameter count
    mock_load_fcs.return_value = (mock_fcs_metadata, mock_fcs_events)
    
    with tempfile.NamedTemporaryFile(suffix='.fcs', delete=False) as tmp:
        # Write some content to avoid "empty file" error
        tmp.write(b"dummy content")
        tmp_path = Path(tmp.name)
    
    try:
        issues = validate_fcs(tmp_path)
        assert len(issues) == 0
    finally:
        tmp_path.unlink()


def test_validate_fcs_file_not_found():
    """Test FCS validation with non-existent file."""
    issues = validate_fcs("nonexistent.fcs")
    assert len(issues) == 1
    assert "File not found" in issues[0]


@patch('amprenta_rag.flow_cytometry.fcs_parser.load_fcs')
def test_validate_fcs_dimension_mismatch(mock_load_fcs, mock_fcs_metadata):
    """Test FCS validation with dimension mismatch."""
    # Create mismatched data
    wrong_events = np.random.random((5000, 3))  # Wrong event count and params
    mock_load_fcs.return_value = (mock_fcs_metadata, wrong_events)
    
    with tempfile.NamedTemporaryFile(suffix='.fcs', delete=False) as tmp:
        # Write some content to avoid "empty file" error
        tmp.write(b"dummy content")
        tmp_path = Path(tmp.name)
    
    try:
        issues = validate_fcs(tmp_path)
        assert len(issues) >= 2  # Both event and parameter count mismatches
        assert any("Event count mismatch" in issue for issue in issues)
        assert any("Parameter count mismatch" in issue for issue in issues)
    finally:
        tmp_path.unlink()


def test_parquet_roundtrip(mock_fcs_events):
    """Test Parquet save and load roundtrip."""
    param_names = ['FSC-A', 'SSC-A', 'CD3-FITC', 'CD4-PE']
    metadata = {'test_key': 'test_value'}
    
    with tempfile.NamedTemporaryFile(suffix='.parquet', delete=False) as tmp:
        tmp_path = Path(tmp.name)
    
    try:
        # Save to Parquet
        save_events_parquet(mock_fcs_events, param_names, tmp_path, metadata)
        assert tmp_path.exists()
        
        # Load from Parquet
        loaded_df = load_events_parquet(tmp_path)
        
        # Verify data integrity
        assert loaded_df.shape == mock_fcs_events.shape
        assert list(loaded_df.columns) == param_names
        np.testing.assert_array_almost_equal(loaded_df.values, mock_fcs_events)
        
    finally:
        tmp_path.unlink()


def test_parquet_dimension_mismatch():
    """Test Parquet save with mismatched dimensions."""
    events = np.random.random((100, 3))
    param_names = ['A', 'B']  # Wrong number of names
    
    with tempfile.NamedTemporaryFile(suffix='.parquet', delete=False) as tmp:
        tmp_path = Path(tmp.name)
    
    try:
        with pytest.raises(ValueError, match="Parameter count mismatch"):
            save_events_parquet(events, param_names, tmp_path)
    finally:
        if tmp_path.exists():
            tmp_path.unlink()


def test_logicle_transform_basic():
    """Test basic logicle transformation."""
    data = np.array([100, 1000, 10000], dtype=np.float32)
    T, M, W, A = 10000, 4.5, 0.5, 0.0
    
    result = logicle_transform(data, T, M, W, A)
    
    # Basic functionality tests
    assert result.shape == data.shape
    assert np.all(np.isfinite(result))
    # Transformation should change the values
    assert not np.allclose(result, data)
    # Should return numeric results
    assert result.dtype in [np.float32, np.float64]


def test_logicle_transform_parameter_validation():
    """Test logicle parameter validation."""
    data = np.array([100, 1000])
    
    # Test invalid T
    with pytest.raises(ValueError, match="T must be positive"):
        logicle_transform(data, T=-100, M=4, W=0.5, A=0)
    
    # Test invalid M  
    with pytest.raises(ValueError, match="M must be positive"):
        logicle_transform(data, T=1000, M=-1, W=0.5, A=0)
    
    # Test invalid W
    with pytest.raises(ValueError, match="W must be in range"):
        logicle_transform(data, T=1000, M=4, W=3, A=0)  # W >= M/2
    
    # Test invalid A
    with pytest.raises(ValueError, match="A must be non-negative"):
        logicle_transform(data, T=1000, M=4, W=0.5, A=-1)


def test_arcsinh_transform():
    """Test arcsinh transformation."""
    data = np.array([0, 150, 300, 1500], dtype=np.float32)
    cofactor = 150
    
    result = arcsinh_transform(data, cofactor)
    expected = np.arcsinh(data / cofactor)
    
    np.testing.assert_array_almost_equal(result, expected)


def test_arcsinh_transform_invalid_cofactor():
    """Test arcsinh with invalid cofactor."""
    data = np.array([100, 200])
    
    with pytest.raises(ValueError, match="Cofactor must be positive"):
        arcsinh_transform(data, cofactor=-150)


def test_apply_compensation():
    """Test spectral compensation application."""
    # Create test data and spillover matrix
    data = np.array([
        [1000, 200, 50, 10],
        [100, 1500, 100, 20],
        [50, 100, 2000, 150]
    ], dtype=np.float32)
    
    spillover = np.array([
        [1.0, 0.1, 0.05, 0.02],
        [0.1, 1.0, 0.15, 0.03],
        [0.05, 0.15, 1.0, 0.2],
        [0.02, 0.03, 0.2, 1.0]
    ])
    
    compensated = apply_compensation(data, spillover)
    
    assert compensated.shape == data.shape
    assert np.all(np.isfinite(compensated))
    # Compensation matrix inversion should produce reasonable values
    # Just check that compensation was applied (values changed)
    assert not np.allclose(compensated, data)


def test_apply_compensation_dimension_mismatch():
    """Test compensation with mismatched dimensions."""
    data = np.array([[100, 200, 300]], dtype=np.float32)
    spillover = np.array([[1.0, 0.1], [0.1, 1.0]])  # Wrong size
    
    with pytest.raises(ValueError, match="Data channels.*must match spillover matrix"):
        apply_compensation(data, spillover)


def test_apply_compensation_non_square_matrix():
    """Test compensation with non-square spillover matrix."""
    data = np.array([[100, 200]], dtype=np.float32)
    spillover = np.array([[1.0, 0.1, 0.05], [0.1, 1.0, 0.15]])  # Not square
    
    with pytest.raises(ValueError, match="Spillover matrix must be square"):
        apply_compensation(data, spillover)


def test_auto_logicle_params():
    """Test automatic logicle parameter estimation."""
    # Create test data with positive and negative values
    np.random.seed(42)
    positive_data = np.random.lognormal(8, 1, 5000)
    negative_data = -np.random.exponential(100, 500)
    data = np.concatenate([positive_data, negative_data])
    
    T, M, W, A = auto_logicle_params(data)
    
    # Validate parameter bounds
    assert T > 0
    assert M > 0
    assert 0 < W < M/2
    assert A >= 0
    
    # Test with positive-only data
    T2, M2, W2, A2 = auto_logicle_params(positive_data)
    assert A2 == 0.0  # No negative data


def test_auto_logicle_params_empty_data():
    """Test auto parameter estimation with empty data."""
    with pytest.raises(ValueError, match="Data array cannot be empty"):
        auto_logicle_params(np.array([]))


def test_subsample_events():
    """Test event subsampling."""
    np.random.seed(42)
    events = np.random.random((10000, 4))
    max_events = 5000
    
    subsampled = subsample_events(events, max_events, random_state=42)
    
    assert subsampled.shape == (max_events, 4)
    assert subsampled.shape[1] == events.shape[1]
    
    # Should be a subset of original data
    for i in range(subsampled.shape[0]):
        found = False
        for j in range(events.shape[0]):
            if np.allclose(subsampled[i], events[j]):
                found = True
                break
        assert found, f"Subsampled row {i} not found in original data"


def test_subsample_events_no_subsampling_needed():
    """Test subsampling when no subsampling is needed."""
    events = np.random.random((1000, 4))
    max_events = 5000
    
    result = subsample_events(events, max_events)
    
    np.testing.assert_array_equal(result, events)


def test_subsample_events_invalid_max_events():
    """Test subsampling with invalid max_events."""
    events = np.random.random((100, 4))
    
    with pytest.raises(ValueError, match="max_events must be positive"):
        subsample_events(events, max_events=-100)
