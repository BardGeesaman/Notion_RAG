"""Tests for drift monitoring functionality."""

from __future__ import annotations

import numpy as np

from amprenta_rag.ml.monitoring.drift import (
    compute_fingerprint_aggregates,
    compute_psi,
    get_drift_status,
)


def test_compute_psi_identical_distributions() -> None:
    """Test PSI computation with identical distributions."""
    reference_bins = [0.25, 0.25, 0.25, 0.25]
    current_bins = [0.25, 0.25, 0.25, 0.25]
    
    psi = compute_psi(reference_bins, current_bins)
    
    # PSI should be very close to 0 for identical distributions
    assert abs(psi) < 1e-6, f"Expected PSI ≈ 0, got {psi}"


def test_compute_psi_different_distributions() -> None:
    """Test PSI computation with different distributions."""
    reference_bins = [0.4, 0.3, 0.2, 0.1]
    current_bins = [0.1, 0.2, 0.3, 0.4]  # Reversed distribution
    
    psi = compute_psi(reference_bins, current_bins)
    
    # PSI should be positive and significant for different distributions
    assert psi > 0.1, f"Expected PSI > 0.1 for different distributions, got {psi}"
    assert psi < 5.0, f"Expected reasonable PSI value, got {psi}"


def test_compute_psi_with_zeros() -> None:
    """Test PSI computation handles zero bins correctly."""
    reference_bins = [0.0, 0.5, 0.5, 0.0]
    current_bins = [0.1, 0.4, 0.4, 0.1]
    
    psi = compute_psi(reference_bins, current_bins)
    
    # Should handle zeros without crashing and return reasonable value
    assert isinstance(psi, float)
    assert not np.isnan(psi)
    assert not np.isinf(psi)


def test_compute_psi_mismatched_lengths() -> None:
    """Test PSI computation raises error for mismatched bin lengths."""
    reference_bins = [0.25, 0.25, 0.25, 0.25]
    current_bins = [0.33, 0.33, 0.33]  # Different length
    
    try:
        compute_psi(reference_bins, current_bins)
        assert False, "Should have raised ValueError for mismatched lengths"
    except ValueError as e:
        assert "same length" in str(e)


def test_get_drift_status_thresholds() -> None:
    """Test drift status determination based on PSI thresholds."""
    # Test good status (PSI < 0.1)
    assert get_drift_status(0.05) == "good"
    assert get_drift_status(0.09) == "good"
    
    # Test warning status (0.1 <= PSI < 0.25)
    assert get_drift_status(0.1) == "warning"
    assert get_drift_status(0.15) == "warning"
    assert get_drift_status(0.24) == "warning"
    
    # Test alert status (PSI >= 0.25)
    assert get_drift_status(0.25) == "alert"
    assert get_drift_status(0.5) == "alert"
    assert get_drift_status(1.0) == "alert"


def test_get_drift_status_edge_cases() -> None:
    """Test drift status with edge case values."""
    # Test boundary values
    assert get_drift_status(0.0) == "good"
    assert get_drift_status(0.099999) == "good"
    assert get_drift_status(0.100001) == "warning"
    assert get_drift_status(0.249999) == "warning"
    assert get_drift_status(0.250001) == "alert"


def test_compute_fingerprint_aggregates() -> None:
    """Test fingerprint aggregate computation."""
    # Create test fingerprint data
    np.random.seed(42)
    fp_array = np.random.choice([True, False], size=(10, 20), p=[0.3, 0.7])
    centroid = np.random.choice([True, False], size=20, p=[0.3, 0.7])
    
    aggregates = compute_fingerprint_aggregates(fp_array, centroid)
    
    # Check that all expected keys are present
    expected_keys = ["pct_bits_set", "tanimoto_to_centroid", "hamming_distance_mean", "hamming_distance_std"]
    for key in expected_keys:
        assert key in aggregates, f"Missing key: {key}"
    
    # Check that values are reasonable
    assert 0 <= aggregates["pct_bits_set"] <= 1, "pct_bits_set should be between 0 and 1"
    assert 0 <= aggregates["tanimoto_to_centroid"] <= 1, "tanimoto should be between 0 and 1"
    assert aggregates["hamming_distance_mean"] >= 0, "hamming distance should be non-negative"
    assert aggregates["hamming_distance_std"] >= 0, "hamming std should be non-negative"


def test_compute_fingerprint_aggregates_empty() -> None:
    """Test fingerprint aggregates with empty input."""
    empty_fp = np.array([]).reshape(0, 10)
    centroid = np.ones(10, dtype=bool)
    
    aggregates = compute_fingerprint_aggregates(empty_fp, centroid)
    
    # Should return zero values for empty input
    assert aggregates["pct_bits_set"] == 0.0
    assert aggregates["tanimoto_to_centroid"] == 0.0
    assert aggregates["hamming_distance_mean"] == 0.0
    assert aggregates["hamming_distance_std"] == 0.0


def test_compute_fingerprint_aggregates_single() -> None:
    """Test fingerprint aggregates with single fingerprint."""
    fp_array = np.array([[True, False, True, False, True]], dtype=bool)
    centroid = np.array([False, True, True, False, False], dtype=bool)
    
    aggregates = compute_fingerprint_aggregates(fp_array, centroid)
    
    # Should compute reasonable values for single fingerprint
    assert 0 <= aggregates["pct_bits_set"] <= 1
    assert 0 <= aggregates["tanimoto_to_centroid"] <= 1
    assert aggregates["hamming_distance_std"] == 0.0  # Only one sample, so std = 0


def test_compute_fingerprint_aggregates_identical() -> None:
    """Test fingerprint aggregates when all fingerprints are identical to centroid."""
    centroid = np.array([True, False, True, False], dtype=bool)
    fp_array = np.array([centroid, centroid, centroid], dtype=bool)
    
    aggregates = compute_fingerprint_aggregates(fp_array, centroid)
    
    # Should have perfect similarity
    assert aggregates["tanimoto_to_centroid"] == 1.0
    assert aggregates["hamming_distance_mean"] == 0.0
    assert aggregates["hamming_distance_std"] == 0.0
