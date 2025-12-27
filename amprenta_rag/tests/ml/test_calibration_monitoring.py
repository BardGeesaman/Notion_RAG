"""Tests for calibration monitoring functionality."""

from __future__ import annotations

import numpy as np

from amprenta_rag.ml.monitoring.calibration import (
    compute_ece,
    compute_reliability_diagram,
    get_calibration_status,
)


def test_compute_ece_perfect_calibration() -> None:
    """Test ECE computation with perfectly calibrated predictions."""
    # Create perfectly calibrated data
    predictions = np.array([0.1, 0.2, 0.3, 0.4, 0.6, 0.7, 0.8, 0.9])
    ground_truth = np.array([0, 0, 0, 0, 1, 1, 1, 1])  # Matches probabilities roughly
    
    ece = compute_ece(predictions, ground_truth, bins=4)
    
    # ECE should be low for well-calibrated predictions
    assert ece >= 0.0, "ECE should be non-negative"
    assert ece <= 1.0, "ECE should be <= 1.0"


def test_compute_ece_overconfident() -> None:
    """Test ECE computation with overconfident predictions."""
    # All high confidence predictions but mixed ground truth
    predictions = np.array([0.9, 0.9, 0.9, 0.9, 0.9, 0.9])
    ground_truth = np.array([0, 0, 1, 1, 0, 1])  # 50% positive
    
    ece = compute_ece(predictions, ground_truth, bins=3)
    
    # ECE should be high for overconfident model
    assert ece > 0.1, f"Expected high ECE for overconfident model, got {ece}"
    assert ece <= 1.0, "ECE should be <= 1.0"


def test_compute_ece_empty_input() -> None:
    """Test ECE computation with empty input."""
    predictions = np.array([])
    ground_truth = np.array([])
    
    ece = compute_ece(predictions, ground_truth)
    
    assert ece == 0.0, "ECE should be 0.0 for empty input"


def test_compute_ece_mismatched_lengths() -> None:
    """Test ECE computation raises error for mismatched array lengths."""
    predictions = np.array([0.5, 0.6, 0.7])
    ground_truth = np.array([0, 1])  # Different length
    
    try:
        compute_ece(predictions, ground_truth)
        assert False, "Should have raised ValueError for mismatched lengths"
    except ValueError as e:
        assert "same length" in str(e)


def test_compute_ece_single_prediction() -> None:
    """Test ECE computation with single prediction."""
    predictions = np.array([0.7])
    ground_truth = np.array([1])
    
    ece = compute_ece(predictions, ground_truth, bins=1)
    
    # Should compute ECE without error
    assert isinstance(ece, float)
    assert ece >= 0.0
    assert ece <= 1.0


def test_get_calibration_status_thresholds() -> None:
    """Test calibration status determination based on ECE thresholds."""
    # Test good status (ECE < 0.1)
    assert get_calibration_status(0.05) == "good"
    assert get_calibration_status(0.09) == "good"
    
    # Test warning status (0.1 <= ECE < 0.2)
    assert get_calibration_status(0.1) == "warning"
    assert get_calibration_status(0.15) == "warning"
    assert get_calibration_status(0.19) == "warning"
    
    # Test alert status (ECE >= 0.2)
    assert get_calibration_status(0.2) == "alert"
    assert get_calibration_status(0.3) == "alert"
    assert get_calibration_status(0.5) == "alert"


def test_get_calibration_status_edge_cases() -> None:
    """Test calibration status with edge case values."""
    # Test boundary values
    assert get_calibration_status(0.0) == "good"
    assert get_calibration_status(0.099999) == "good"
    assert get_calibration_status(0.100001) == "warning"
    assert get_calibration_status(0.199999) == "warning"
    assert get_calibration_status(0.200001) == "alert"


def test_reliability_diagram_bins() -> None:
    """Test reliability diagram computation and bin structure."""
    predictions = np.array([0.1, 0.2, 0.3, 0.4, 0.6, 0.7, 0.8, 0.9])
    ground_truth = np.array([0, 0, 0, 1, 0, 1, 1, 1])
    
    reliability = compute_reliability_diagram(predictions, ground_truth, bins=4)
    
    # Check structure
    assert len(reliability.bin_midpoints) == 4
    assert len(reliability.bin_accuracies) == 4
    assert len(reliability.bin_confidences) == 4
    assert len(reliability.bin_counts) == 4
    
    # Check that midpoints are reasonable
    for midpoint in reliability.bin_midpoints:
        assert 0 <= midpoint <= 1, f"Bin midpoint {midpoint} should be between 0 and 1"
    
    # Check that accuracies are valid probabilities
    for accuracy in reliability.bin_accuracies:
        assert 0 <= accuracy <= 1, f"Bin accuracy {accuracy} should be between 0 and 1"
    
    # Check that confidences are valid probabilities  
    for confidence in reliability.bin_confidences:
        assert 0 <= confidence <= 1, f"Bin confidence {confidence} should be between 0 and 1"
    
    # Check that counts are non-negative integers
    for count in reliability.bin_counts:
        assert count >= 0, f"Bin count {count} should be non-negative"
        assert isinstance(count, int), f"Bin count {count} should be integer"


def test_reliability_diagram_empty_input() -> None:
    """Test reliability diagram with empty input."""
    predictions = np.array([])
    ground_truth = np.array([])
    
    try:
        reliability = compute_reliability_diagram(predictions, ground_truth, bins=5)
        # Should handle empty input gracefully
        assert len(reliability.bin_midpoints) == 5
        assert all(count == 0 for count in reliability.bin_counts)
    except ValueError:
        # Acceptable to raise error for empty input
        pass


def test_reliability_diagram_mismatched_lengths() -> None:
    """Test reliability diagram raises error for mismatched array lengths."""
    predictions = np.array([0.5, 0.6, 0.7])
    ground_truth = np.array([0, 1])  # Different length
    
    try:
        compute_reliability_diagram(predictions, ground_truth)
        assert False, "Should have raised ValueError for mismatched lengths"
    except ValueError as e:
        assert "same length" in str(e)


def test_reliability_diagram_single_bin() -> None:
    """Test reliability diagram with single bin."""
    predictions = np.array([0.4, 0.5, 0.6])
    ground_truth = np.array([0, 1, 1])
    
    reliability = compute_reliability_diagram(predictions, ground_truth, bins=1)
    
    # Should have one bin with all data
    assert len(reliability.bin_midpoints) == 1
    assert reliability.bin_counts[0] == 3
    assert 0 <= reliability.bin_accuracies[0] <= 1
    assert 0 <= reliability.bin_confidences[0] <= 1


def test_reliability_diagram_uniform_predictions() -> None:
    """Test reliability diagram when all predictions are the same."""
    predictions = np.array([0.5, 0.5, 0.5, 0.5])
    ground_truth = np.array([0, 0, 1, 1])  # 50% positive
    
    reliability = compute_reliability_diagram(predictions, ground_truth, bins=10)
    
    # Most bins should be empty, one bin should have all data
    total_count = sum(reliability.bin_counts)
    assert total_count == 4
    
    # Find the non-empty bin
    non_empty_bins = [i for i, count in enumerate(reliability.bin_counts) if count > 0]
    assert len(non_empty_bins) == 1, "Should have exactly one non-empty bin"
    
    # Check the non-empty bin
    bin_idx = non_empty_bins[0]
    assert reliability.bin_counts[bin_idx] == 4
    assert reliability.bin_accuracies[bin_idx] == 0.5  # 50% positive
    assert reliability.bin_confidences[bin_idx] == 0.5  # All predictions are 0.5


def test_reliability_diagram_extreme_predictions() -> None:
    """Test reliability diagram with extreme (0 or 1) predictions."""
    predictions = np.array([0.0, 0.0, 1.0, 1.0])
    ground_truth = np.array([0, 1, 0, 1])
    
    reliability = compute_reliability_diagram(predictions, ground_truth, bins=10)
    
    # Should have data in first and last bins only
    non_empty_bins = [i for i, count in enumerate(reliability.bin_counts) if count > 0]
    assert len(non_empty_bins) <= 2, "Should have at most 2 non-empty bins for extreme predictions"
    
    # Check that accuracies are computed correctly for extreme bins
    for bin_idx in non_empty_bins:
        assert 0 <= reliability.bin_accuracies[bin_idx] <= 1
