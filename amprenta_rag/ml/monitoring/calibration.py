"""Calibration monitoring service for ML classification models."""

from __future__ import annotations

from dataclasses import dataclass
from datetime import datetime, timezone
from typing import List, Optional
from uuid import UUID

import numpy as np

from amprenta_rag.database.models import MLModel, ModelMonitoringLog
from amprenta_rag.database.session import db_session


@dataclass
class ReliabilityDiagram:
    """Data for reliability diagram visualization."""
    
    bin_midpoints: List[float]  # [0.05, 0.15, ..., 0.95]
    bin_accuracies: List[float]  # Actual accuracy per bin
    bin_confidences: List[float]  # Mean predicted confidence per bin
    bin_counts: List[int]  # Number of samples per bin


@dataclass
class CalibrationReport:
    """Report containing calibration analysis results for a classification model."""
    
    model_id: UUID
    model_name: str
    status: str  # "good", "warning", "alert", "no_data", "not_applicable"
    ece: Optional[float]  # Expected Calibration Error
    reliability: Optional[ReliabilityDiagram]
    n_predictions_with_truth: int
    computed_at: datetime
    threshold_warning: float = 0.1
    threshold_alert: float = 0.2


def compute_ece(predictions: np.ndarray, ground_truth: np.ndarray, bins: int = 10) -> float:
    """
    Compute Expected Calibration Error (ECE) for binary classification.
    
    Args:
        predictions: Array of predicted probabilities (0-1)
        ground_truth: Array of binary ground truth labels (0/1)
        bins: Number of bins for calibration analysis
        
    Returns:
        ECE score (lower is better, 0 = perfectly calibrated)
    """
    if len(predictions) != len(ground_truth):
        raise ValueError("Predictions and ground truth must have the same length")
    
    if len(predictions) == 0:
        return 0.0
    
    # Create bins for confidence levels
    bin_boundaries = np.linspace(0, 1, bins + 1)
    bin_lowers = bin_boundaries[:-1]
    bin_uppers = bin_boundaries[1:]
    
    ece = 0.0
    
    for bin_lower, bin_upper in zip(bin_lowers, bin_uppers):
        # Find predictions in this confidence bin
        in_bin = (predictions > bin_lower) & (predictions <= bin_upper)
        prop_in_bin = in_bin.mean()
        
        if prop_in_bin > 0:
            # Compute accuracy and confidence for this bin
            accuracy_in_bin = ground_truth[in_bin].mean()
            avg_confidence_in_bin = predictions[in_bin].mean()
            
            # Add weighted difference to ECE
            ece += np.abs(avg_confidence_in_bin - accuracy_in_bin) * prop_in_bin
    
    return float(ece)


def compute_reliability_diagram(predictions: np.ndarray, ground_truth: np.ndarray, bins: int = 10) -> ReliabilityDiagram:
    """
    Compute reliability diagram data for calibration visualization.
    
    Args:
        predictions: Array of predicted probabilities (0-1)
        ground_truth: Array of binary ground truth labels (0/1)
        bins: Number of bins for calibration analysis
        
    Returns:
        ReliabilityDiagram with bin data for visualization
    """
    if len(predictions) != len(ground_truth):
        raise ValueError("Predictions and ground truth must have the same length")
    
    # Create bins for confidence levels
    bin_boundaries = np.linspace(0, 1, bins + 1)
    bin_lowers = bin_boundaries[:-1]
    bin_uppers = bin_boundaries[1:]
    
    # Calculate bin midpoints
    bin_midpoints = (bin_lowers + bin_uppers) / 2
    
    bin_accuracies = []
    bin_confidences = []
    bin_counts = []
    
    for bin_lower, bin_upper in zip(bin_lowers, bin_uppers):
        # Find predictions in this confidence bin
        in_bin = (predictions > bin_lower) & (predictions <= bin_upper)
        
        if in_bin.sum() > 0:
            # Compute accuracy and confidence for this bin
            accuracy_in_bin = ground_truth[in_bin].mean()
            avg_confidence_in_bin = predictions[in_bin].mean()
            count_in_bin = in_bin.sum()
        else:
            # Empty bin
            accuracy_in_bin = 0.0
            avg_confidence_in_bin = (bin_lower + bin_upper) / 2  # Use bin midpoint
            count_in_bin = 0
        
        bin_accuracies.append(float(accuracy_in_bin))
        bin_confidences.append(float(avg_confidence_in_bin))
        bin_counts.append(int(count_in_bin))
    
    return ReliabilityDiagram(
        bin_midpoints=bin_midpoints.tolist(),
        bin_accuracies=bin_accuracies,
        bin_confidences=bin_confidences,
        bin_counts=bin_counts
    )


def get_calibration_status(ece: float) -> str:
    """
    Determine calibration status based on ECE score.
    
    Args:
        ece: Expected Calibration Error score
        
    Returns:
        Status string: "good", "warning", or "alert"
    """
    if ece < 0.1:
        return "good"
    elif ece < 0.2:
        return "warning"
    else:
        return "alert"


def get_calibration_report(model_id: UUID) -> CalibrationReport:
    """
    Generate calibration report for a classification model.
    
    Args:
        model_id: UUID of the model to analyze
        
    Returns:
        CalibrationReport containing calibration analysis results
    """
    with db_session() as db:
        # Get model information
        model = db.query(MLModel).filter(MLModel.id == model_id).first()
        if not model:
            raise ValueError(f"Model with ID {model_id} not found")
        
        # Check if model is applicable for calibration analysis
        if model.model_type != "admet_classification":
            return CalibrationReport(
                model_id=model_id,
                model_name=model.name,
                status="not_applicable",
                ece=None,
                reliability=None,
                n_predictions_with_truth=0,
                computed_at=datetime.now(timezone.utc)
            )
        
        # Query for predictions with ground truth
        logs_with_truth = (
            db.query(ModelMonitoringLog)
            .filter(
                ModelMonitoringLog.model_id == model_id,
                ModelMonitoringLog.ground_truth.isnot(None)
            )
            .all()
        )
        
        if not logs_with_truth:
            return CalibrationReport(
                model_id=model_id,
                model_name=model.name,
                status="no_data",
                ece=None,
                reliability=None,
                n_predictions_with_truth=0,
                computed_at=datetime.now(timezone.utc)
            )
        
        # Extract predictions and ground truth
        predictions = []
        ground_truth = []
        
        for log in logs_with_truth:
            if log.prediction is not None and log.ground_truth is not None:
                # Ensure prediction is a probability (0-1)
                pred_value = float(log.prediction)
                if 0.0 <= pred_value <= 1.0:
                    predictions.append(pred_value)
                    # Convert ground truth to binary (0/1)
                    gt_value = float(log.ground_truth)
                    # Assume ground truth > 0.5 means positive class
                    binary_gt = 1.0 if gt_value > 0.5 else 0.0
                    ground_truth.append(binary_gt)
        
        if len(predictions) == 0:
            return CalibrationReport(
                model_id=model_id,
                model_name=model.name,
                status="no_data",
                ece=None,
                reliability=None,
                n_predictions_with_truth=0,
                computed_at=datetime.now(timezone.utc)
            )
        
        # Convert to numpy arrays
        predictions_array = np.array(predictions)
        ground_truth_array = np.array(ground_truth)
        
        # Compute ECE and reliability diagram
        ece = compute_ece(predictions_array, ground_truth_array)
        reliability = compute_reliability_diagram(predictions_array, ground_truth_array)
        
        # Determine status
        status = get_calibration_status(ece)
        
        return CalibrationReport(
            model_id=model_id,
            model_name=model.name,
            status=status,
            ece=ece,
            reliability=reliability,
            n_predictions_with_truth=len(predictions),
            computed_at=datetime.now(timezone.utc)
        )
