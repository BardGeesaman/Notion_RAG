"""Drift detection service for ML model monitoring."""

from __future__ import annotations

from dataclasses import dataclass
from datetime import datetime, timedelta, timezone
from typing import Dict, List, Tuple
from uuid import UUID

import numpy as np

from amprenta_rag.database.models import MLModel, ModelMonitoringLog
from amprenta_rag.database.session import db_session


@dataclass
class DriftReport:
    """Report containing drift detection results for a model."""
    
    model_id: UUID
    model_name: str
    status: str  # "good", "warning", "alert", "no_baseline"
    psi_scores: Dict[str, float]  # {feature_name: psi_value}
    fp_aggregate_drift: Dict[str, float]  # {metric: psi}
    window_hours: int
    n_predictions: int
    computed_at: datetime
    threshold_warning: float = 0.1
    threshold_alert: float = 0.25


def compute_psi(reference_bins: List[float], current_bins: List[float]) -> float:
    """
    Compute Population Stability Index (PSI) between reference and current distributions.
    
    Args:
        reference_bins: Normalized bin counts from baseline/reference period
        current_bins: Normalized bin counts from current period
        
    Returns:
        PSI score (higher values indicate more drift)
    """
    if len(reference_bins) != len(current_bins):
        raise ValueError("Reference and current bins must have the same length")
    
    epsilon = 1e-6  # Small value to handle zero bins
    psi = 0.0
    
    for ref_pct, curr_pct in zip(reference_bins, current_bins):
        # Add epsilon to avoid log(0)
        ref_pct_adj = max(ref_pct, epsilon)
        curr_pct_adj = max(curr_pct, epsilon)
        
        psi += (curr_pct_adj - ref_pct_adj) * np.log(curr_pct_adj / ref_pct_adj)
    
    return float(psi)


def compute_feature_bins(values: List[float], n_bins: int = 10) -> Tuple[List[float], List[float]]:
    """
    Compute histogram bins and normalized counts for feature values.
    
    Args:
        values: List of feature values
        n_bins: Number of bins to create
        
    Returns:
        Tuple of (bin_edges, normalized_bin_counts)
    """
    if not values:
        return [], []
    
    # Convert to numpy array
    values_array = np.array(values)
    
    # Compute histogram
    counts, bin_edges = np.histogram(values_array, bins=n_bins)
    
    # Normalize counts to get percentages
    total_count = len(values)
    normalized_counts = counts / total_count if total_count > 0 else counts
    
    return bin_edges.tolist(), normalized_counts.tolist()


def get_drift_status(psi: float) -> str:
    """
    Determine drift status based on PSI score.
    
    Args:
        psi: Population Stability Index score
        
    Returns:
        Status string: "good", "warning", or "alert"
    """
    if psi < 0.1:
        return "good"
    elif psi < 0.25:
        return "warning"
    else:
        return "alert"


def compute_fingerprint_aggregates(fp_array: np.ndarray, centroid: np.ndarray) -> Dict[str, float]:
    """
    Compute aggregate metrics for molecular fingerprint drift detection.
    
    Args:
        fp_array: Array of molecular fingerprints (n_samples, n_bits)
        centroid: Centroid fingerprint from baseline period
        
    Returns:
        Dictionary with aggregate metrics
    """
    if len(fp_array) == 0:
        return {
            "pct_bits_set": 0.0,
            "tanimoto_to_centroid": 0.0,
            "hamming_distance_mean": 0.0,
            "hamming_distance_std": 0.0
        }
    
    # Percentage of bits set (average across all fingerprints)
    pct_bits_set = float(np.mean(np.sum(fp_array, axis=1) / fp_array.shape[1]))
    
    # Tanimoto similarity to centroid (average)
    tanimoto_scores = []
    for fp in fp_array:
        intersection = np.sum(fp & centroid)
        union = np.sum(fp | centroid)
        tanimoto = intersection / union if union > 0 else 0.0
        tanimoto_scores.append(tanimoto)
    
    tanimoto_to_centroid = float(np.mean(tanimoto_scores))
    
    # Hamming distance statistics
    hamming_distances = []
    for fp in fp_array:
        hamming_dist = np.sum(fp != centroid)
        hamming_distances.append(hamming_dist)
    
    hamming_distance_mean = float(np.mean(hamming_distances))
    hamming_distance_std = float(np.std(hamming_distances))
    
    return {
        "pct_bits_set": pct_bits_set,
        "tanimoto_to_centroid": tanimoto_to_centroid,
        "hamming_distance_mean": hamming_distance_mean,
        "hamming_distance_std": hamming_distance_std
    }


def get_feature_drift(model_id: UUID, window_hours: int = 24) -> DriftReport:
    """
    Compute feature drift for a model over a specified time window.
    
    Args:
        model_id: UUID of the model to analyze
        window_hours: Time window in hours for recent predictions
        
    Returns:
        DriftReport containing drift analysis results
    """
    with db_session() as db:
        # Get model information
        model = db.query(MLModel).filter(MLModel.id == model_id).first()
        if not model:
            raise ValueError(f"Model with ID {model_id} not found")
        
        # Check if baseline exists
        drift_baseline = None
        if model.metrics and "drift_baseline" in model.metrics:
            drift_baseline = model.metrics["drift_baseline"]
        
        if not drift_baseline:
            return DriftReport(
                model_id=model_id,
                model_name=model.name,
                status="no_baseline",
                psi_scores={},
                fp_aggregate_drift={},
                window_hours=window_hours,
                n_predictions=0,
                computed_at=datetime.now(timezone.utc)
            )
        
        # Get recent predictions within time window
        cutoff_time = datetime.now(timezone.utc) - timedelta(hours=window_hours)
        recent_logs = (
            db.query(ModelMonitoringLog)
            .filter(
                ModelMonitoringLog.model_id == model_id,
                ModelMonitoringLog.created_at >= cutoff_time
            )
            .all()
        )
        
        if not recent_logs:
            return DriftReport(
                model_id=model_id,
                model_name=model.name,
                status="good",
                psi_scores={},
                fp_aggregate_drift={},
                window_hours=window_hours,
                n_predictions=0,
                computed_at=datetime.now(timezone.utc)
            )
        
        # Extract feature summaries from recent logs
        current_features = {}
        fingerprints = []
        
        for log in recent_logs:
            if log.feature_summary:
                # Collect individual feature values
                for feature_name, feature_value in log.feature_summary.items():
                    if feature_name == "fingerprint":
                        # Handle fingerprint separately
                        if isinstance(feature_value, list):
                            fingerprints.append(np.array(feature_value, dtype=bool))
                    else:
                        # Regular numerical features
                        if feature_name not in current_features:
                            current_features[feature_name] = []
                        if isinstance(feature_value, (int, float)):
                            current_features[feature_name].append(float(feature_value))
        
        # Compute PSI for each feature
        psi_scores = {}
        for feature_name, current_values in current_features.items():
            if feature_name in drift_baseline and "bins" in drift_baseline[feature_name]:
                baseline_bins = drift_baseline[feature_name]["bins"]
                
                # Compute current bins using same edges as baseline
                if "bin_edges" in drift_baseline[feature_name]:
                    bin_edges = drift_baseline[feature_name]["bin_edges"]
                    counts, _ = np.histogram(current_values, bins=bin_edges)
                    current_bins = (counts / len(current_values)).tolist() if len(current_values) > 0 else []
                else:
                    # Fallback: compute new bins
                    _, current_bins = compute_feature_bins(current_values)
                
                if len(current_bins) == len(baseline_bins):
                    psi = compute_psi(baseline_bins, current_bins)
                    psi_scores[feature_name] = psi
        
        # Compute fingerprint aggregate drift
        fp_aggregate_drift = {}
        if fingerprints and "fingerprint_centroid" in drift_baseline:
            baseline_centroid = np.array(drift_baseline["fingerprint_centroid"], dtype=bool)
            fp_array = np.array(fingerprints)
            
            # Compute current aggregates
            current_aggregates = compute_fingerprint_aggregates(fp_array, baseline_centroid)
            
            # Compare with baseline aggregates
            if "fingerprint_aggregates" in drift_baseline:
                baseline_aggregates = drift_baseline["fingerprint_aggregates"]
                
                for metric_name, current_value in current_aggregates.items():
                    if metric_name in baseline_aggregates:
                        baseline_value = baseline_aggregates[metric_name]
                        # Compute PSI-like metric for aggregates (simplified)
                        if baseline_value > 0:
                            drift_ratio = abs(current_value - baseline_value) / baseline_value
                            fp_aggregate_drift[metric_name] = drift_ratio
        
        # Determine overall status (worst case across all features)
        all_psi_scores = list(psi_scores.values()) + list(fp_aggregate_drift.values())
        if not all_psi_scores:
            overall_status = "good"
        else:
            max_psi = max(all_psi_scores)
            overall_status = get_drift_status(max_psi)
        
        return DriftReport(
            model_id=model_id,
            model_name=model.name,
            status=overall_status,
            psi_scores=psi_scores,
            fp_aggregate_drift=fp_aggregate_drift,
            window_hours=window_hours,
            n_predictions=len(recent_logs),
            computed_at=datetime.now(timezone.utc)
        )
