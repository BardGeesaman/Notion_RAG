"""Model monitoring API endpoints."""

from __future__ import annotations

from datetime import datetime, timezone
from typing import Dict
from uuid import UUID

from fastapi import APIRouter, Depends, HTTPException
from sqlalchemy.orm import Session

from amprenta_rag.api.schemas import (
    CalibrationReportSchema,
    DriftReportSchema,
    HealthReportSchema,
    MonitoringFeedbackRequest,
    MonitoringLogRequest,
)
from amprenta_rag.database.base import get_db
from amprenta_rag.database.models import MLModel, ModelMonitoringLog
from amprenta_rag.ml.monitoring.calibration import get_calibration_report
from amprenta_rag.ml.monitoring.drift import get_feature_drift

router = APIRouter(prefix="/monitoring", tags=["monitoring"])


@router.post(
    "/log",
    summary="Log model prediction for monitoring",
    response_model=Dict[str, str],
)
def log_prediction(
    request: MonitoringLogRequest,
    db: Session = Depends(get_db),
) -> Dict[str, str]:
    """Log a model prediction for monitoring and drift detection."""
    
    # Verify model exists
    model = db.query(MLModel).filter(MLModel.id == request.model_id).first()
    if not model:
        raise HTTPException(status_code=404, detail="Model not found")
    
    # Create monitoring log entry
    log_entry = ModelMonitoringLog(
        model_id=request.model_id,
        model_version=request.model_version,
        prediction_id=request.prediction_id,
        prediction=request.prediction,
        ground_truth=None,  # Will be updated via feedback
        input_hash=request.input_hash,
        feature_summary=request.feature_summary,
        created_at=datetime.now(timezone.utc),
    )
    
    db.add(log_entry)
    db.commit()
    
    return {
        "status": "logged",
        "prediction_id": str(request.prediction_id),
    }


@router.post(
    "/feedback",
    summary="Provide ground truth feedback for a prediction",
    response_model=Dict[str, str],
)
def provide_feedback(
    request: MonitoringFeedbackRequest,
    db: Session = Depends(get_db),
) -> Dict[str, str]:
    """Update a logged prediction with ground truth feedback."""
    
    # Find existing log entry
    log_entry = (
        db.query(ModelMonitoringLog)
        .filter(ModelMonitoringLog.prediction_id == request.prediction_id)
        .first()
    )
    
    if not log_entry:
        raise HTTPException(status_code=404, detail="Prediction ID not found")
    
    # Update with ground truth
    log_entry.ground_truth = request.ground_truth
    db.commit()
    
    return {
        "status": "updated",
        "prediction_id": str(request.prediction_id),
    }


@router.get(
    "/models/{model_id}/drift",
    summary="Get drift analysis report for a model",
    response_model=DriftReportSchema,
)
def get_drift_report(
    model_id: UUID,
    window_hours: int = 24,
    db: Session = Depends(get_db),
) -> DriftReportSchema:
    """Get feature drift analysis for a model over a specified time window."""
    
    try:
        drift_report = get_feature_drift(model_id, window_hours)
    except ValueError as e:
        raise HTTPException(status_code=404, detail=str(e))
    
    return DriftReportSchema(
        model_id=drift_report.model_id,
        model_name=drift_report.model_name,
        status=drift_report.status,
        psi_scores=drift_report.psi_scores,
        fp_aggregate_drift=drift_report.fp_aggregate_drift,
        window_hours=drift_report.window_hours,
        n_predictions=drift_report.n_predictions,
        computed_at=drift_report.computed_at,
    )


@router.get(
    "/models/{model_id}/calibration",
    summary="Get calibration analysis report for a classification model",
    response_model=CalibrationReportSchema,
)
def get_calibration_analysis(
    model_id: UUID,
    db: Session = Depends(get_db),
) -> CalibrationReportSchema:
    """Get calibration analysis for a classification model."""
    
    try:
        calibration_report = get_calibration_report(model_id)
    except ValueError as e:
        raise HTTPException(status_code=404, detail=str(e))
    
    return CalibrationReportSchema(
        model_id=calibration_report.model_id,
        model_name=calibration_report.model_name,
        status=calibration_report.status,
        ece=calibration_report.ece,
        n_predictions_with_truth=calibration_report.n_predictions_with_truth,
        computed_at=calibration_report.computed_at,
    )


@router.get(
    "/models/{model_id}/health",
    summary="Get overall health report combining drift and calibration",
    response_model=HealthReportSchema,
)
def get_model_health(
    model_id: UUID,
    window_hours: int = 24,
    db: Session = Depends(get_db),
) -> HealthReportSchema:
    """Get comprehensive model health report combining drift and calibration analysis."""
    
    # Get drift report
    try:
        drift_report = get_feature_drift(model_id, window_hours)
    except ValueError as e:
        raise HTTPException(status_code=404, detail=str(e))
    
    # Get calibration report
    try:
        calibration_report = get_calibration_report(model_id)
    except ValueError as e:
        raise HTTPException(status_code=404, detail=str(e))
    
    # Determine overall status (worst case)
    status_priority = {
        "good": 0,
        "no_baseline": 1,
        "no_data": 1,
        "not_applicable": 1,
        "warning": 2,
        "alert": 3,
    }
    
    drift_priority = status_priority.get(drift_report.status, 1)
    calibration_priority = status_priority.get(calibration_report.status, 1)
    
    if max(drift_priority, calibration_priority) >= 3:
        overall_status = "alert"
    elif max(drift_priority, calibration_priority) >= 2:
        overall_status = "warning"
    else:
        overall_status = "good"
    
    # Calculate max PSI across all features
    psi_max = None
    if drift_report.psi_scores:
        all_psi_values = list(drift_report.psi_scores.values())
        if drift_report.fp_aggregate_drift:
            all_psi_values.extend(drift_report.fp_aggregate_drift.values())
        if all_psi_values:
            psi_max = max(all_psi_values)
    
    return HealthReportSchema(
        model_id=model_id,
        model_name=drift_report.model_name,
        overall_status=overall_status,
        drift_status=drift_report.status,
        calibration_status=calibration_report.status,
        psi_max=psi_max,
        ece=calibration_report.ece,
        last_checked=datetime.now(timezone.utc),
    )
