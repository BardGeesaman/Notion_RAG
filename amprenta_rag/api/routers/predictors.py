"""API router for assay outcome predictors."""
from __future__ import annotations

import time
from typing import Optional
from uuid import UUID

from fastapi import APIRouter, HTTPException, Query

from amprenta_rag.api.schemas import (
    AssayModelResponse,
    AssayModelsListResponse,
    AssayPredictionRequest,
    AssayPredictionResponse,
    AssayPredictorTrainRequest,
    AssayPredictorTrainResponse,
    PredictionResultResponse,
    TrainingDataStatsResponse,
)

router = APIRouter()


@router.post("/train", response_model=AssayPredictorTrainResponse)
async def train_assay_predictor_endpoint(
    request: AssayPredictorTrainRequest,
) -> AssayPredictorTrainResponse:
    """
    Train a program-specific ML model for assay outcome prediction.
    
    Collects training data from BiochemicalResult/HTSCampaign tables,
    extracts molecular features, trains a Random Forest classifier,
    and stores the model in the registry.
    """
    try:
        from amprenta_rag.analysis.assay_predictor import train_assay_predictor
        
        # Train the model
        result = train_assay_predictor(
            program_id=request.program_id,
            assay_type=request.assay_type,
            features=request.features,
            min_actives=request.min_actives,
            min_inactives=request.min_inactives,
        )
        
        # Convert training stats to response format
        training_stats_response = TrainingDataStatsResponse(
            total_compounds=result.training_stats.total_compounds,
            active_compounds=result.training_stats.active_compounds,
            inactive_compounds=result.training_stats.inactive_compounds,
            activity_rate=result.training_stats.activity_rate,
            feature_count=result.training_stats.feature_count,
            data_quality_score=result.training_stats.data_quality_score,
        )
        
        return AssayPredictorTrainResponse(
            model_id=result.model_id,
            program_id=result.program_id,
            assay_type=result.assay_type,
            model_performance=result.model_performance,
            training_stats=training_stats_response,
            feature_names=result.feature_names,
            training_time_seconds=result.training_time_seconds,
            success=result.success,
            error_message=result.error_message,
        )
        
    except Exception as e:
        raise HTTPException(
            status_code=500,
            detail=f"Model training failed: {str(e)}"
        )


@router.post("/{model_id}/predict", response_model=AssayPredictionResponse)
async def predict_assay_outcome_endpoint(
    model_id: UUID,
    request: AssayPredictionRequest,
) -> AssayPredictionResponse:
    """
    Predict assay outcomes using a trained model.
    
    Loads the specified model from the registry, extracts features
    for the provided compounds, and returns predictions with
    probabilities and confidence scores.
    """
    try:
        from amprenta_rag.analysis.assay_predictor import predict_assay_outcome
        
        start_time = time.time()
        
        # Make predictions
        results = predict_assay_outcome(
            model_id=model_id,
            compound_smiles=request.compound_smiles,
        )
        
        # Convert results to response format
        prediction_responses = []
        for result in results:
            prediction_responses.append(PredictionResultResponse(
                compound_smiles=result.compound_smiles,
                prediction=result.prediction,
                probability_active=result.probability_active,
                confidence=result.confidence,
                feature_vector=result.feature_vector,
            ))
        
        processing_time = time.time() - start_time
        
        return AssayPredictionResponse(
            model_id=model_id,
            predictions=prediction_responses,
            total_predictions=len(prediction_responses),
            processing_time_seconds=processing_time,
        )
        
    except Exception as e:
        raise HTTPException(
            status_code=500,
            detail=f"Prediction failed: {str(e)}"
        )


@router.get("", response_model=AssayModelsListResponse)
async def list_assay_models_endpoint(
    program_id: Optional[UUID] = Query(None, description="Filter models by program ID"),
) -> AssayModelsListResponse:
    """
    List available assay prediction models.
    
    Optionally filter by program ID to show only models trained
    for a specific program.
    """
    try:
        from amprenta_rag.analysis.assay_predictor import list_assay_models
        
        # Get model list
        models_data = list_assay_models(program_id=program_id)
        
        # Convert to response format
        model_responses = []
        for model_data in models_data:
            model_responses.append(AssayModelResponse(
                model_id=model_data["model_id"],
                name=model_data["name"],
                version=model_data["version"],
                program_id=model_data["program_id"],
                assay_type=model_data["assay_type"],
                created_at=model_data["created_at"],
                performance=model_data["performance"],
                training_stats=model_data["training_stats"],
            ))
        
        return AssayModelsListResponse(
            models=model_responses,
            total_models=len(model_responses),
        )
        
    except Exception as e:
        raise HTTPException(
            status_code=500,
            detail=f"Failed to list models: {str(e)}"
        )
