"""Program-specific ML models for assay outcome prediction."""
from __future__ import annotations

import logging
import pickle
from dataclasses import dataclass
from typing import Any, Dict, List, Optional, Tuple
from uuid import UUID, uuid4

import numpy as np
from sklearn.ensemble import RandomForestClassifier
from sklearn.metrics import accuracy_score, precision_score, recall_score, roc_auc_score
from sklearn.model_selection import train_test_split
from sqlalchemy.orm import Session

from amprenta_rag.database.models import BiochemicalResult, Compound, HTSCampaign, HTSResult, MLModel, Program
from amprenta_rag.database.session import db_session
from amprenta_rag.ml.registry import get_registry

logger = logging.getLogger(__name__)


@dataclass
class TrainingDataStats:
    """Statistics about training data."""
    
    total_compounds: int
    active_compounds: int
    inactive_compounds: int
    activity_rate: float
    feature_count: int
    data_quality_score: float


@dataclass
class TrainedModel:
    """Result of model training."""
    
    model_id: UUID
    program_id: UUID
    assay_type: str
    model_performance: Dict[str, float]
    training_stats: TrainingDataStats
    feature_names: List[str]
    training_time_seconds: float
    success: bool
    error_message: Optional[str] = None


@dataclass
class PredictionResult:
    """Result of assay prediction."""
    
    compound_smiles: str
    prediction: str  # "active" or "inactive"
    probability_active: float
    confidence: float  # Based on prediction probability distance from 0.5
    feature_vector: Optional[List[float]] = None


def train_assay_predictor(
    program_id: UUID,
    assay_type: str,
    features: Optional[List[str]] = None,
    min_actives: int = 50,
    min_inactives: int = 50,
) -> TrainedModel:
    """
    Train a program-specific ML model for assay outcome prediction.
    
    Args:
        program_id: UUID of the program to train model for
        assay_type: Type of assay to predict (e.g., "biochemical", "hts", "cell_viability")
        features: Optional list of specific features to use (defaults to molecular descriptors)
        min_actives: Minimum number of active compounds required for training
        min_inactives: Minimum number of inactive compounds required for training
        
    Returns:
        TrainedModel with training results and model metadata
    """
    import time
    start_time = time.time()
    
    logger.info(f"Starting assay predictor training for program {program_id}, assay_type {assay_type}")
    
    try:
        with db_session() as db:
            # Validate program exists
            program = db.query(Program).filter(Program.id == program_id).first()
            if not program:
                return TrainedModel(
                    model_id=uuid4(),
                    program_id=program_id,
                    assay_type=assay_type,
                    model_performance={},
                    training_stats=TrainingDataStats(0, 0, 0, 0.0, 0, 0.0),
                    feature_names=[],
                    training_time_seconds=time.time() - start_time,
                    success=False,
                    error_message="Program not found",
                )
            
            # Collect training data based on assay type
            training_data = _collect_training_data(db, program_id, assay_type)
            
            if not training_data:
                return TrainedModel(
                    model_id=uuid4(),
                    program_id=program_id,
                    assay_type=assay_type,
                    model_performance={},
                    training_stats=TrainingDataStats(0, 0, 0, 0.0, 0, 0.0),
                    feature_names=[],
                    training_time_seconds=time.time() - start_time,
                    success=False,
                    error_message="No training data found for the specified program and assay type",
                )
            
            # Prepare features and labels
            X, y, feature_names, training_stats = _prepare_training_data(training_data, features)
            
            # Check minimum data requirements
            n_actives = np.sum(y == 1)
            n_inactives = np.sum(y == 0)
            
            if n_actives < min_actives or n_inactives < min_inactives:
                return TrainedModel(
                    model_id=uuid4(),
                    program_id=program_id,
                    assay_type=assay_type,
                    model_performance={},
                    training_stats=training_stats,
                    feature_names=feature_names,
                    training_time_seconds=time.time() - start_time,
                    success=False,
                    error_message=f"Insufficient training data: {n_actives} actives, {n_inactives} inactives (need {min_actives}/{min_inactives})",
                )
            
            # Train model
            model, performance = _train_model(X, y)
            
            # Register model in the registry
            model_metadata = {
                "program_id": str(program_id),
                "assay_type": assay_type,
                "feature_names": feature_names,
                "training_stats": {
                    "total_compounds": training_stats.total_compounds,
                    "active_compounds": training_stats.active_compounds,
                    "inactive_compounds": training_stats.inactive_compounds,
                    "activity_rate": training_stats.activity_rate,
                    "feature_count": training_stats.feature_count,
                    "data_quality_score": training_stats.data_quality_score,
                },
                "performance": performance,
            }
            
            registry = get_registry()
            registered_model = registry.register_model(
                name=f"{program.name}_{assay_type}_predictor",
                version="1.0",
                model_type="assay_predictor",
                framework="sklearn",
                model_object=model,
                metrics=performance,
                description=f"Assay predictor for {program.name} {assay_type} assays",
            )
            
            # Store additional metadata
            registered_model.metadata = model_metadata
            db.commit()
            
            model_id = registered_model.id
            
            training_time = time.time() - start_time
            
            logger.info(
                f"Successfully trained assay predictor {model_id} for program {program_id}: "
                f"AUC={performance.get('roc_auc', 0):.3f}, training_time={training_time:.2f}s"
            )
            
            return TrainedModel(
                model_id=model_id,
                program_id=program_id,
                assay_type=assay_type,
                model_performance=performance,
                training_stats=training_stats,
                feature_names=feature_names,
                training_time_seconds=training_time,
                success=True,
            )
            
    except Exception as e:
        training_time = time.time() - start_time
        error_msg = f"Failed to train assay predictor: {str(e)}"
        logger.error(f"Error training assay predictor for program {program_id}: {e}", exc_info=True)
        
        return TrainedModel(
            model_id=uuid4(),
            program_id=program_id,
            assay_type=assay_type,
            model_performance={},
            training_stats=TrainingDataStats(0, 0, 0, 0.0, 0, 0.0),
            feature_names=[],
            training_time_seconds=training_time,
            success=False,
            error_message=error_msg,
        )


def predict_assay_outcome(
    model_id: UUID,
    compound_smiles: List[str],
) -> List[PredictionResult]:
    """
    Predict assay outcomes using a trained model.
    
    Args:
        model_id: UUID of the trained model
        compound_smiles: List of SMILES strings to predict
        
    Returns:
        List of PredictionResult objects
    """
    logger.info(f"Predicting assay outcomes for {len(compound_smiles)} compounds using model {model_id}")
    
    results = []
    
    try:
        with db_session() as db:
            # Load model from registry
            ml_model = db.query(MLModel).filter(MLModel.id == model_id).first()
            if not ml_model:
                logger.error(f"Model {model_id} not found in registry")
                return [
                    PredictionResult(
                        compound_smiles=smiles,
                        prediction="inactive",
                        probability_active=0.0,
                        confidence=0.0,
                    )
                    for smiles in compound_smiles
                ]
            
            # Load the actual model object
            model = pickle.loads(ml_model.model_data)
            
            # Get feature names from metadata
            metadata = ml_model.metadata or {}
            feature_names = metadata.get("feature_names", [])
            
            # Extract features for all compounds
            X = _extract_features_for_prediction(compound_smiles, feature_names)
            
            if X is None or X.shape[1] == 0:
                logger.error("Failed to extract features for prediction")
                return [
                    PredictionResult(
                        compound_smiles=smiles,
                        prediction="inactive",
                        probability_active=0.0,
                        confidence=0.0,
                    )
                    for smiles in compound_smiles
                ]
            
            # Make predictions
            try:
                probabilities = model.predict_proba(X)[:, 1]  # Probability of active class
                predictions = model.predict(X)
            except Exception as e:
                logger.error(f"Error during model prediction: {e}")
                return [
                    PredictionResult(
                        compound_smiles=smiles,
                        prediction="inactive",
                        probability_active=0.0,
                        confidence=0.0,
                    )
                    for smiles in compound_smiles
                ]
            
            # Convert to results
            for i, smiles in enumerate(compound_smiles):
                prob_active = float(probabilities[i])
                prediction = "active" if predictions[i] == 1 else "inactive"
                
                # Confidence based on distance from decision boundary (0.5)
                confidence = abs(prob_active - 0.5) * 2  # Scale to 0-1
                
                results.append(PredictionResult(
                    compound_smiles=smiles,
                    prediction=prediction,
                    probability_active=prob_active,
                    confidence=confidence,
                    feature_vector=X[i].tolist() if X is not None else None,
                ))
            
            logger.info(f"Successfully predicted outcomes for {len(results)} compounds")
            
    except Exception as e:
        logger.error(f"Error in predict_assay_outcome: {e}", exc_info=True)
        # Return default predictions on error
        results = [
            PredictionResult(
                compound_smiles=smiles,
                prediction="inactive",
                probability_active=0.0,
                confidence=0.0,
            )
            for smiles in compound_smiles
        ]
    
    return results


def list_assay_models(program_id: Optional[UUID] = None) -> List[Dict[str, Any]]:
    """
    List available assay prediction models.
    
    Args:
        program_id: Optional program ID to filter models
        
    Returns:
        List of model metadata dictionaries
    """
    try:
        with db_session() as db:
            query = db.query(MLModel).filter(MLModel.model_type == "assay_predictor")
            
            if program_id:
                # Filter by program_id in metadata
                models = query.all()
                filtered_models = []
                for model in models:
                    metadata = model.metadata or {}
                    if metadata.get("program_id") == str(program_id):
                        filtered_models.append(model)
                models = filtered_models
            else:
                models = query.all()
            
            results = []
            for model in models:
                metadata = model.metadata or {}
                performance = model.performance_metrics or {}
                
                results.append({
                    "model_id": model.id,
                    "name": model.name,
                    "version": model.version,
                    "program_id": metadata.get("program_id"),
                    "assay_type": metadata.get("assay_type"),
                    "created_at": model.created_at,
                    "performance": performance,
                    "training_stats": metadata.get("training_stats", {}),
                })
            
            return results
            
    except Exception as e:
        logger.error(f"Error listing assay models: {e}", exc_info=True)
        return []


def _collect_training_data(db: Session, program_id: UUID, assay_type: str) -> List[Tuple[str, bool]]:
    """Collect training data (SMILES, activity) from database."""
    training_data = []
    
    try:
        if assay_type.lower() in ["biochemical", "enzymatic"]:
            # Query BiochemicalResult table
            results = (
                db.query(BiochemicalResult, Compound.smiles)
                .join(Compound, BiochemicalResult.compound_id == Compound.id)
                .filter(BiochemicalResult.program_id == program_id)
                .filter(Compound.smiles.isnot(None))
                .all()
            )
            
            for result, smiles in results:
                # Define activity threshold (e.g., IC50 < 10 µM = active)
                is_active = False
                if result.ic50 is not None:
                    is_active = result.ic50 < 10.0  # 10 µM threshold
                elif result.ki is not None:
                    is_active = result.ki < 10.0
                elif result.ec50 is not None:
                    is_active = result.ec50 < 10.0
                
                training_data.append((smiles, is_active))
        
        elif assay_type.lower() in ["hts", "screening"]:
            # Query HTSResult table
            results = (
                db.query(HTSResult, Compound.smiles, HTSCampaign.program_id)
                .join(Compound, HTSResult.compound_id == Compound.id)
                .join(HTSCampaign, HTSResult.campaign_id == HTSCampaign.id)
                .filter(HTSCampaign.program_id == program_id)
                .filter(Compound.smiles.isnot(None))
                .filter(HTSResult.hit_flag.isnot(None))
                .all()
            )
            
            for result, smiles, _ in results:
                is_active = bool(result.hit_flag)
                training_data.append((smiles, is_active))
        
        else:
            logger.warning(f"Unknown assay type: {assay_type}")
        
        logger.info(f"Collected {len(training_data)} training examples for {assay_type}")
        return training_data
        
    except Exception as e:
        logger.error(f"Error collecting training data: {e}", exc_info=True)
        return []


def _prepare_training_data(
    training_data: List[Tuple[str, bool]],
    feature_names: Optional[List[str]] = None,
) -> Tuple[np.ndarray, np.ndarray, List[str], TrainingDataStats]:
    """Prepare features and labels for training."""
    smiles_list = [smiles for smiles, _ in training_data]
    labels = np.array([int(activity) for _, activity in training_data])
    
    # Extract features using ADMET predictor's feature extraction
    X = _extract_features_for_prediction(smiles_list, feature_names)
    
    if X is None or X.shape[1] == 0:
        # Return empty arrays if feature extraction failed
        return np.array([]), np.array([]), [], TrainingDataStats(0, 0, 0, 0.0, 0, 0.0)
    
    # Remove invalid samples (NaN features)
    valid_mask = ~np.isnan(X).any(axis=1)
    X = X[valid_mask]
    labels = labels[valid_mask]
    
    # Calculate training statistics
    n_actives = np.sum(labels == 1)
    n_inactives = np.sum(labels == 0)
    total_compounds = len(labels)
    activity_rate = n_actives / total_compounds if total_compounds > 0 else 0.0
    
    # Simple data quality score based on balance and size
    balance_score = min(n_actives, n_inactives) / max(n_actives, n_inactives) if max(n_actives, n_inactives) > 0 else 0.0
    size_score = min(total_compounds / 1000, 1.0)  # Normalize to 1000 compounds
    data_quality_score = (balance_score + size_score) / 2
    
    training_stats = TrainingDataStats(
        total_compounds=total_compounds,
        active_compounds=n_actives,
        inactive_compounds=n_inactives,
        activity_rate=activity_rate,
        feature_count=X.shape[1],
        data_quality_score=data_quality_score,
    )
    
    # Generate feature names if not provided
    if feature_names is None:
        feature_names = [f"feature_{i}" for i in range(X.shape[1])]
    
    return X, labels, feature_names, training_stats


def _extract_features_for_prediction(smiles_list: List[str], feature_names: Optional[List[str]] = None) -> Optional[np.ndarray]:
    """Extract features using ADMET predictor's feature extraction (per P1-2)."""
    try:
        # Import here to avoid circular imports
        from amprenta_rag.ml.admet.predictor import ADMETPredictor
        
        # Create a temporary ADMET predictor instance to use its feature extraction
        admet_predictor = ADMETPredictor()
        
        # Extract features using the existing method
        features_list = []
        for smiles in smiles_list:
            try:
                features = admet_predictor._get_features(smiles)
                if features is not None:
                    features_list.append(features)
                else:
                    # Use zero vector for invalid SMILES
                    features_list.append(np.zeros(200))  # Assuming 200 features
            except Exception as e:
                logger.warning(f"Failed to extract features for SMILES {smiles}: {e}")
                features_list.append(np.zeros(200))
        
        if not features_list:
            return None
        
        X = np.array(features_list)
        
        # Handle feature selection if specific features requested
        if feature_names is not None and len(feature_names) < X.shape[1]:
            # For now, just use the first N features
            # In production, this would map feature names to indices
            X = X[:, :len(feature_names)]
        
        return X
        
    except Exception as e:
        logger.error(f"Error extracting features: {e}", exc_info=True)
        return None


def _train_model(X: np.ndarray, y: np.ndarray) -> Tuple[RandomForestClassifier, Dict[str, float]]:
    """Train a random forest classifier and evaluate performance."""
    # Split data for training and validation
    X_train, X_test, y_train, y_test = train_test_split(
        X, y, test_size=0.2, random_state=42, stratify=y
    )
    
    # Train Random Forest model
    model = RandomForestClassifier(
        n_estimators=100,
        max_depth=10,
        min_samples_split=5,
        min_samples_leaf=2,
        random_state=42,
        n_jobs=-1,
    )
    
    model.fit(X_train, y_train)
    
    # Evaluate performance
    y_pred = model.predict(X_test)
    y_pred_proba = model.predict_proba(X_test)[:, 1]
    
    performance = {
        "accuracy": float(accuracy_score(y_test, y_pred)),
        "precision": float(precision_score(y_test, y_pred, zero_division=0)),
        "recall": float(recall_score(y_test, y_pred, zero_division=0)),
        "roc_auc": float(roc_auc_score(y_test, y_pred_proba)),
        "training_samples": len(X_train),
        "test_samples": len(X_test),
    }
    
    return model, performance


__all__ = [
    "TrainedModel",
    "PredictionResult",
    "TrainingDataStats",
    "train_assay_predictor",
    "predict_assay_outcome",
    "list_assay_models",
]
