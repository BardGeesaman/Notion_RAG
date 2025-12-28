"""Tests for assay predictor functionality."""
import pytest
from unittest.mock import MagicMock, patch
import numpy as np
from uuid import uuid4

from amprenta_rag.analysis.assay_predictor import (
    train_assay_predictor,
    predict_assay_outcome,
    list_assay_models,
    TrainedModel,
    PredictionResult,
    TrainingDataStats,
)


class TestAssayPredictor:
    """Test assay predictor functionality."""

    def test_train_assay_predictor_success(self):
        """Test successful assay predictor training."""
        program_id = uuid4()
        mock_program = MagicMock()
        mock_program.id = program_id
        mock_program.name = "Test Program"
        
        mock_registered_model = MagicMock()
        mock_registered_model.id = uuid4()
        
        with patch('amprenta_rag.analysis.assay_predictor.db_session') as mock_session, \
             patch('amprenta_rag.analysis.assay_predictor._collect_training_data') as mock_collect, \
             patch('amprenta_rag.analysis.assay_predictor._prepare_training_data') as mock_prepare, \
             patch('amprenta_rag.analysis.assay_predictor._train_model') as mock_train, \
             patch('amprenta_rag.analysis.assay_predictor.get_registry') as mock_get_registry:
            
            # Mock database session
            mock_db = MagicMock()
            mock_db.query.return_value.filter.return_value.first.return_value = mock_program
            mock_session.return_value.__enter__.return_value = mock_db
            
            # Mock training data collection
            training_data = [("CCO", True), ("CC(=O)O", False)] * 30  # 60 compounds
            mock_collect.return_value = training_data
            
            # Mock data preparation with enough actives/inactives
            X = np.random.rand(60, 200)
            y = np.array([1] * 30 + [0] * 30)  # 30 actives, 30 inactives
            feature_names = [f"feat_{i}" for i in range(200)]
            stats = TrainingDataStats(60, 30, 30, 0.5, 200, 0.8)
            mock_prepare.return_value = (X, y, feature_names, stats)
            
            # Mock model training
            mock_model = MagicMock()
            performance = {"accuracy": 0.85, "roc_auc": 0.87}
            mock_train.return_value = (mock_model, performance)
            
            # Mock registry
            mock_registry = MagicMock()
            mock_registry.register_model.return_value = mock_registered_model
            mock_get_registry.return_value = mock_registry
            
            result = train_assay_predictor(
                program_id=program_id,
                assay_type="biochemical",
                min_actives=25,
                min_inactives=25
            )
            
            assert isinstance(result, TrainedModel)
            assert result.program_id == program_id
            assert result.assay_type == "biochemical"
            assert result.success is True

    def test_train_assay_predictor_program_not_found(self):
        """Test assay predictor training with program not found."""
        program_id = uuid4()
        
        with patch('amprenta_rag.analysis.assay_predictor.db_session') as mock_session:
            # Mock database session - program not found
            mock_db = MagicMock()
            mock_db.query.return_value.filter.return_value.first.return_value = None
            mock_session.return_value.__enter__.return_value = mock_db
            
            result = train_assay_predictor(
                program_id=program_id,
                assay_type="biochemical"
            )
            
            assert isinstance(result, TrainedModel)
            assert result.success is False
            assert "Program not found" in result.error_message

    def test_train_assay_predictor_no_training_data(self):
        """Test assay predictor training with no training data."""
        program_id = uuid4()
        mock_program = MagicMock()
        mock_program.id = program_id
        
        with patch('amprenta_rag.analysis.assay_predictor.db_session') as mock_session, \
             patch('amprenta_rag.analysis.assay_predictor._collect_training_data') as mock_collect:
            
            # Mock database session - program found
            mock_db = MagicMock()
            mock_db.query.return_value.filter.return_value.first.return_value = mock_program
            mock_session.return_value.__enter__.return_value = mock_db
            
            # Mock empty training data
            mock_collect.return_value = []
            
            result = train_assay_predictor(
                program_id=program_id,
                assay_type="biochemical"
            )
            
            assert isinstance(result, TrainedModel)
            assert result.success is False
            assert "No training data" in result.error_message

    def test_train_assay_predictor_insufficient_data(self):
        """Test assay predictor training with insufficient data."""
        program_id = uuid4()
        mock_program = MagicMock()
        mock_program.id = program_id
        
        with patch('amprenta_rag.analysis.assay_predictor.db_session') as mock_session, \
             patch('amprenta_rag.analysis.assay_predictor._collect_training_data') as mock_collect, \
             patch('amprenta_rag.analysis.assay_predictor._prepare_training_data') as mock_prepare:
            
            # Mock database session
            mock_db = MagicMock()
            mock_db.query.return_value.filter.return_value.first.return_value = mock_program
            mock_session.return_value.__enter__.return_value = mock_db
            
            # Mock some training data
            training_data = [("CCO", True), ("CC(=O)O", False)]
            mock_collect.return_value = training_data
            
            # Mock insufficient actives/inactives (only 5 each, need 50)
            X = np.random.rand(10, 200)
            y = np.array([1] * 5 + [0] * 5)
            feature_names = [f"feat_{i}" for i in range(200)]
            stats = TrainingDataStats(10, 5, 5, 0.5, 200, 0.3)
            mock_prepare.return_value = (X, y, feature_names, stats)
            
            result = train_assay_predictor(
                program_id=program_id,
                assay_type="biochemical",
                min_actives=50,
                min_inactives=50
            )
            
            assert isinstance(result, TrainedModel)
            assert result.success is False
            assert "Insufficient" in result.error_message

    def test_predict_assay_outcome_success(self):
        """Test successful assay outcome prediction."""
        model_id = uuid4()
        smiles_list = ["CCO", "CC(=O)O", "CCC"]
        
        with patch('amprenta_rag.analysis.assay_predictor.db_session') as mock_session, \
             patch('amprenta_rag.analysis.assay_predictor.pickle.loads') as mock_pickle_loads, \
             patch('amprenta_rag.analysis.assay_predictor._extract_features_for_prediction') as mock_extract:
            
            # Create mock model
            mock_model = MagicMock()
            mock_model.predict.return_value = np.array([1, 0, 1])
            mock_model.predict_proba.return_value = np.array([[0.2, 0.8], [0.7, 0.3], [0.1, 0.9]])
            mock_pickle_loads.return_value = mock_model
            
            # Mock MLModel from database
            mock_ml_model = MagicMock()
            mock_ml_model.model_data = b"mock_bytes"
            mock_ml_model.metadata = {"feature_names": ["feat1", "feat2", "feat3"]}
            
            # Mock database session
            mock_db = MagicMock()
            mock_db.query.return_value.filter.return_value.first.return_value = mock_ml_model
            mock_session.return_value.__enter__.return_value = mock_db
            
            # Mock feature extraction - returns numpy array only
            feature_matrix = np.array([[1, 0, 1], [0, 1, 0], [1, 1, 0]])
            mock_extract.return_value = feature_matrix
            
            results = predict_assay_outcome(model_id, smiles_list)
            
            assert len(results) == 3
            assert all(isinstance(r, PredictionResult) for r in results)
            
            # Check first prediction (active, 0.8 prob)
            assert results[0].compound_smiles == "CCO"
            assert results[0].prediction == "active"
            assert results[0].probability_active == 0.8
            assert results[0].confidence > 0.0

    def test_predict_assay_outcome_model_not_found(self):
        """Test prediction with model not found."""
        model_id = uuid4()
        smiles_list = ["CCO"]
        
        with patch('amprenta_rag.analysis.assay_predictor.db_session') as mock_session:
            # Mock database session - model not found
            mock_db = MagicMock()
            mock_db.query.return_value.filter.return_value.first.return_value = None
            mock_session.return_value.__enter__.return_value = mock_db
            
            results = predict_assay_outcome(model_id, smiles_list)
            
            assert len(results) == 1
            # Real function returns "inactive" with 0.0 probability on model not found
            assert results[0].prediction == "inactive"
            assert results[0].probability_active == 0.0

    def test_list_assay_models_success(self):
        """Test successful listing of assay models."""
        with patch('amprenta_rag.analysis.assay_predictor.db_session') as mock_session:
            # Create mock models
            mock_model_1 = MagicMock()
            mock_model_1.id = uuid4()
            mock_model_1.name = "Test Model 1"
            mock_model_1.version = "1.0"
            mock_model_1.created_at = "2025-01-01T00:00:00Z"
            mock_model_1.metadata = {"program_id": str(uuid4()), "assay_type": "biochemical"}
            mock_model_1.performance_metrics = {"accuracy": 0.85}
            
            mock_model_2 = MagicMock()
            mock_model_2.id = uuid4()
            mock_model_2.name = "Test Model 2"
            mock_model_2.version = "1.0"
            mock_model_2.created_at = "2025-01-02T00:00:00Z"
            mock_model_2.metadata = {"program_id": str(uuid4()), "assay_type": "hts"}
            mock_model_2.performance_metrics = {"accuracy": 0.9}
            
            # Mock database session
            mock_db = MagicMock()
            mock_db.query.return_value.filter.return_value.all.return_value = [mock_model_1, mock_model_2]
            mock_session.return_value.__enter__.return_value = mock_db
            
            models = list_assay_models()
            
            assert len(models) == 2
            assert models[0]["name"] == "Test Model 1"
            assert models[1]["name"] == "Test Model 2"

    def test_list_assay_models_empty(self):
        """Test listing when no models exist."""
        with patch('amprenta_rag.analysis.assay_predictor.db_session') as mock_session:
            # Mock database session - no models
            mock_db = MagicMock()
            mock_db.query.return_value.filter.return_value.all.return_value = []
            mock_session.return_value.__enter__.return_value = mock_db
            
            models = list_assay_models()
            
            assert len(models) == 0


class TestDataClasses:
    """Test data class functionality."""

    def test_training_data_stats_creation(self):
        """Test TrainingDataStats data class creation."""
        stats = TrainingDataStats(
            total_compounds=100,
            active_compounds=40,
            inactive_compounds=60,
            activity_rate=0.4,
            feature_count=200,
            data_quality_score=0.85
        )
        
        assert stats.total_compounds == 100
        assert stats.activity_rate == 0.4
        assert stats.data_quality_score == 0.85

    def test_trained_model_creation(self):
        """Test TrainedModel data class creation."""
        model_id = uuid4()
        program_id = uuid4()
        stats = TrainingDataStats(100, 40, 60, 0.4, 200, 0.85)
        performance = {"accuracy": 0.9, "precision": 0.85, "recall": 0.88}
        
        model = TrainedModel(
            model_id=model_id,
            program_id=program_id,
            assay_type="biochemical",
            model_performance=performance,
            training_stats=stats,
            feature_names=["feat1", "feat2"],
            training_time_seconds=120.5,
            success=True
        )
        
        assert model.model_id == model_id
        assert model.program_id == program_id
        assert model.success is True
        assert model.model_performance["accuracy"] == 0.9

    def test_prediction_result_creation(self):
        """Test PredictionResult data class creation."""
        result = PredictionResult(
            compound_smiles="CCO",
            prediction="active",
            probability_active=0.85,
            confidence=0.75,
            feature_vector=[1.0, 2.0, 3.0]
        )
        
        assert result.compound_smiles == "CCO"
        assert result.prediction == "active"
        assert result.probability_active == 0.85
        assert result.confidence == 0.75
        assert len(result.feature_vector) == 3


@pytest.mark.integration
class TestAssayPredictorIntegration:
    """Integration tests for assay predictor (requires RDKit)."""

    def test_real_feature_extraction(self):
        """Test real feature extraction with RDKit."""
        from amprenta_rag.analysis.assay_predictor import _extract_features_for_prediction
        
        smiles_list = ["CCO", "CC(=O)O", "c1ccccc1"]
        
        feature_matrix = _extract_features_for_prediction(smiles_list)
        
        assert feature_matrix is not None
        assert feature_matrix.shape[0] == 3  # 3 compounds
        assert feature_matrix.shape[1] > 0   # Some features extracted
