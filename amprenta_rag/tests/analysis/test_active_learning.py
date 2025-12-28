"""Tests for active learning functionality."""
import pytest
from unittest.mock import MagicMock, patch
import numpy as np
from uuid import uuid4

from amprenta_rag.analysis.active_learning import (
    suggest_next_compounds,
    calculate_acquisition_scores,
    SuggestionResult,
    AcquisitionScore,
    _calculate_uncertainty_scores,
    _calculate_diversity_scores,
    _get_morgan_fingerprint,
    get_compound_suggestions_for_program,
    evaluate_strategy_performance,
)


class TestActiveLearning:
    """Test active learning functionality."""

    def test_suggest_next_compounds_uncertainty_success(self):
        """Test successful compound suggestion using uncertainty strategy."""
        screened = [
            {"compound_id": uuid4(), "smiles": "CCO", "activity": True},
            {"compound_id": uuid4(), "smiles": "CC(=O)O", "activity": False},
        ]
        
        candidates = [
            {"compound_id": uuid4(), "smiles": "CCC"},
            {"compound_id": uuid4(), "smiles": "CCCC"},
            {"compound_id": uuid4(), "smiles": "CCCCC"},
        ]
        
        model_id = uuid4()
        
        with patch('amprenta_rag.analysis.active_learning.calculate_acquisition_scores') as mock_scores:
            mock_scores.return_value = [
                AcquisitionScore(candidates[0]["compound_id"], 0.9, "High uncertainty"),
                AcquisitionScore(candidates[1]["compound_id"], 0.6, "Medium uncertainty"),
                AcquisitionScore(candidates[2]["compound_id"], 0.3, "Low uncertainty"),
            ]
            
            suggestions = suggest_next_compounds(
                screened=screened,
                candidates=candidates,
                strategy="uncertainty",
                batch_size=2,
                model_id=model_id
            )
            
            assert len(suggestions) == 2
            assert all(isinstance(s, SuggestionResult) for s in suggestions)
            
            # Should be ranked by acquisition score (highest first)
            assert suggestions[0].acquisition_score == 0.9
            assert suggestions[1].acquisition_score == 0.6
            assert suggestions[0].rank == 1
            assert suggestions[1].rank == 2

    def test_suggest_next_compounds_diversity_success(self):
        """Test successful compound suggestion using diversity strategy."""
        screened = [
            {"compound_id": uuid4(), "smiles": "CCO", "activity": True},
        ]
        
        candidates = [
            {"compound_id": uuid4(), "smiles": "c1ccccc1"},  # Aromatic - diverse
            {"compound_id": uuid4(), "smiles": "CCO"},       # Same as screened - not diverse
        ]
        
        with patch('amprenta_rag.analysis.active_learning.calculate_acquisition_scores') as mock_scores:
            mock_scores.return_value = [
                AcquisitionScore(candidates[0]["compound_id"], 0.8, "High diversity"),
                AcquisitionScore(candidates[1]["compound_id"], 0.2, "Low diversity"),
            ]
            
            suggestions = suggest_next_compounds(
                screened=screened,
                candidates=candidates,
                strategy="diversity",
                batch_size=2
            )
            
            assert len(suggestions) == 2
            assert suggestions[0].acquisition_score == 0.8  # Most diverse first
            assert suggestions[1].acquisition_score == 0.2

    def test_suggest_next_compounds_empty_candidates(self):
        """Test compound suggestion with empty candidates list."""
        screened = [{"compound_id": uuid4(), "smiles": "CCO", "activity": True}]
        candidates = []
        
        suggestions = suggest_next_compounds(
            screened=screened,
            candidates=candidates,
            strategy="diversity"
        )
        
        assert len(suggestions) == 0

    def test_calculate_acquisition_scores_uncertainty(self):
        """Test acquisition score calculation for uncertainty strategy."""
        candidates = [
            {"compound_id": uuid4(), "smiles": "CCC"},
            {"compound_id": uuid4(), "smiles": "CCCC"},
        ]
        model_id = uuid4()
        
        with patch('amprenta_rag.analysis.active_learning._calculate_uncertainty_scores') as mock_uncertainty:
            mock_uncertainty.return_value = [
                AcquisitionScore(candidates[0]["compound_id"], 0.8, "High uncertainty"),
                AcquisitionScore(candidates[1]["compound_id"], 0.4, "Low uncertainty"),
            ]
            
            scores = calculate_acquisition_scores(
                candidates=candidates,
                screened=[],
                strategy="uncertainty",
                model_id=model_id
            )
            
            assert len(scores) == 2
            assert scores[0].score == 0.8
            assert scores[1].score == 0.4

    def test_calculate_acquisition_scores_diversity(self):
        """Test acquisition score calculation for diversity strategy."""
        candidates = [{"compound_id": uuid4(), "smiles": "CCC"}]
        screened = [{"compound_id": uuid4(), "smiles": "CCO", "activity": True}]
        
        with patch('amprenta_rag.analysis.active_learning._calculate_diversity_scores') as mock_diversity:
            mock_diversity.return_value = [
                AcquisitionScore(candidates[0]["compound_id"], 0.7, "Moderate diversity"),
            ]
            
            scores = calculate_acquisition_scores(
                candidates=candidates,
                screened=screened,
                strategy="diversity"
            )
            
            assert len(scores) == 1
            assert scores[0].score == 0.7

    def test_calculate_acquisition_scores_unknown_strategy(self):
        """Test acquisition score calculation with unknown strategy."""
        candidates = [{"compound_id": uuid4(), "smiles": "CCC"}]
        
        scores = calculate_acquisition_scores(
            candidates=candidates,
            screened=[],
            strategy="unknown_strategy"
        )
        
        assert len(scores) == 0

    def test_calculate_uncertainty_scores_success(self):
        """Test uncertainty score calculation with model predictions."""
        candidates = [
            {"compound_id": uuid4(), "smiles": "CCC"},
            {"compound_id": uuid4(), "smiles": "CCCC"},
        ]
        model_id = uuid4()
        
        # Mock the prediction results
        mock_prediction1 = MagicMock()
        mock_prediction1.probability_active = 0.5  # Maximum uncertainty
        mock_prediction2 = MagicMock()
        mock_prediction2.probability_active = 0.9  # Low uncertainty
        
        with patch('amprenta_rag.analysis.active_learning.db_session') as mock_db_session, \
             patch('amprenta_rag.analysis.assay_predictor.predict_assay_outcome') as mock_predict:
            
            # Mock database
            mock_db = MagicMock()
            mock_db_session.return_value.__enter__.return_value = mock_db
            mock_model = MagicMock()
            mock_db.query.return_value.filter.return_value.first.return_value = mock_model
            
            # Mock predictions
            mock_predict.return_value = [mock_prediction1, mock_prediction2]
            
            scores = _calculate_uncertainty_scores(candidates, model_id)
            
            assert len(scores) == 2
            # First compound should have higher uncertainty (closer to 0.5 probability)
            assert scores[0].score > scores[1].score

    def test_calculate_uncertainty_scores_no_model(self):
        """Test uncertainty score calculation without model ID."""
        candidates = [{"compound_id": uuid4(), "smiles": "CCC"}]
        
        scores = _calculate_uncertainty_scores(candidates, model_id=None)
        
        assert len(scores) == 0

    def test_calculate_uncertainty_scores_model_not_found(self):
        """Test uncertainty score calculation with model not found."""
        candidates = [{"compound_id": uuid4(), "smiles": "CCC"}]
        model_id = uuid4()
        
        with patch('amprenta_rag.analysis.active_learning.db_session') as mock_db_session:
            mock_db = MagicMock()
            mock_db_session.return_value.__enter__.return_value = mock_db
            mock_db.query.return_value.filter.return_value.first.return_value = None
            
            scores = _calculate_uncertainty_scores(candidates, model_id)
            
            assert len(scores) == 0

    def test_calculate_uncertainty_scores_prediction_failure(self):
        """Test uncertainty score calculation with prediction failure."""
        candidates = [{"compound_id": uuid4(), "smiles": "CCC"}]
        model_id = uuid4()
        
        with patch('amprenta_rag.analysis.active_learning.db_session') as mock_db_session, \
             patch('amprenta_rag.analysis.assay_predictor.predict_assay_outcome') as mock_predict:
            
            mock_db = MagicMock()
            mock_db_session.return_value.__enter__.return_value = mock_db
            mock_model = MagicMock()
            mock_db.query.return_value.filter.return_value.first.return_value = mock_model
            
            # Mock prediction failure
            mock_predict.side_effect = Exception("Prediction failed")
            
            scores = _calculate_uncertainty_scores(candidates, model_id)
            
            assert len(scores) == 0

    def test_calculate_uncertainty_scores_no_probability_attribute(self):
        """Test uncertainty score calculation when predictions lack probability_active."""
        candidates = [{"compound_id": uuid4(), "smiles": "CCC"}]
        model_id = uuid4()
        
        mock_prediction = MagicMock()
        del mock_prediction.probability_active  # Remove the attribute
        
        with patch('amprenta_rag.analysis.active_learning.db_session') as mock_db_session, \
             patch('amprenta_rag.analysis.assay_predictor.predict_assay_outcome') as mock_predict:
            
            mock_db = MagicMock()
            mock_db_session.return_value.__enter__.return_value = mock_db
            mock_model = MagicMock()
            mock_db.query.return_value.filter.return_value.first.return_value = mock_model
            
            mock_predict.return_value = [mock_prediction]
            
            scores = _calculate_uncertainty_scores(candidates, model_id)
            
            assert len(scores) == 0

    def test_calculate_diversity_scores_success(self):
        """Test diversity score calculation with fingerprint similarity."""
        candidates = [
            {"compound_id": uuid4(), "smiles": "c1ccccc1"},    # Benzene - aromatic
            {"compound_id": uuid4(), "smiles": "CCO"},         # Ethanol - aliphatic
        ]
        screened = [
            {"compound_id": uuid4(), "smiles": "CCCO", "activity": True},  # Propanol - similar to ethanol
        ]
        
        with patch('amprenta_rag.analysis.active_learning._get_morgan_fingerprint') as mock_fp, \
             patch('amprenta_rag.analysis.active_learning.DataStructs') as mock_ds:
            
            # Mock fingerprints
            screened_fp = MagicMock()
            benzene_fp = MagicMock()
            ethanol_fp = MagicMock()
            
            mock_fp.side_effect = [screened_fp, benzene_fp, ethanol_fp]
            
            # Mock Tanimoto similarities
            # Benzene vs propanol = low similarity (high diversity)
            # Ethanol vs propanol = high similarity (low diversity)
            mock_ds.TanimotoSimilarity.side_effect = [0.2, 0.8]
            
            scores = _calculate_diversity_scores(candidates, screened)
            
            assert len(scores) == 2
            # Benzene should have higher diversity score (1 - 0.2 = 0.8)
            # Ethanol should have lower diversity score (1 - 0.8 = 0.2)
            assert abs(scores[0].score - 0.8) < 0.01
            assert abs(scores[1].score - 0.2) < 0.01

    def test_calculate_diversity_scores_no_screened(self):
        """Test diversity score calculation with no screened compounds."""
        candidates = [{"compound_id": uuid4(), "smiles": "CCC"}]
        screened = []
        
        scores = _calculate_diversity_scores(candidates, screened)
        
        assert len(scores) == 1
        # Should assign random score when no screened compounds
        assert 0.0 <= scores[0].score <= 1.0

    def test_get_morgan_fingerprint_success(self):
        """Test successful Morgan fingerprint generation."""
        smiles = "CCO"  # Ethanol
        
        # This would require RDKit, so we'll test the interface
        with patch('amprenta_rag.analysis.active_learning.Chem') as mock_chem, \
             patch('amprenta_rag.analysis.active_learning.rdMolDescriptors') as mock_desc:
            
            mock_mol = MagicMock()
            mock_chem.MolFromSmiles.return_value = mock_mol
            mock_fp = MagicMock()
            mock_desc.GetMorganFingerprintAsBitVect.return_value = mock_fp
            
            result = _get_morgan_fingerprint(smiles)
            
            assert result == mock_fp
            mock_chem.MolFromSmiles.assert_called_once_with(smiles)

    def test_get_morgan_fingerprint_invalid_smiles(self):
        """Test Morgan fingerprint generation with invalid SMILES."""
        smiles = "INVALID"
        
        with patch('amprenta_rag.analysis.active_learning.Chem') as mock_chem:
            mock_chem.MolFromSmiles.return_value = None  # Invalid SMILES
            
            result = _get_morgan_fingerprint(smiles)
            
            assert result is None

    def test_evaluate_strategy_performance(self):
        """Test strategy performance evaluation."""
        suggestions = [
            SuggestionResult(uuid4(), "CCO", 0.9, "uncertainty", 1, "High uncertainty"),
            SuggestionResult(uuid4(), "CC(=O)O", 0.7, "uncertainty", 2, "Medium uncertainty"),
            SuggestionResult(uuid4(), "CCC", 0.5, "uncertainty", 3, "Low uncertainty"),
        ]
        
        true_activities = {
            suggestions[0].compound_id: True,   # Hit
            suggestions[1].compound_id: False,  # Miss
            suggestions[2].compound_id: True,   # Hit
        }
        
        performance = evaluate_strategy_performance(suggestions, true_activities)
        
        assert performance["hit_rate"] == 2/3  # 2 hits out of 3 suggestions
        assert performance["total_suggested"] == 3
        assert performance["total_hits"] == 2
        assert performance["coverage"] == 1.0  # All suggestions have known activity


class TestDataClasses:
    """Test data class functionality."""

    def test_suggestion_result_creation(self):
        """Test SuggestionResult data class creation."""
        compound_id = uuid4()
        result = SuggestionResult(
            compound_id=compound_id,
            smiles="CCO",
            acquisition_score=0.85,
            strategy_used="uncertainty",
            rank=1,
            explanation="High uncertainty prediction"
        )
        
        assert result.compound_id == compound_id
        assert result.smiles == "CCO"
        assert result.acquisition_score == 0.85
        assert result.strategy_used == "uncertainty"
        assert result.rank == 1

    def test_acquisition_score_creation(self):
        """Test AcquisitionScore data class creation."""
        compound_id = uuid4()
        score = AcquisitionScore(
            compound_id=compound_id,
            score=0.75,
            explanation="Moderate diversity from screened compounds"
        )
        
        assert score.compound_id == compound_id
        assert score.score == 0.75
        assert "diversity" in score.explanation


@pytest.mark.integration
class TestActiveLearningIntegration:
    """Integration tests for active learning (requires RDKit and database)."""

    def test_real_fingerprint_generation(self):
        """Test real fingerprint generation with RDKit."""
        from amprenta_rag.analysis.active_learning import _get_morgan_fingerprint
        
        smiles_list = ["CCO", "CC(=O)O", "c1ccccc1"]
        
        fingerprints = [_get_morgan_fingerprint(smiles) for smiles in smiles_list]
        
        assert all(fp is not None for fp in fingerprints)
        # Test that different molecules have different fingerprints
        assert fingerprints[0] != fingerprints[2]  # Ethanol vs benzene should differ

    @pytest.mark.skip(reason="Requires database setup")
    def test_real_program_suggestions(self):
        """Test real compound suggestions for a program."""
        program_id = uuid4()
        
        suggestions = get_compound_suggestions_for_program(
            program_id=program_id,
            strategy="diversity",
            batch_size=5
        )
        
        assert isinstance(suggestions, list)
        if suggestions:  # If compounds exist for the program
            assert all(isinstance(s, SuggestionResult) for s in suggestions)
            assert len(suggestions) <= 5
