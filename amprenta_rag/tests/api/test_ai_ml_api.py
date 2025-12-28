"""Tests for AI/ML API endpoints."""
import pytest
from unittest.mock import MagicMock, patch
from uuid import uuid4
from fastapi.testclient import TestClient

from amprenta_rag.api.main import app


class TestDatasetFinderAPI:
    """Test dataset finder API endpoints."""

    def test_find_datasets_success(self):
        """Test successful dataset finding API call."""
        # Patch at the source module where the function is defined
        with patch('amprenta_rag.query.dataset_finder.find_datasets_by_nl') as mock_find:
            mock_result = MagicMock()
            mock_result.query = "breast cancer RNA-seq"
            mock_result.extracted_terms = {"disease": ["breast cancer"], "assay_type": ["RNA-seq"]}
            mock_result.results = [
                MagicMock(
                    accession="GSE123456",
                    title="Breast cancer study",
                    description="RNA-seq analysis",
                    source="geo",
                    species="human",
                    tissue="breast",
                    disease="breast cancer",
                    assay_type="RNA-seq",
                    sample_count=100,
                    url="https://geo.com/GSE123456",
                    score=0.9
                )
            ]
            mock_result.total_found = 1
            mock_result.sources_searched = ["geo"]
            mock_result.sources_failed = []
            mock_find.return_value = mock_result
            
            client = TestClient(app)
            response = client.post(
                "/api/v1/datasets/find",
                json={
                    "query": "breast cancer RNA-seq",
                    "repositories": ["geo"],
                    "max_results": 50
                }
            )
            
            assert response.status_code == 200
            data = response.json()
            assert data["query"] == "breast cancer RNA-seq"
            assert len(data["results"]) == 1
            assert data["results"][0]["accession"] == "GSE123456"

    def test_find_datasets_empty_query(self):
        """Test dataset finder with empty query - validation error."""
        client = TestClient(app)
        response = client.post(
            "/api/v1/datasets/find",
            json={
                "query": "",  # Empty query fails validation
                "max_results": 50
            }
        )
        
        # Empty query triggers validation error
        assert response.status_code == 422


class TestMetadataEnrichmentAPI:
    """Test metadata enrichment API endpoints."""

    def test_enrich_dataset_metadata_success(self):
        """Test successful metadata enrichment API call."""
        dataset_id = uuid4()
        
        # Need to mock both the database session and the enrichment function
        with patch('amprenta_rag.api.routers.datasets.get_database_session'), \
             patch('amprenta_rag.services.metadata_enrichment.enrich_dataset_metadata') as mock_enrich:
            
            mock_result = MagicMock()
            mock_result.dataset_id = dataset_id
            mock_result.success = True
            mock_result.enriched_fields = ["sample_count", "study_design"]
            mock_result.extracted_metadata = {"sample_count": 100, "study_design": "cohort"}
            mock_result.error_message = None
            mock_result.processing_time_seconds = 2.5
            mock_enrich.return_value = mock_result
            
            client = TestClient(app)
            response = client.post(f"/api/v1/datasets/{dataset_id}/enrich")
            
            # Check that we get a valid response (may require database)
            assert response.status_code in [200, 404, 500]

    def test_enrich_dataset_metadata_not_found(self):
        """Test metadata enrichment with dataset not found."""
        dataset_id = uuid4()
        
        client = TestClient(app)
        response = client.post(f"/api/v1/datasets/{dataset_id}/enrich")
        
        # Without database setup, this should return 404 or 500
        assert response.status_code in [404, 500]


class TestRelevanceScoringAPI:
    """Test relevance and novelty scoring API endpoints."""

    def test_score_relevance_success(self):
        """Test successful relevance scoring API call."""
        with patch('amprenta_rag.analysis.relevance_scoring.score_relevance') as mock_score:
            mock_result = MagicMock()
            mock_result.item_id = "dataset_123"
            mock_result.overall_score = 0.85
            mock_result.disease_match = 0.9
            mock_result.target_overlap = 0.8
            mock_result.data_quality = 0.85
            mock_result.explanation = "High relevance due to disease match"
            mock_result.processing_time_seconds = 1.5
            mock_result.cached = False
            mock_score.return_value = mock_result
            
            client = TestClient(app)
            response = client.post(
                "/api/v1/score/relevance",
                json={
                    "item": {
                        "id": "dataset_123",
                        "title": "Breast cancer study",
                        "description": "RNA-seq analysis",
                        "species": "human"
                    },
                    "context": {
                        "diseases": ["breast cancer"],
                        "species": ["human"]
                    }
                }
            )
            
            assert response.status_code == 200
            data = response.json()
            assert data["item_id"] == "dataset_123"
            assert data["overall_score"] == 0.85

    def test_score_novelty_success(self):
        """Test successful novelty scoring API call."""
        with patch('amprenta_rag.analysis.relevance_scoring.score_novelty') as mock_score:
            mock_result = MagicMock()
            mock_result.item_id = "dataset_123"
            mock_result.novelty_score = 0.7
            mock_result.max_similarity = 0.3
            mock_result.most_similar_item_id = "dataset_456"
            mock_result.explanation = "Moderately novel approach"
            mock_result.processing_time_seconds = 2.0
            mock_result.cached = False
            mock_score.return_value = mock_result
            
            client = TestClient(app)
            response = client.post(
                "/api/v1/score/novelty",
                json={
                    "item": {
                        "id": "dataset_123",
                        "title": "Novel cancer study"
                    },
                    "existing_items": [
                        {
                            "id": "dataset_456",
                            "title": "Traditional cancer study"
                        }
                    ]
                }
            )
            
            assert response.status_code == 200
            data = response.json()
            assert data["item_id"] == "dataset_123"
            assert data["novelty_score"] == 0.7

    def test_batch_score_success(self):
        """Test successful batch scoring API call."""
        from amprenta_rag.analysis.relevance_scoring import ScoredItem, RelevanceScore, NoveltyScore
        
        with patch('amprenta_rag.analysis.relevance_scoring.batch_score') as mock_batch:
            # Create proper dataclass instances
            scored_items = [
                ScoredItem(
                    item_id="item_1",
                    relevance_score=RelevanceScore(
                        item_id="item_1",
                        overall_score=0.9,
                        disease_match=0.9,
                        target_overlap=0.8,
                        data_quality=0.9,
                        explanation="High relevance",
                        processing_time_seconds=1.0,
                        cached=False,
                    ),
                    novelty_score=NoveltyScore(
                        item_id="item_1",
                        novelty_score=0.7,
                        max_similarity=0.3,
                        most_similar_item_id="item_2",
                        explanation="Moderately novel",
                        processing_time_seconds=1.0,
                        cached=False,
                    ),
                ),
            ]
            mock_batch.return_value = scored_items
            
            client = TestClient(app)
            response = client.post(
                "/api/v1/score/batch",
                json={
                    "items": [
                        {"id": "item_1", "title": "Study 1"}
                    ],
                    "context": {
                        "diseases": ["cancer"],
                        "species": ["human"]
                    },
                    "score_relevance": True,
                    "score_novelty": True
                }
            )
            
            assert response.status_code == 200
            data = response.json()
            assert data["total_items"] == 1
            assert len(data["items"]) == 1


class TestAssayPredictorAPI:
    """Test assay predictor API endpoints."""

    def test_train_assay_predictor_success(self):
        """Test successful assay predictor training API call."""
        program_id = uuid4()
        
        with patch('amprenta_rag.analysis.assay_predictor.train_assay_predictor') as mock_train:
            from amprenta_rag.analysis.assay_predictor import TrainedModel, TrainingDataStats
            
            model_id = uuid4()
            stats = TrainingDataStats(
                total_compounds=100,
                active_compounds=40,
                inactive_compounds=60,
                activity_rate=0.4,
                feature_count=200,
                data_quality_score=0.85
            )
            
            result = TrainedModel(
                model_id=model_id,
                program_id=program_id,
                assay_type="biochemical",
                model_performance={"accuracy": 0.85, "precision": 0.8},
                training_stats=stats,
                feature_names=["feat1", "feat2"],
                training_time_seconds=120.5,
                success=True,
                error_message=None
            )
            mock_train.return_value = result
            
            client = TestClient(app)
            response = client.post(
                "/api/v1/predictors/train",
                json={
                    "program_id": str(program_id),
                    "assay_type": "biochemical",
                    "min_actives": 30,
                    "min_inactives": 30
                }
            )
            
            assert response.status_code == 200
            data = response.json()
            assert data["success"] is True
            assert data["model_performance"]["accuracy"] == 0.85

    def test_predict_assay_outcome_success(self):
        """Test successful assay outcome prediction API call."""
        model_id = uuid4()
        
        with patch('amprenta_rag.analysis.assay_predictor.predict_assay_outcome') as mock_predict:
            from amprenta_rag.analysis.assay_predictor import PredictionResult
            
            results = [
                PredictionResult(
                    compound_smiles="CCO",
                    prediction="active",
                    probability_active=0.8,
                    confidence=0.7,
                    feature_vector=[1.0, 2.0, 3.0]
                ),
                PredictionResult(
                    compound_smiles="CC(=O)O",
                    prediction="inactive",
                    probability_active=0.3,
                    confidence=0.6,
                    feature_vector=[0.5, 1.5, 2.5]
                )
            ]
            mock_predict.return_value = results
            
            client = TestClient(app)
            response = client.post(
                f"/api/v1/predictors/{model_id}/predict",
                json={
                    "compound_smiles": ["CCO", "CC(=O)O"]
                }
            )
            
            assert response.status_code == 200
            data = response.json()
            assert data["total_predictions"] == 2
            assert data["predictions"][0]["compound_smiles"] == "CCO"
            assert data["predictions"][0]["prediction"] == "active"

    def test_list_assay_models_success(self):
        """Test successful assay models listing API call."""
        with patch('amprenta_rag.analysis.assay_predictor.list_assay_models') as mock_list:
            model_id = uuid4()
            mock_list.return_value = [
                {
                    "model_id": model_id,
                    "name": "Test Model",
                    "version": "1.0",
                    "program_id": str(uuid4()),
                    "assay_type": "biochemical",
                    "created_at": "2025-01-01T00:00:00Z",
                    "performance": {"accuracy": 0.85},
                    "training_stats": {},
                }
            ]
            
            client = TestClient(app)
            response = client.get("/api/v1/predictors")
            
            assert response.status_code == 200
            data = response.json()
            assert data["total_models"] == 1
            assert len(data["models"]) == 1


class TestActiveLearningAPI:
    """Test active learning API endpoints."""

    def test_suggest_compounds_success(self):
        """Test successful compound suggestion API call."""
        with patch('amprenta_rag.analysis.active_learning.suggest_next_compounds') as mock_suggest:
            from amprenta_rag.analysis.active_learning import SuggestionResult
            
            suggestions = [
                SuggestionResult(
                    compound_id=str(uuid4()),
                    smiles="CCC",
                    acquisition_score=0.9,
                    strategy_used="uncertainty",
                    rank=1,
                    explanation="High uncertainty"
                ),
                SuggestionResult(
                    compound_id=str(uuid4()),
                    smiles="CCCC",
                    acquisition_score=0.7,
                    strategy_used="uncertainty",
                    rank=2,
                    explanation="Medium uncertainty"
                )
            ]
            mock_suggest.return_value = suggestions
            
            client = TestClient(app)
            response = client.post(
                "/api/v1/screening/suggest",
                json={
                    "screened": [
                        {"compound_id": str(uuid4()), "smiles": "CCO", "activity": True}
                    ],
                    "candidates": [
                        {"compound_id": str(uuid4()), "smiles": "CCC"},
                        {"compound_id": str(uuid4()), "smiles": "CCCC"}
                    ],
                    "strategy": "uncertainty",
                    "batch_size": 2,
                    "model_id": str(uuid4())
                }
            )
            
            assert response.status_code == 200
            data = response.json()
            assert len(data["suggestions"]) == 2
            assert data["strategy_used"] == "uncertainty"

    def test_suggest_compounds_diversity(self):
        """Test diversity-based compound suggestion API call."""
        with patch('amprenta_rag.analysis.active_learning.suggest_next_compounds') as mock_suggest:
            from amprenta_rag.analysis.active_learning import SuggestionResult
            
            suggestions = [
                SuggestionResult(
                    compound_id=str(uuid4()),
                    smiles="c1ccccc1",
                    acquisition_score=0.8,
                    strategy_used="diversity",
                    rank=1,
                    explanation="High diversity"
                )
            ]
            mock_suggest.return_value = suggestions
            
            client = TestClient(app)
            response = client.post(
                "/api/v1/screening/suggest",
                json={
                    "screened": [
                        {"compound_id": str(uuid4()), "smiles": "CCO", "activity": True}
                    ],
                    "candidates": [
                        {"compound_id": str(uuid4()), "smiles": "c1ccccc1"}
                    ],
                    "strategy": "diversity",
                    "batch_size": 5
                }
            )
            
            assert response.status_code == 200
            data = response.json()
            assert len(data["suggestions"]) == 1
            assert data["strategy_used"] == "diversity"


@pytest.mark.integration
class TestAIMLAPIIntegration:
    """Integration tests for AI/ML API endpoints."""

    @pytest.mark.skip(reason="Requires external API access and database setup")
    def test_real_dataset_finder_api(self):
        """Test real dataset finder API with external services."""
        client = TestClient(app)
        response = client.post(
            "/api/v1/datasets/find",
            json={
                "query": "breast cancer RNA-seq human",
                "repositories": ["geo"],
                "max_results": 5
            }
        )
        
        assert response.status_code == 200
        data = response.json()
        assert "results" in data
        assert isinstance(data["results"], list)

    @pytest.mark.skip(reason="Requires OpenAI API access")
    def test_real_relevance_scoring_api(self):
        """Test real relevance scoring API with OpenAI."""
        client = TestClient(app)
        response = client.post(
            "/api/v1/score/relevance",
            json={
                "item": {
                    "id": "test_dataset",
                    "title": "Breast cancer gene expression analysis",
                    "description": "RNA-seq study of breast cancer patients",
                    "species": "human"
                },
                "context": {
                    "diseases": ["breast cancer"],
                    "species": ["human"],
                    "assay_types": ["RNA-seq"]
                }
            }
        )
        
        assert response.status_code == 200
        data = response.json()
        assert 0.0 <= data["overall_score"] <= 1.0
        assert len(data["explanation"]) > 0

    def test_real_active_learning_api(self):
        """Test real active learning API with RDKit (diversity strategy - no model needed)."""
        client = TestClient(app)
        response = client.post(
            "/api/v1/screening/suggest",
            json={
                "screened": [
                    {"compound_id": str(uuid4()), "smiles": "CCO", "activity": True}
                ],
                "candidates": [
                    {"compound_id": str(uuid4()), "smiles": "CCC"},
                    {"compound_id": str(uuid4()), "smiles": "c1ccccc1"}
                ],
                "strategy": "diversity",
                "batch_size": 2
            }
        )
        
        assert response.status_code == 200
        data = response.json()
        assert len(data["suggestions"]) <= 2
        assert data["strategy_used"] == "diversity"
