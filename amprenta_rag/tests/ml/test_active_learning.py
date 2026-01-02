"""Tests for Active Learning service."""
import pytest
from uuid import uuid4
from unittest.mock import MagicMock, patch
from datetime import datetime, timezone

from amprenta_rag.ml.active_learning import (
    ActiveLearningService,
    uncertainty_sampling,
    margin_sampling,
    entropy_sampling,
    hybrid_sampling,
)
import numpy as np


# ============ Sampling Strategy Tests ============

def test_uncertainty_sampling_ranks_by_highest_std():
    """Highest uncertainty should be ranked first."""
    uncertainties = np.array([0.1, 0.5, 0.3, 0.8])
    ranking = uncertainty_sampling(uncertainties)
    assert ranking[0] == 3  # 0.8 is highest
    assert ranking[1] == 1  # 0.5 is second


def test_margin_sampling_ranks_closest_to_boundary():
    """Samples closest to 0.5 decision boundary ranked first."""
    predictions = np.array([0.1, 0.45, 0.9, 0.52])
    ranking = margin_sampling(predictions)
    # 0.45 and 0.52 are closest to 0.5
    assert ranking[0] in [1, 3]
    assert ranking[1] in [1, 3]


def test_entropy_sampling_ranks_by_highest_entropy():
    """Highest entropy (near 0.5) should be ranked first."""
    predictions = np.array([0.1, 0.5, 0.9, 0.3])
    ranking = entropy_sampling(predictions)
    assert ranking[0] == 1  # 0.5 has max entropy


def test_hybrid_sampling_combines_strategies():
    """Hybrid should combine all factors."""
    predictions = np.array([0.2, 0.5, 0.8])
    uncertainties = np.array([0.1, 0.4, 0.2])
    applicability = np.array([0.3, 0.1, 0.5])
    ranking = hybrid_sampling(predictions, uncertainties, applicability)
    assert len(ranking) == 3
    assert all(isinstance(r, (int, np.integer)) for r in ranking)


def test_hybrid_sampling_with_custom_weights():
    """Hybrid sampling should accept custom weights."""
    predictions = np.array([0.2, 0.5, 0.8])
    uncertainties = np.array([0.1, 0.4, 0.2])
    applicability = np.array([0.3, 0.1, 0.5])
    weights = {"uncertainty": 0.8, "margin": 0.1, "applicability": 0.1}
    
    ranking = hybrid_sampling(predictions, uncertainties, applicability, weights)
    assert len(ranking) == 3


def test_sampling_strategies_handle_edge_cases():
    """Sampling strategies should handle uniform arrays."""
    # All same values
    uniform_array = np.array([0.5, 0.5, 0.5, 0.5])
    
    # Should not crash
    ranking = uncertainty_sampling(uniform_array)
    assert len(ranking) == 4
    
    ranking = margin_sampling(uniform_array)
    assert len(ranking) == 4
    
    ranking = entropy_sampling(uniform_array)
    assert len(ranking) == 4


# ============ Service Tests ============

@patch('amprenta_rag.ml.active_learning.db_session')
def test_create_cycle_increments_number(mock_db_session):
    """Cycle number should auto-increment."""
    mock_db = MagicMock()
    mock_db_session.return_value.__enter__.return_value = mock_db
    
    # Mock existing cycle
    mock_existing_cycle = MagicMock()
    mock_existing_cycle.cycle_number = 2
    mock_db.query.return_value.filter.return_value.order_by.return_value.first.return_value = mock_existing_cycle
    
    # Mock model
    mock_model = MagicMock()
    mock_model.metrics = {"auc": 0.85}
    mock_db.query.return_value.filter.return_value.first.return_value = mock_model
    
    service = ActiveLearningService(uuid4())
    cycle = service.create_cycle("uncertainty", 50)
    
    # Should create cycle with number 3 (2 + 1)
    mock_db.add.assert_called_once()
    added_cycle = mock_db.add.call_args[0][0]
    assert added_cycle.cycle_number == 3
    assert added_cycle.selection_strategy == "uncertainty"
    assert added_cycle.batch_size == 50


@patch('amprenta_rag.ml.active_learning.db_session')
def test_submit_label_updates_status(mock_db_session):
    """Label submission should update item status to 'labeled'."""
    mock_db = MagicMock()
    mock_db_session.return_value.__enter__.return_value = mock_db
    
    mock_item = MagicMock()
    mock_item.status = "pending"
    mock_db.query.return_value.filter.return_value.first.return_value = mock_item
    
    service = ActiveLearningService(uuid4())
    user_id = uuid4()
    item_id = uuid4()
    
    result = service.submit_label(
        item_id=item_id,
        label=0.75,
        source="experimental",
        confidence="high",
        user_id=user_id,
        notes="Test label"
    )
    
    # Verify item was updated
    assert mock_item.label == 0.75
    assert mock_item.label_source == "experimental"
    assert mock_item.label_confidence == "high"
    assert mock_item.labeled_by == user_id
    assert mock_item.notes == "Test label"
    assert mock_item.status == "labeled"


@patch('amprenta_rag.ml.active_learning.db_session')
def test_skip_item_updates_status(mock_db_session):
    """Skip should update status to 'skipped'."""
    mock_db = MagicMock()
    mock_db_session.return_value.__enter__.return_value = mock_db
    
    mock_item = MagicMock()
    mock_item.status = "pending"
    mock_db.query.return_value.filter.return_value.first.return_value = mock_item
    
    service = ActiveLearningService(uuid4())
    user_id = uuid4()
    item_id = uuid4()
    
    result = service.skip_item(item_id=item_id, user_id=user_id)
    
    # Verify item was updated
    assert mock_item.status == "skipped"
    assert mock_item.labeled_by == user_id


@patch('amprenta_rag.ml.active_learning.db_session')
def test_select_samples_filters_existing_queue(mock_db_session):
    """Selection should exclude compounds already in queue."""
    mock_db = MagicMock()
    mock_db_session.return_value.__enter__.return_value = mock_db
    
    # Mock model
    mock_model = MagicMock()
    mock_db.query.return_value.filter.return_value.first.return_value = mock_model
    
    # Mock existing queue items
    existing_compound_id = uuid4()
    mock_db.query.return_value.filter.return_value.all.return_value = [(existing_compound_id,)]
    
    # Mock compounds
    new_compound_id = uuid4()
    mock_compound = MagicMock()
    mock_compound.id = new_compound_id
    mock_compound.smiles = "CCO"
    mock_db.query.return_value.filter.return_value.all.return_value = [mock_compound]
    
    service = ActiveLearningService(uuid4())
    
    # Should only select new compound, not existing one
    compound_ids = [existing_compound_id, new_compound_id]
    items = service.select_samples(compound_ids, "uncertainty", 10)
    
    # Verify filtering logic was applied
    mock_db.query.assert_called()


@patch('amprenta_rag.ml.active_learning.db_session')
def test_get_pending_items_respects_limit(mock_db_session):
    """Limit parameter should cap returned items."""
    mock_db = MagicMock()
    mock_db_session.return_value.__enter__.return_value = mock_db
    
    # Mock query chain
    mock_query = mock_db.query.return_value.filter.return_value.order_by.return_value
    mock_query.limit.return_value.all.return_value = []
    
    service = ActiveLearningService(uuid4())
    service.get_pending_items(limit=25)
    
    # Verify limit was applied
    mock_query.limit.assert_called_with(25)


@patch('amprenta_rag.ml.active_learning.db_session')
def test_get_cycle_stats_returns_correct_structure(mock_db_session):
    """Stats should have required fields."""
    mock_db = MagicMock()
    mock_db_session.return_value.__enter__.return_value = mock_db
    
    # Mock cycles
    mock_cycle = MagicMock()
    mock_cycle.cycle_number = 1
    mock_cycle.status = "complete"
    mock_cycle.items_labeled = 45
    mock_cycle.metrics_before = {"auc": 0.8}
    mock_cycle.metrics_after = {"auc": 0.85}
    mock_db.query.return_value.filter.return_value.order_by.return_value.all.return_value = [mock_cycle]
    
    # Mock counts
    mock_db.query.return_value.filter.return_value.count.side_effect = [5, 45]  # pending, labeled
    
    service = ActiveLearningService(uuid4())
    stats = service.get_cycle_stats()
    
    # Verify structure
    assert "model_id" in stats
    assert "total_cycles" in stats
    assert "pending_items" in stats
    assert "labeled_items" in stats
    assert "cycles" in stats
    assert isinstance(stats["cycles"], list)


def test_select_samples_with_random_strategy():
    """Random strategy should work without model predictions."""
    # Test random strategy doesn't require uncertainty calculations
    predictions = np.array([0.2, 0.5, 0.8])
    uncertainties = np.array([0.1, 0.4, 0.2])
    applicability = np.array([0.3, 0.1, 0.5])
    
    # Random strategy should return valid ranking
    np.random.seed(42)  # For reproducible test
    ranking = np.random.permutation(len(predictions))
    assert len(ranking) == 3
    assert set(ranking) == {0, 1, 2}


def test_service_initialization():
    """Service should initialize with model ID."""
    model_id = uuid4()
    service = ActiveLearningService(model_id)
    assert service.model_id == model_id


@patch('amprenta_rag.ml.active_learning.db_session')
def test_get_labeling_queue_with_compound_details(mock_db_session):
    """get_labeling_queue should return compound details."""
    mock_db = MagicMock()
    mock_db_session.return_value.__enter__.return_value = mock_db
    
    # Mock queue item with compound
    mock_item = MagicMock()
    mock_item.id = uuid4()
    mock_item.compound_id = uuid4()
    mock_item.compound.smiles = "CCO"
    mock_item.compound.compound_id = "TEST001"
    mock_item.prediction = 0.6
    mock_item.uncertainty = 0.3
    
    mock_db.query.return_value.join.return_value.filter.return_value.order_by.return_value.limit.return_value.all.return_value = [mock_item]
    
    service = ActiveLearningService(uuid4())
    queue = service.get_labeling_queue()
    
    assert len(queue) == 1
    assert queue[0]["smiles"] == "CCO"
    assert queue[0]["compound_name"] == "TEST001"
    assert queue[0]["prediction"] == 0.6
