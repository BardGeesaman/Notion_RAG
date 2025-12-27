"""Tests for multi-omics visualization analysis service."""

from __future__ import annotations

from types import SimpleNamespace
from uuid import uuid4

from amprenta_rag.analysis.multi_omics_viz import (
    _build_unified_feature_index,
    _normalize_feature_name,
    compute_alluvial_data,
    compute_upset_data,
)


def test_normalize_feature_name() -> None:
    """Test feature name normalization."""
    assert _normalize_feature_name("TP53") == "tp53"
    assert _normalize_feature_name("TP53_HUMAN") == "tp53"
    assert _normalize_feature_name("tp53_human") == "tp53"
    assert _normalize_feature_name("  BRCA1  ") == "brca1"
    assert _normalize_feature_name("BRCA1_MOUSE") == "brca1"
    assert _normalize_feature_name("") == ""
    assert _normalize_feature_name("AKT1_HUMAN_V2") == "akt1_human_v2"  # Only removes trailing _HUMAN/_MOUSE


def test_build_unified_feature_index() -> None:
    """Test building unified feature index from features."""
    # Mock features with different names but same normalized key
    features = [
        SimpleNamespace(
            name="TP53", 
            external_ids={"hgnc_symbol": "TP53", "uniprot_id": "P04637"}
        ),
        SimpleNamespace(
            name="tp53_human", 
            external_ids={}
        ),
        SimpleNamespace(
            name="BRCA1", 
            external_ids={"hgnc_symbol": "BRCA1"}
        ),
        SimpleNamespace(
            name="different_gene", 
            external_ids=None
        ),
    ]
    
    index = _build_unified_feature_index(features)
    
    # TP53 and tp53_human should be grouped together
    assert "tp53" in index
    assert len(index["tp53"]) == 2
    
    # BRCA1 should have its own group
    assert "brca1" in index
    assert len(index["brca1"]) == 1
    
    # different_gene should have its own group
    assert "different_gene" in index
    assert len(index["different_gene"]) == 1
    
    # External IDs should create additional keys (normalized to lowercase)
    assert "p04637" in index  # uniprot_id
    assert len(index["p04637"]) == 1


def test_compute_alluvial_data_structure(monkeypatch) -> None:
    """Test compute_alluvial_data returns correct structure."""
    # Mock dataset and feature data
    dataset1_id = uuid4()
    dataset2_id = uuid4()
    
    mock_datasets = [
        SimpleNamespace(id=dataset1_id, name="Dataset A", omics_type="transcriptomics"),
        SimpleNamespace(id=dataset2_id, name="Dataset B", omics_type="proteomics"),
    ]
    
    mock_features = [
        SimpleNamespace(name="TP53", external_ids={}),
        SimpleNamespace(name="BRCA1", external_ids={}),
        SimpleNamespace(name="AKT1", external_ids={}),
        SimpleNamespace(name="PIK3CA", external_ids={}),
        SimpleNamespace(name="EGFR", external_ids={}),
    ]
    
    # Mock database query results
    def mock_query_side_effect(*args, **kwargs):
        mock_query = SimpleNamespace()
        
        def mock_join(*args, **kwargs):
            return mock_query
        
        def mock_filter(*args, **kwargs):
            return mock_query
        
        def mock_all():
            if "Dataset" in str(args[0]):
                return mock_datasets
            elif "Feature" in str(args[0]):
                return mock_features
            return []
        
        mock_query.join = mock_join
        mock_query.filter = mock_filter
        mock_query.all = mock_all
        return mock_query
    
    mock_db = SimpleNamespace()
    mock_db.query = mock_query_side_effect
    
    result = compute_alluvial_data([dataset1_id, dataset2_id], mock_db)
    
    # Verify structure
    assert "nodes" in result
    assert "links" in result
    assert isinstance(result["nodes"], list)
    assert isinstance(result["links"], list)
    
    # Verify node structure
    if result["nodes"]:
        node = result["nodes"][0]
        assert "id" in node
        assert "label" in node
        assert "color" in node
    
    # Verify link structure
    if result["links"]:
        link = result["links"][0]
        assert "source" in link
        assert "target" in link
        assert "value" in link
        assert isinstance(link["source"], int)
        assert isinstance(link["target"], int)
        assert isinstance(link["value"], int)


def test_compute_upset_data_structure(monkeypatch) -> None:
    """Test compute_upset_data returns correct structure."""
    # Mock dataset and feature data
    dataset_ids = [uuid4(), uuid4(), uuid4()]
    
    mock_datasets = [
        SimpleNamespace(id=dataset_ids[0], name="Dataset A", omics_type="transcriptomics"),
        SimpleNamespace(id=dataset_ids[1], name="Dataset B", omics_type="proteomics"),
        SimpleNamespace(id=dataset_ids[2], name="Dataset C", omics_type="metabolomics"),
    ]
    
    mock_features = [
        SimpleNamespace(name="TP53", external_ids={}),
        SimpleNamespace(name="BRCA1", external_ids={}),
        SimpleNamespace(name="AKT1", external_ids={}),
    ]
    
    # Mock database query results
    def mock_query_side_effect(*args, **kwargs):
        mock_query = SimpleNamespace()
        
        def mock_join(*args, **kwargs):
            return mock_query
        
        def mock_filter(*args, **kwargs):
            return mock_query
        
        def mock_all():
            if "Dataset" in str(args[0]):
                return mock_datasets
            elif "Feature" in str(args[0]):
                return mock_features
            return []
        
        mock_query.join = mock_join
        mock_query.filter = mock_filter
        mock_query.all = mock_all
        return mock_query
    
    mock_db = SimpleNamespace()
    mock_db.query = mock_query_side_effect
    
    result = compute_upset_data(dataset_ids, mock_db)
    
    # Verify structure
    assert "sets" in result
    assert "intersections" in result
    assert "matrix" in result
    assert isinstance(result["sets"], list)
    assert isinstance(result["intersections"], list)
    assert isinstance(result["matrix"], list)
    
    # Verify sets structure
    if result["sets"]:
        set_item = result["sets"][0]
        assert "id" in set_item
        assert "name" in set_item
        assert "omics_type" in set_item
        assert "color" in set_item
        assert "size" in set_item
    
    # Verify intersections structure
    if result["intersections"]:
        intersection = result["intersections"][0]
        assert "sets" in intersection
        assert "count" in intersection
        assert "bitmask" in intersection
        assert isinstance(intersection["sets"], list)
        assert isinstance(intersection["count"], int)
        assert isinstance(intersection["bitmask"], int)
    
    # Verify matrix structure (list of dicts with presence vectors)
    if result["matrix"]:
        matrix_item = result["matrix"][0]
        assert "key" in matrix_item
        assert "presence" in matrix_item
        assert isinstance(matrix_item["presence"], list)
        assert len(matrix_item["presence"]) == len(dataset_ids)
        assert all(isinstance(x, int) and x in (0, 1) for x in matrix_item["presence"])


def test_compute_alluvial_data_empty_datasets() -> None:
    """Test alluvial computation with empty dataset list."""
    mock_db = SimpleNamespace()
    mock_db.query = lambda *args: SimpleNamespace(join=lambda *args: SimpleNamespace(filter=lambda *args: SimpleNamespace(all=lambda: [])))
    
    result = compute_alluvial_data([], mock_db)
    
    assert result == {"nodes": [], "links": []}


def test_compute_upset_data_empty_datasets() -> None:
    """Test upset computation with empty dataset list."""
    mock_db = SimpleNamespace()
    mock_db.query = lambda *args: SimpleNamespace(join=lambda *args: SimpleNamespace(filter=lambda *args: SimpleNamespace(all=lambda: [])))
    
    result = compute_upset_data([], mock_db)
    
    assert result == {"sets": [], "intersections": [], "matrix": []}
