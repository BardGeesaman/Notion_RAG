"""Tests for interactive plate analysis functions."""

from __future__ import annotations

from types import SimpleNamespace
from uuid import uuid4

from amprenta_rag.analysis.hts_qc import detect_plate_format, get_plate_heatmap_data


def test_get_plate_heatmap_data_include_all(monkeypatch) -> None:
    """Test get_plate_heatmap_data with include_all parameter."""
    from amprenta_rag.analysis import hts_qc
    
    # Mock database session and results
    test_campaign_id = uuid4()
    
    # Create mock results - some hits, some non-hits
    mock_results = [
        SimpleNamespace(
            well_position="A01",
            normalized_value=0.8,
            z_score=2.5,
            hit_flag=True,
            compound_id=uuid4(),
            result_id="result1"
        ),
        SimpleNamespace(
            well_position="A02", 
            normalized_value=0.2,
            z_score=0.5,
            hit_flag=False,
            compound_id=uuid4(),
            result_id="result2"
        ),
        SimpleNamespace(
            well_position="A03",
            normalized_value=0.9,
            z_score=3.0,
            hit_flag=True,
            compound_id=uuid4(),
            result_id="result3"
        )
    ]
    
    def mock_db_session():
        db = SimpleNamespace()
        
        def mock_query(model):
            query_obj = SimpleNamespace()
            
            def mock_filter(condition):
                # Simulate filtering for campaign_id
                filtered_obj = SimpleNamespace()
                
                def mock_filter_hit_flag(hit_condition):
                    # Return only hits when hit_flag filter is applied
                    hit_results = [r for r in mock_results if r.hit_flag]
                    return SimpleNamespace(all=lambda: hit_results)
                
                filtered_obj.filter = mock_filter_hit_flag
                filtered_obj.all = lambda: mock_results  # All results when no hit filter
                return filtered_obj
            
            query_obj.filter = mock_filter
            return query_obj
        
        db.query = mock_query
        return db
    
    # Mock the db_session context manager
    class MockDBSession:
        def __enter__(self):
            return mock_db_session()
        def __exit__(self, *args):
            pass
    
    monkeypatch.setattr(hts_qc, "db_session", MockDBSession)
    
    # Test include_all=False (hits only)
    wells_hits = get_plate_heatmap_data(test_campaign_id, include_all=False)
    assert len(wells_hits) == 2  # Only the 2 hits
    hit_positions = [w.well_position for w in wells_hits]
    assert "A01" in hit_positions
    assert "A03" in hit_positions
    assert "A02" not in hit_positions
    
    # Test include_all=True (all wells)
    wells_all = get_plate_heatmap_data(test_campaign_id, include_all=True)
    assert len(wells_all) == 3  # All 3 wells
    all_positions = [w.well_position for w in wells_all]
    assert "A01" in all_positions
    assert "A02" in all_positions
    assert "A03" in all_positions


def test_get_plate_heatmap_data_default_behavior(monkeypatch) -> None:
    """Test that default behavior preserves backward compatibility."""
    from amprenta_rag.analysis import hts_qc
    
    test_campaign_id = uuid4()
    
    # Mock with mixed hit/non-hit results
    mock_results = [
        SimpleNamespace(
            well_position="A01",
            normalized_value=0.8,
            z_score=2.5,
            hit_flag=True,
            compound_id=uuid4(),
            result_id="result1"
        ),
        SimpleNamespace(
            well_position="A02",
            normalized_value=0.2,
            z_score=0.5,
            hit_flag=False,
            compound_id=uuid4(),
            result_id="result2"
        )
    ]
    
    def mock_db_session():
        db = SimpleNamespace()
        
        def mock_query(model):
            query_obj = SimpleNamespace()
            
            def mock_filter(condition):
                filtered_obj = SimpleNamespace()
                
                def mock_filter_hit_flag(hit_condition):
                    # Return only hits when hit_flag filter is applied
                    hit_results = [r for r in mock_results if r.hit_flag]
                    return SimpleNamespace(all=lambda: hit_results)
                
                filtered_obj.filter = mock_filter_hit_flag
                return filtered_obj
            
            query_obj.filter = mock_filter
            return query_obj
        
        db.query = mock_query
        return db
    
    class MockDBSession:
        def __enter__(self):
            return mock_db_session()
        def __exit__(self, *args):
            pass
    
    monkeypatch.setattr(hts_qc, "db_session", MockDBSession)
    
    # Test default behavior (no include_all parameter)
    wells_default = get_plate_heatmap_data(test_campaign_id)
    assert len(wells_default) == 1  # Only hits by default
    assert wells_default[0].well_position == "A01"
    assert wells_default[0].hit_flag is True


def test_detect_plate_format_96() -> None:
    """Test 96-well plate format detection."""
    result = detect_plate_format(7, 11)
    assert result == "96-well"


def test_detect_plate_format_384() -> None:
    """Test 384-well plate format detection."""
    result = detect_plate_format(15, 23)
    assert result == "384-well"


def test_detect_plate_format_1536() -> None:
    """Test 1536-well plate format detection."""
    result = detect_plate_format(31, 47)
    assert result == "1536-well"


def test_detect_plate_format_custom() -> None:
    """Test custom plate format detection."""
    result = detect_plate_format(10, 10)
    assert result == "Custom (11×11)"


def test_detect_plate_format_edge_cases() -> None:
    """Test edge cases for plate format detection."""
    # Single well
    result = detect_plate_format(0, 0)
    assert result == "Custom (1×1)"
    
    # Rectangular custom format
    result = detect_plate_format(5, 7)
    assert result == "Custom (6×8)"
    
    # Large custom format
    result = detect_plate_format(20, 30)
    assert result == "Custom (21×31)"


def test_detect_plate_format_near_standard() -> None:
    """Test formats that are close to but not exactly standard."""
    # Almost 96-well but not quite
    result = detect_plate_format(7, 10)
    assert result == "Custom (8×11)"
    
    # Almost 384-well but not quite
    result = detect_plate_format(15, 22)
    assert result == "Custom (16×23)"
    
    # Almost 1536-well but not quite
    result = detect_plate_format(31, 46)
    assert result == "Custom (32×47)"
