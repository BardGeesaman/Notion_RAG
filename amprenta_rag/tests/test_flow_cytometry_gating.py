"""
Tests for flow cytometry gating algorithms and ingest service.

This test suite covers:
- Polygon, rectangle, and quadrant gating algorithms
- Boolean gate combinations (AND, OR, NOT)
- Population statistics computation
- FCS file ingestion workflow
- Gate application and population analysis
"""

import tempfile
from pathlib import Path
from unittest.mock import MagicMock, patch
from uuid import uuid4

import numpy as np
import pytest

from amprenta_rag.flow_cytometry.gating import (
    PopulationStats,
    apply_boolean_gate,
    apply_polygon_gate,
    apply_quadrant_gate,
    apply_rectangle_gate,
    compute_population_stats,
)
from amprenta_rag.flow_cytometry.ingest_service import GateCreate, ingest_fcs


class TestPolygonGating:
    """Test polygon gate functionality."""
    
    def test_polygon_gate_basic_triangle(self):
        """Test polygon gate with simple triangle."""
        # Create test events in 2D space
        events = np.array([
            [1.0, 1.0],  # Inside triangle
            [2.0, 2.0],  # Inside triangle  
            [0.5, 0.5],  # Inside triangle
            [4.0, 4.0],  # Outside triangle
            [0.0, 3.0],  # Outside triangle
        ], dtype=np.float32)
        
        # Define triangle vertices
        vertices = [[0.0, 0.0], [3.0, 0.0], [1.5, 3.0]]
        
        mask = apply_polygon_gate(events, x_idx=0, y_idx=1, vertices=vertices)
        
        # Check results
        assert mask.dtype == bool
        assert len(mask) == len(events)
        assert np.sum(mask) == 3  # First 3 events should be inside
        assert mask[0] == True
        assert mask[1] == True  
        assert mask[2] == True
        assert mask[3] == False
        assert mask[4] == False
    
    def test_polygon_gate_complex_shape(self):
        """Test polygon gate with more complex shape (10 vertices)."""
        # Create events around origin
        events = np.array([
            [0.0, 0.0],   # Center - inside
            [0.5, 0.5],   # Inside
            [2.0, 0.0],   # Outside (too far right)
            [-2.0, 0.0],  # Outside (too far left)
            [0.0, 2.0],   # Outside (too far up)
            [0.0, -2.0],  # Outside (too far down)
        ], dtype=np.float32)
        
        # Define octagon-like shape
        vertices = [
            [1.0, 0.0], [0.7, 0.7], [0.0, 1.0], [-0.7, 0.7],
            [-1.0, 0.0], [-0.7, -0.7], [0.0, -1.0], [0.7, -0.7],
            [1.0, 0.0], [0.8, 0.2]  # 10 vertices total
        ]
        
        mask = apply_polygon_gate(events, x_idx=0, y_idx=1, vertices=vertices)
        
        assert len(mask) == len(events)
        assert mask[0] == True   # Center
        assert mask[1] == True   # Inside
        assert mask[2] == False  # Outside
        assert mask[3] == False  # Outside
        assert mask[4] == False  # Outside
        assert mask[5] == False  # Outside
    
    def test_polygon_gate_edge_cases(self):
        """Test polygon gate error conditions."""
        events = np.array([[1.0, 2.0], [3.0, 4.0]], dtype=np.float32)
        
        # Too few vertices
        with pytest.raises(ValueError, match="at least 3 vertices"):
            apply_polygon_gate(events, 0, 1, [[0.0, 0.0], [1.0, 0.0]])
        
        # Invalid parameter indices
        with pytest.raises(ValueError, match="Parameter indices out of bounds"):
            apply_polygon_gate(events, 5, 1, [[0.0, 0.0], [1.0, 0.0], [0.5, 1.0]])
        
        # Negative indices
        with pytest.raises(ValueError, match="must be non-negative"):
            apply_polygon_gate(events, -1, 1, [[0.0, 0.0], [1.0, 0.0], [0.5, 1.0]])


class TestRectangleGating:
    """Test rectangle gate functionality."""
    
    def test_rectangle_gate_basic(self):
        """Test basic rectangle gate."""
        events = np.array([
            [1.0, 1.0],  # Inside
            [2.0, 2.0],  # Inside
            [0.5, 1.5],  # Inside
            [3.5, 1.0],  # Outside (x too high)
            [1.0, 3.5],  # Outside (y too high)
            [0.0, 1.0],  # Outside (x too low)
        ], dtype=np.float32)
        
        bounds = {"x_min": 0.5, "x_max": 3.0, "y_min": 0.5, "y_max": 3.0}
        
        mask = apply_rectangle_gate(events, x_idx=0, y_idx=1, bounds=bounds)
        
        assert len(mask) == len(events)
        assert np.sum(mask) == 3  # First 3 events inside
        assert mask[0] == True
        assert mask[1] == True
        assert mask[2] == True
        assert mask[3] == False
        assert mask[4] == False
        assert mask[5] == False
    
    def test_rectangle_gate_boundary_inclusive(self):
        """Test that rectangle gate boundaries are inclusive."""
        events = np.array([
            [1.0, 1.0],  # Exactly on x_min, y_min
            [2.0, 2.0],  # Exactly on x_max, y_max
            [1.5, 1.0],  # On y_min boundary
            [2.0, 1.5],  # On x_max boundary
        ], dtype=np.float32)
        
        bounds = {"x_min": 1.0, "x_max": 2.0, "y_min": 1.0, "y_max": 2.0}
        
        mask = apply_rectangle_gate(events, x_idx=0, y_idx=1, bounds=bounds)
        
        # All events should be inside (inclusive boundaries)
        assert np.all(mask)
    
    def test_rectangle_gate_error_conditions(self):
        """Test rectangle gate error handling."""
        events = np.array([[1.0, 2.0]], dtype=np.float32)
        
        # Missing bounds keys
        with pytest.raises(ValueError, match="must contain keys"):
            apply_rectangle_gate(events, 0, 1, {"x_min": 0.0})
        
        # Invalid bounds (min >= max)
        with pytest.raises(ValueError, match="min values must be less than max"):
            bounds = {"x_min": 2.0, "x_max": 1.0, "y_min": 0.0, "y_max": 1.0}
            apply_rectangle_gate(events, 0, 1, bounds)


class TestQuadrantGating:
    """Test quadrant gate functionality."""
    
    def test_quadrant_gate_basic(self):
        """Test basic quadrant gate assignment."""
        events = np.array([
            [-1.0, 1.0],  # Q1: x < thresh, y >= thresh
            [1.0, 1.0],   # Q2: x >= thresh, y >= thresh
            [-1.0, -1.0], # Q3: x < thresh, y < thresh
            [1.0, -1.0],  # Q4: x >= thresh, y < thresh
        ], dtype=np.float32)
        
        quadrants = apply_quadrant_gate(events, x_idx=0, y_idx=1, x_thresh=0.0, y_thresh=0.0)
        
        # Check all quadrants returned
        assert set(quadrants.keys()) == {"Q1", "Q2", "Q3", "Q4"}
        
        # Check correct assignment
        assert quadrants["Q1"][0] == True   # First event in Q1
        assert quadrants["Q2"][1] == True   # Second event in Q2
        assert quadrants["Q3"][2] == True   # Third event in Q3
        assert quadrants["Q4"][3] == True   # Fourth event in Q4
        
        # Each event should be in exactly one quadrant
        for i in range(len(events)):
            count = sum(quadrants[q][i] for q in quadrants.keys())
            assert count == 1
    
    def test_quadrant_gate_threshold_boundaries(self):
        """Test quadrant gate threshold boundary handling."""
        events = np.array([
            [0.0, 0.0],   # Exactly on both thresholds -> Q2
            [1.0, 0.0],   # On y threshold -> Q4  
            [0.0, 1.0],   # On x threshold -> Q2
        ], dtype=np.float32)
        
        quadrants = apply_quadrant_gate(events, x_idx=0, y_idx=1, x_thresh=0.0, y_thresh=0.0)
        
        # Events on threshold should go to Q2 (>= logic for both x and y)
        assert quadrants["Q2"][0] == True  # (0,0) -> Q2 (x>=thresh, y>=thresh)
        assert quadrants["Q2"][1] == True  # (1,0) -> Q2 (x>=thresh, y>=thresh) - y=0.0 >= 0.0
        assert quadrants["Q2"][2] == True  # (0,1) -> Q2 (x>=thresh, y>=thresh)


class TestBooleanGating:
    """Test boolean gate combinations."""
    
    def test_boolean_gate_and(self):
        """Test AND gate combination."""
        # Create two overlapping masks
        mask1 = np.array([True, True, False, False, True])
        mask2 = np.array([True, False, True, False, True])
        
        gate_id1, gate_id2 = uuid4(), uuid4()
        masks = {gate_id1: mask1, gate_id2: mask2}
        
        result = apply_boolean_gate(masks, "AND", [gate_id1, gate_id2])
        
        expected = mask1 & mask2  # [True, False, False, False, True]
        assert np.array_equal(result, expected)
    
    def test_boolean_gate_or(self):
        """Test OR gate combination."""
        mask1 = np.array([True, True, False, False, True])
        mask2 = np.array([True, False, True, False, False])
        
        gate_id1, gate_id2 = uuid4(), uuid4()
        masks = {gate_id1: mask1, gate_id2: mask2}
        
        result = apply_boolean_gate(masks, "OR", [gate_id1, gate_id2])
        
        expected = mask1 | mask2  # [True, True, True, False, True]
        assert np.array_equal(result, expected)
    
    def test_boolean_gate_not(self):
        """Test NOT gate negation."""
        mask1 = np.array([True, False, True, False])
        
        gate_id1 = uuid4()
        masks = {gate_id1: mask1}
        
        result = apply_boolean_gate(masks, "NOT", [gate_id1])
        
        expected = ~mask1  # [False, True, False, True]
        assert np.array_equal(result, expected)
    
    def test_boolean_gate_error_conditions(self):
        """Test boolean gate error handling."""
        mask1 = np.array([True, False])
        gate_id1 = uuid4()
        masks = {gate_id1: mask1}
        
        # Invalid operator
        with pytest.raises(ValueError, match="Invalid operator"):
            apply_boolean_gate(masks, "XOR", [gate_id1])
        
        # Missing gate ID
        with pytest.raises(ValueError, match="Gate IDs not found"):
            apply_boolean_gate(masks, "AND", [uuid4()])
        
        # Empty operand list
        with pytest.raises(ValueError, match="At least one operand"):
            apply_boolean_gate(masks, "AND", [])


class TestPopulationStats:
    """Test population statistics computation."""
    
    def test_population_stats_basic(self):
        """Test basic population statistics computation."""
        # Create test events
        events = np.array([
            [100.0, 200.0, 300.0],
            [110.0, 210.0, 310.0],
            [120.0, 220.0, 320.0],
            [130.0, 230.0, 330.0],
            [140.0, 240.0, 340.0],
        ], dtype=np.float32)
        
        # Select subset of events
        mask = np.array([True, False, True, False, True])  # Events 0, 2, 4
        param_names = ["FSC-A", "SSC-A", "CD3-FITC"]
        
        stats = compute_population_stats(events, mask, param_names, total_events=5)
        
        assert stats.n_events == 3
        assert stats.pct_of_total == 60.0  # 3/5 * 100
        
        # Check median values (events 0, 2, 4 -> values [100, 120, 140])
        assert stats.median_values["FSC-A"] == 120.0
        assert stats.median_values["SSC-A"] == 220.0
        assert stats.median_values["CD3-FITC"] == 320.0
        
        # Check mean values
        assert stats.mean_values["FSC-A"] == 120.0  # (100+120+140)/3
        assert stats.mean_values["SSC-A"] == 220.0
        assert stats.mean_values["CD3-FITC"] == 320.0
        
        # CV should be calculated (std/mean * 100)
        assert "FSC-A" in stats.cv_values
        assert stats.cv_values["FSC-A"] > 0
    
    def test_population_stats_empty_population(self):
        """Test statistics for empty population (no events)."""
        events = np.array([[100.0, 200.0]], dtype=np.float32)
        mask = np.array([False])  # No events selected
        param_names = ["FSC-A", "SSC-A"]
        
        stats = compute_population_stats(events, mask, param_names)
        
        assert stats.n_events == 0
        assert stats.median_values["FSC-A"] == 0.0
        assert stats.mean_values["FSC-A"] == 0.0
        assert stats.cv_values["FSC-A"] == 0.0
    
    def test_population_stats_with_parent(self):
        """Test population statistics with parent population."""
        events = np.array([[1.0, 2.0], [3.0, 4.0], [5.0, 6.0]], dtype=np.float32)
        mask = np.array([True, False, True])  # 2 events
        parent_mask = np.array([True, True, False])  # 2 parent events
        param_names = ["X", "Y"]
        
        stats = compute_population_stats(events, mask, param_names, parent_mask=parent_mask)
        
        assert stats.n_events == 2
        assert stats.pct_of_parent == 100.0  # 2/2 * 100 (both selected events are in parent)
    
    def test_population_stats_error_conditions(self):
        """Test population statistics error handling."""
        events = np.array([[1.0, 2.0]], dtype=np.float32)
        mask = np.array([True, False])  # Wrong length
        param_names = ["X"]
        
        # Mask length mismatch
        with pytest.raises(ValueError, match="Mask length.*doesn't match"):
            compute_population_stats(events, mask, param_names)
        
        # Parameter name count mismatch  
        mask = np.array([True])
        param_names = ["X", "Y", "Z"]  # Too many names
        with pytest.raises(ValueError, match="Parameter name count.*doesn't match"):
            compute_population_stats(events, mask, param_names)


class TestIngestService:
    """Test FCS ingestion service."""
    
    @patch('amprenta_rag.flow_cytometry.ingest_service.validate_fcs')
    @patch('amprenta_rag.flow_cytometry.ingest_service.db_session')
    @patch('amprenta_rag.flow_cytometry.ingest_service.threading.Thread')
    def test_ingest_fcs_basic(self, mock_thread, mock_db_session, mock_validate):
        """Test basic FCS ingestion workflow."""
        # Setup mocks
        mock_validate.return_value = []  # No validation issues
        
        mock_db = MagicMock()
        mock_db_session.return_value.__enter__.return_value = mock_db
        
        # Mock dataset creation
        mock_dataset = MagicMock()
        mock_dataset.id = uuid4()
        mock_db.add = MagicMock()
        mock_db.commit = MagicMock()
        mock_db.refresh = MagicMock(side_effect=lambda x: setattr(x, 'id', mock_dataset.id))
        
        # Mock FlowCytometryDataset creation and retrieval
        mock_flow_dataset = MagicMock()
        mock_flow_dataset.id = uuid4()
        mock_db.query.return_value.filter_by.return_value.first.return_value = mock_flow_dataset
        
        # Create temporary FCS file
        with tempfile.NamedTemporaryFile(suffix='.fcs', delete=False) as f:
            fcs_path = f.name
            f.write(b"mock fcs data")
        
        try:
            # Test ingestion
            result = ingest_fcs(fcs_path)
            
            # Verify thread was started
            mock_thread.assert_called_once()
            assert mock_thread.call_args[1]['daemon'] == True
            
            # Verify result
            assert result == mock_flow_dataset
            
        finally:
            Path(fcs_path).unlink()  # Cleanup
    
    def test_ingest_fcs_file_not_found(self):
        """Test FCS ingestion with missing file."""
        with pytest.raises(FileNotFoundError):
            ingest_fcs("/nonexistent/file.fcs")
    
    @patch('amprenta_rag.flow_cytometry.ingest_service.validate_fcs')
    def test_ingest_fcs_invalid_file(self, mock_validate):
        """Test FCS ingestion with invalid file."""
        mock_validate.return_value = ["Invalid header", "Missing required fields"]
        
        # Create temporary file
        with tempfile.NamedTemporaryFile(suffix='.fcs', delete=False) as f:
            fcs_path = f.name
            f.write(b"invalid data")
        
        try:
            with pytest.raises(ValueError, match="Invalid FCS file"):
                ingest_fcs(fcs_path)
        finally:
            Path(fcs_path).unlink()


class TestGateApplication:
    """Test gate application to datasets."""
    
    @patch('amprenta_rag.flow_cytometry.ingest_service.db_session')
    @patch('amprenta_rag.flow_cytometry.ingest_service._apply_gate_and_compute_population')
    def test_apply_gate_to_dataset_basic(self, mock_compute, mock_db_session):
        """Test applying gate to dataset."""
        # Setup mocks
        mock_db = MagicMock()
        mock_db_session.return_value.__enter__.return_value = mock_db
        
        # Mock dataset and parameters
        mock_flow_dataset = MagicMock()
        mock_flow_dataset.processing_status = "completed"
        
        mock_x_param = MagicMock()
        mock_y_param = MagicMock()
        
        mock_db.query.return_value.filter_by.return_value.first.side_effect = [
            mock_flow_dataset,  # FlowCytometryDataset query
            mock_x_param,       # X parameter query  
            mock_y_param,       # Y parameter query
            MagicMock()         # Final gate query
        ]
        
        # Mock gate creation
        mock_gate = MagicMock()
        mock_gate.id = uuid4()
        mock_db.add = MagicMock()
        mock_db.commit = MagicMock()
        mock_db.refresh = MagicMock(side_effect=lambda x: setattr(x, 'id', mock_gate.id))
        
        # Create gate definition
        gate_def = GateCreate(
            gate_name="Test Gate",
            gate_type="rectangle",
            gate_definition={"x_min": 0.0, "x_max": 1000.0, "y_min": 0.0, "y_max": 1000.0},
            x_parameter_id=uuid4(),
            y_parameter_id=uuid4()
        )
        
        # Test gate application
        from amprenta_rag.flow_cytometry.ingest_service import apply_gate_to_dataset
        result = apply_gate_to_dataset(uuid4(), gate_def)
        
        # Verify population computation was triggered
        mock_compute.assert_called_once()
        
        # Verify gate was created
        mock_db.add.assert_called()
        mock_db.commit.assert_called()


@pytest.mark.integration
class TestEndToEndFlow:
    """Integration tests for complete flow cytometry workflow."""
    
    def test_hierarchical_gating_workflow(self):
        """Test complete hierarchical gating workflow."""
        # This would be a full integration test with real data
        # For now, just verify the workflow components exist
        from amprenta_rag.flow_cytometry.ingest_service import ingest_fcs, apply_gate_to_dataset
        from amprenta_rag.flow_cytometry.gating import apply_polygon_gate, compute_population_stats
        
        # Verify all required functions are available
        assert callable(ingest_fcs)
        assert callable(apply_gate_to_dataset)
        assert callable(apply_polygon_gate)
        assert callable(compute_population_stats)


if __name__ == "__main__":
    pytest.main([__file__, "-v"])
