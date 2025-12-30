"""Unit tests for CellPose service and image analysis pipeline."""

from __future__ import annotations

import numpy as np
import pytest
from unittest.mock import MagicMock, patch, Mock
from typing import Tuple

from amprenta_rag.imaging.cellpose_service import CellPoseService
from amprenta_rag.imaging.tiling import TileManager, create_test_image, create_test_mask
from amprenta_rag.imaging.feature_extraction import (
    FeatureExtractor, 
    CellMorphologyFeatures, 
    CellIntensityFeatures,
    WellAggregatedFeatures
)


class TestCellPoseService:
    """Test CellPose service functionality."""

    def test_cellpose_service_init(self):
        """Test initializing CellPose service with different model types."""
        # Test default initialization
        service = CellPoseService()
        assert service.model_type == "cyto"
        assert service.use_gpu is True
        assert service.tile_manager.tile_size == 1024
        assert service.tile_manager.overlap == 128
        
        # Test custom initialization
        service_custom = CellPoseService(
            model_type="nuclei", 
            gpu=False, 
            tile_size=512, 
            overlap=64
        )
        assert service_custom.model_type == "nuclei"
        assert service_custom.use_gpu is False
        assert service_custom.tile_manager.tile_size == 512
        assert service_custom.tile_manager.overlap == 64

    @patch('amprenta_rag.imaging.cellpose_service.models')
    def test_segment_small_image(self, mock_models):
        """Test segmenting small image that doesn't require tiling."""
        # Mock CellPose model
        mock_model = MagicMock()
        mock_model.eval.return_value = (
            np.array([[0, 1, 1], [0, 1, 1], [0, 0, 0]]),  # masks
            np.zeros((3, 3, 2)),  # flows
            None,  # styles
            [30]   # diams
        )
        mock_models.Cellpose.return_value = mock_model
        
        service = CellPoseService()
        service._model = mock_model  # Set mock model directly
        
        # Create small test image
        image = np.random.randint(0, 255, (100, 100), dtype=np.uint8)
        
        # Segment image
        masks, flows = service.segment(image, diameter=30)
        
        # Verify model was called
        mock_model.eval.assert_called_once()
        assert masks.shape == (3, 3)
        assert flows.shape == (3, 3, 2)

    @patch('amprenta_rag.imaging.cellpose_service.models')
    def test_segment_returns_masks(self, mock_models):
        """Test that segmentation returns properly formatted masks."""
        # Mock model with realistic output
        mock_model = MagicMock()
        test_masks = np.array([
            [0, 0, 1, 1],
            [0, 2, 1, 1], 
            [2, 2, 0, 3],
            [2, 2, 3, 3]
        ], dtype=np.int32)
        
        mock_model.eval.return_value = (
            test_masks,
            np.zeros((4, 4, 2)),
            None,
            [25]
        )
        mock_models.Cellpose.return_value = mock_model
        
        service = CellPoseService()
        service._model = mock_model
        
        image = np.random.randint(0, 255, (64, 64), dtype=np.uint8)
        masks, flows = service.segment(image)
        
        # Check mask properties
        unique_labels = np.unique(masks)
        assert 0 in unique_labels  # Background
        assert len(unique_labels) > 1  # At least one cell
        assert masks.dtype == np.int32

    @patch('amprenta_rag.imaging.cellpose_service.models')
    def test_segment_cell_count(self, mock_models):
        """Test counting cells in segmentation results."""
        # Create mock masks with known cell count
        mock_masks = np.array([
            [0, 1, 1, 0],
            [0, 1, 1, 0],
            [0, 0, 0, 2],
            [3, 3, 2, 2]
        ], dtype=np.int32)
        
        mock_model = MagicMock()
        mock_model.eval.return_value = (mock_masks, np.zeros((4, 4, 2)), None, [20])
        mock_models.Cellpose.return_value = mock_model
        
        service = CellPoseService()
        service._model = mock_model
        
        image = np.random.randint(0, 255, (50, 50), dtype=np.uint8)
        masks, _ = service.segment(image)
        
        # Count cells
        cell_count = service.count_cells(masks)
        assert cell_count == 3  # Labels 1, 2, 3 (excluding background 0)

    def test_tiling_split_image(self):
        """Test splitting large image into tiles."""
        tile_manager = TileManager(tile_size=100, overlap=20)
        
        # Create test image
        image = create_test_image(250, 300, num_cells=5)
        
        # Split into tiles
        tiles = tile_manager.split_image(image)
        
        # Verify tiling
        assert len(tiles) > 1  # Should create multiple tiles
        
        # Check tile properties
        for tile in tiles:
            assert tile.image.shape[0] <= 100  # Height <= tile_size
            assert tile.image.shape[1] <= 100  # Width <= tile_size
            assert tile.x_end > tile.x_start
            assert tile.y_end > tile.y_start

    def test_tiling_stitch_masks(self):
        """Test stitching segmentation masks from tiles."""
        tile_manager = TileManager(tile_size=100, overlap=20)
        
        # Create original mask
        original_mask = create_test_mask(150, 200, num_cells=8)
        
        # Split into tiles
        tiles = tile_manager.split_image(original_mask)
        
        # Create mock tile masks (using original data for simplicity)
        mask_tiles = []
        for tile in tiles:
            mask_tiles.append((tile.image, tile))
        
        # Stitch back together
        stitched_mask = tile_manager.stitch_masks(mask_tiles, (150, 200))
        
        # Verify stitched result
        assert stitched_mask.shape == (150, 200)
        assert stitched_mask.dtype == np.int32
        
        # Should have reasonable number of cells
        unique_labels = np.unique(stitched_mask)
        cell_count = len(unique_labels[unique_labels > 0])
        assert cell_count > 0

    @patch('amprenta_rag.imaging.cellpose_service.models')
    def test_segment_large_image_tiled(self, mock_models):
        """Test segmenting large image using tiling approach."""
        # Mock model
        mock_model = MagicMock()
        
        def mock_eval(image, **kwargs):
            # Return simple mask for each tile
            h, w = image.shape[:2]
            mask = np.zeros((h, w), dtype=np.int32)
            # Add a single cell in center
            center_y, center_x = h // 2, w // 2
            mask[center_y-5:center_y+5, center_x-5:center_x+5] = 1
            flows = np.zeros((h, w, 2), dtype=np.float32)
            return mask, flows, None, [25]
        
        mock_model.eval.side_effect = mock_eval
        mock_models.Cellpose.return_value = mock_model
        
        service = CellPoseService(tile_size=512, overlap=64)
        service._model = mock_model
        
        # Create large test image
        large_image = create_test_image(2500, 3000, num_cells=20)
        
        # Segment using tiling
        masks, flows = service.segment(large_image)
        
        # Verify results
        assert masks.shape == (2500, 3000)
        assert flows.shape == (2500, 3000, 2)
        
        # Should have detected some cells
        cell_count = service.count_cells(masks)
        assert cell_count > 0

    @patch('amprenta_rag.imaging.cellpose_service.models')
    def test_gpu_fallback_on_oom(self, mock_models):
        """Test CPU fallback when GPU runs out of memory."""
        # Mock GPU model that raises OOM error
        mock_gpu_model = MagicMock()
        mock_gpu_model.eval.side_effect = RuntimeError("CUDA out of memory")
        
        # Mock CPU model that works
        mock_cpu_model = MagicMock()
        mock_cpu_model.eval.return_value = (
            np.array([[0, 1], [1, 1]]),
            np.zeros((2, 2, 2)),
            None,
            [20]
        )
        
        # Configure mock to return different models based on gpu parameter
        def create_model(model_type, gpu):
            if gpu:
                return mock_gpu_model
            else:
                return mock_cpu_model
        
        mock_models.Cellpose.side_effect = create_model
        
        service = CellPoseService(gpu=True)  # Start with GPU enabled
        service._model = mock_gpu_model
        
        image = np.random.randint(0, 255, (100, 100), dtype=np.uint8)
        
        # Should automatically fall back to CPU
        masks, flows = service.segment(image)
        
        # Verify fallback worked
        assert masks.shape == (2, 2)
        mock_cpu_model.eval.assert_called_once()

    def test_estimate_memory(self):
        """Test GPU memory estimation function."""
        service = CellPoseService()
        
        # Test different image sizes
        small_image_shape = (512, 512)
        large_image_shape = (2048, 2048, 3)
        
        small_memory = service.estimate_memory_mb(small_image_shape)
        large_memory = service.estimate_memory_mb(large_image_shape)
        
        # Larger image should require more memory
        assert large_memory > small_memory
        
        # Memory should be reasonable (not negative, not extremely large)
        assert small_memory > 0
        assert small_memory < 10000  # Less than 10GB
        assert large_memory > 0
        assert large_memory < 50000  # Less than 50GB


class TestFeatureExtraction:
    """Test feature extraction functionality."""

    def test_feature_extraction_morphology(self):
        """Test extracting morphological features from masks."""
        extractor = FeatureExtractor()
        
        # Create test mask with known cells
        mask = np.array([
            [0, 0, 0, 0, 0],
            [0, 1, 1, 0, 0],
            [0, 1, 1, 0, 0],
            [0, 0, 0, 2, 2],
            [0, 0, 0, 2, 2]
        ], dtype=np.int32)
        
        # Mock regionprops to avoid scikit-image dependency
        mock_region1 = MagicMock()
        mock_region1.area = 4.0
        mock_region1.perimeter = 8.0
        mock_region1.major_axis_length = 2.0
        mock_region1.minor_axis_length = 2.0
        mock_region1.eccentricity = 0.0
        mock_region1.solidity = 1.0
        mock_region1.extent = 1.0
        mock_region1.centroid = (1.5, 1.5)
        
        mock_region2 = MagicMock()
        mock_region2.area = 4.0
        mock_region2.perimeter = 8.0
        mock_region2.major_axis_length = 2.0
        mock_region2.minor_axis_length = 2.0
        mock_region2.eccentricity = 0.0
        mock_region2.solidity = 1.0
        mock_region2.extent = 1.0
        mock_region2.centroid = (3.5, 3.5)
        
        with patch.object(extractor, 'regionprops', return_value=[mock_region1, mock_region2]):
            features = extractor.extract_morphology_features(mask)
        
        # Verify features
        assert len(features) == 2
        
        for feature in features:
            assert isinstance(feature, CellMorphologyFeatures)
            assert feature.area == 4.0
            assert feature.perimeter == 8.0
            assert feature.circularity > 0  # Should be calculated
            assert feature.aspect_ratio == 1.0  # Square cells

    def test_feature_extraction_intensity(self):
        """Test extracting intensity features from multi-channel images."""
        extractor = FeatureExtractor()
        
        # Create test mask
        mask = np.array([
            [0, 1, 1],
            [0, 1, 1],
            [2, 2, 0]
        ], dtype=np.int32)
        
        # Create test images for different channels
        dapi_image = np.array([
            [10, 100, 120],
            [15, 110, 115],
            [200, 180, 5]
        ], dtype=np.uint8)
        
        gfp_image = np.array([
            [5, 50, 60],
            [8, 55, 58],
            [150, 140, 2]
        ], dtype=np.uint8)
        
        images = {"DAPI": dapi_image, "GFP": gfp_image}
        
        features = extractor.extract_intensity_features(mask, images)
        
        # Verify features
        assert len(features) == 2  # Two cells
        
        for feature in features:
            assert isinstance(feature, CellIntensityFeatures)
            assert "DAPI" in feature.channel_features
            assert "GFP" in feature.channel_features
            
            # Check that statistics are calculated
            dapi_stats = feature.channel_features["DAPI"]
            assert "mean" in dapi_stats
            assert "median" in dapi_stats
            assert "std" in dapi_stats
            assert "min" in dapi_stats
            assert "max" in dapi_stats

    def test_aggregate_to_well(self):
        """Test aggregating cell features to well level."""
        extractor = FeatureExtractor()
        
        # Create mock morphology features
        morph_features = [
            CellMorphologyFeatures(
                area=100.0, perimeter=40.0, major_axis_length=12.0,
                minor_axis_length=10.0, eccentricity=0.5, solidity=0.9,
                extent=0.8, centroid_x=50.0, centroid_y=60.0,
                circularity=0.7, aspect_ratio=1.2
            ),
            CellMorphologyFeatures(
                area=120.0, perimeter=45.0, major_axis_length=14.0,
                minor_axis_length=11.0, eccentricity=0.6, solidity=0.85,
                extent=0.75, centroid_x=80.0, centroid_y=90.0,
                circularity=0.75, aspect_ratio=1.3
            )
        ]
        
        # Create mock intensity features
        intensity_features = [
            CellIntensityFeatures(channel_features={
                "DAPI": {"mean": 100.0, "std": 10.0},
                "GFP": {"mean": 50.0, "std": 5.0}
            }),
            CellIntensityFeatures(channel_features={
                "DAPI": {"mean": 120.0, "std": 12.0},
                "GFP": {"mean": 60.0, "std": 6.0}
            })
        ]
        
        # Aggregate features
        well_features = extractor.aggregate_to_well(morph_features, intensity_features)
        
        # Verify aggregation
        assert isinstance(well_features, WellAggregatedFeatures)
        assert well_features.cell_count == 2
        
        # Check morphology aggregation
        assert "area" in well_features.morphology_stats
        area_stats = well_features.morphology_stats["area"]
        assert area_stats["mean"] == 110.0  # (100 + 120) / 2
        assert area_stats["min"] == 100.0
        assert area_stats["max"] == 120.0
        
        # Check intensity aggregation
        assert "DAPI" in well_features.intensity_stats
        dapi_stats = well_features.intensity_stats["DAPI"]["mean"]
        assert dapi_stats["mean"] == 110.0  # (100 + 120) / 2


class TestTilingUtilities:
    """Test tiling utility functions."""

    def test_tile_manager_positions(self):
        """Test tile position calculation."""
        tile_manager = TileManager(tile_size=100, overlap=20)
        
        # Test position calculation for dimension of 250
        positions = tile_manager._calculate_positions(250, 100, 20)
        
        # Should create overlapping tiles
        assert len(positions) >= 3  # At least 3 tiles needed
        
        # Check that positions cover the full dimension
        assert positions[0][0] == 0  # First tile starts at 0
        assert positions[-1][1] >= 250  # Last tile covers end
        
        # Check overlap
        if len(positions) > 1:
            # Second tile should start before first tile ends
            assert positions[1][0] < positions[0][1]

    def test_create_test_utilities(self):
        """Test test image and mask creation utilities."""
        # Test image creation
        image = create_test_image(100, 150, num_cells=5)
        assert image.shape == (100, 150)
        assert image.dtype == np.uint8
        assert np.max(image) > 0  # Should have some non-zero pixels
        
        # Test mask creation
        mask = create_test_mask(100, 150, num_cells=5)
        assert mask.shape == (100, 150)
        assert mask.dtype == np.int32
        
        # Check that mask has labeled cells
        unique_labels = np.unique(mask)
        cell_count = len(unique_labels[unique_labels > 0])
        assert cell_count > 0  # Should have at least some cells
