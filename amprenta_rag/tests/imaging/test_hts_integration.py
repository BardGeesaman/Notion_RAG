"""Unit tests for HTS integration functionality."""

from __future__ import annotations

import uuid
from datetime import datetime
from unittest.mock import MagicMock, patch

import pytest
from sqlalchemy.orm import Session

from amprenta_rag.imaging.hts_integration import HTSImagingIntegration
from amprenta_rag.imaging.aggregation import PlateAggregator
from amprenta_rag.imaging.qc_metrics import QCMetrics
from amprenta_rag.models.chemistry import HTSCampaign, HTSPlate, HTSWell, Compound
from amprenta_rag.imaging.models import MicroscopyImage, CellSegmentation, CellFeature


class TestHTSImagingIntegration:
    """Test HTS imaging integration functionality."""

    def setup_method(self):
        """Set up test fixtures."""
        self.integration = HTSImagingIntegration()
        self.mock_db = MagicMock(spec=Session)
        
        # Sample UUIDs
        self.campaign_id = uuid.uuid4()
        self.plate_id = uuid.uuid4()
        self.well_id = uuid.uuid4()
        self.compound_id = uuid.uuid4()

    def test_link_images_to_plate(self):
        """Test linking images to plate wells."""
        # Mock existing plate
        mock_plate = MagicMock()
        mock_plate.id = self.plate_id
        mock_plate.barcode = "PLATE001"
        self.mock_db.query.return_value.filter.return_value.first.return_value = mock_plate
        
        # Mock well creation
        mock_well = MagicMock()
        mock_well.id = self.well_id
        mock_well.position = "A01"
        
        # Setup query chain for well lookup (returns None first, then returns created well)
        well_query = MagicMock()
        well_query.filter.return_value.first.return_value = None  # No existing well
        
        # Setup image query to return None (no existing image)
        image_query = MagicMock()
        image_query.filter.return_value.first.return_value = None
        
        self.mock_db.query.side_effect = [
            MagicMock(filter=MagicMock(return_value=MagicMock(first=MagicMock(return_value=mock_plate)))),  # Plate query
            well_query,  # Well lookup query
            image_query  # Image lookup query
        ]
        
        # Sample image data
        images = [
            {
                "well_position": "A01",
                "image_path": "/images/plate001_A01_DAPI.tiff",
                "channel": "DAPI",
                "z_slice": 0,
                "timepoint": 0,
                "width": 1024,
                "height": 1024,
                "bit_depth": 16,
                "pixel_size_um": 0.325
            }
        ]
        
        # Test linking
        result = self.integration.link_images_to_plate(
            "PLATE001", images, self.campaign_id, self.mock_db
        )
        
        assert result["plate_barcode"] == "PLATE001"
        assert result["images_processed"] == 1
        assert result["images_linked"] == 1
        assert result["wells_created"] == 1
        assert len(result["errors"]) == 0

    def test_link_images_creates_wells(self):
        """Test that linking images auto-creates missing wells."""
        # Mock refresh to return created objects
        def mock_refresh(obj):
            if hasattr(obj, 'barcode'):  # Plate object
                obj.id = self.plate_id
            elif hasattr(obj, 'position'):  # Well object
                obj.id = self.well_id
        
        self.mock_db.refresh.side_effect = mock_refresh
        
        # Mock queries to return None (no existing records)
        plate_query = MagicMock()
        plate_query.filter.return_value.first.return_value = None
        
        well_query = MagicMock()
        well_query.filter.return_value.first.return_value = None
        
        image_query = MagicMock()
        image_query.filter.return_value.first.return_value = None
        
        self.mock_db.query.side_effect = [
            plate_query,  # No existing plate
            well_query,   # No existing well
            image_query   # No existing image
        ]
        
        images = [
            {
                "well_position": "B12",
                "image_path": "/images/test.tiff",
                "channel": "GFP",
                "compound_id": str(self.compound_id),
                "concentration": 10.0
            }
        ]
        
        result = self.integration.link_images_to_plate(
            "NEWPLATE", images, self.campaign_id, self.mock_db, auto_create_wells=True
        )
        
        assert result["wells_created"] == 1
        assert result["images_linked"] == 1
        
        # Verify plate creation was called
        self.mock_db.add.assert_called()
        self.mock_db.commit.assert_called()

    def test_get_plate_images(self):
        """Test retrieving all images for a plate."""
        # Mock images
        mock_images = [
            MagicMock(channel="DAPI", z_slice=0),
            MagicMock(channel="GFP", z_slice=0),
            MagicMock(channel="RFP", z_slice=0)
        ]
        
        # Mock query chain
        query_mock = MagicMock()
        query_mock.join.return_value.filter.return_value.order_by.return_value.all.return_value = mock_images
        self.mock_db.query.return_value = query_mock
        
        # Test retrieval
        images = self.integration.get_plate_images(self.plate_id, self.mock_db)
        
        assert len(images) == 3
        assert images == mock_images

    def test_get_well_images(self):
        """Test retrieving images for a specific well."""
        # Mock images for a well
        mock_images = [
            MagicMock(channel="DAPI"),
            MagicMock(channel="GFP")
        ]
        
        query_mock = MagicMock()
        query_mock.filter.return_value.order_by.return_value.all.return_value = mock_images
        self.mock_db.query.return_value = query_mock
        
        images = self.integration.get_well_images(self.well_id, self.mock_db)
        
        assert len(images) == 2
        assert images == mock_images

    def test_aggregate_plate_features(self):
        """Test aggregating features across a plate."""
        aggregator = PlateAggregator()
        
        # Mock wells
        mock_wells = [
            MagicMock(id=uuid.uuid4(), position="A01"),
            MagicMock(id=uuid.uuid4(), position="A02")
        ]
        self.mock_db.query.return_value.filter.return_value.all.return_value = mock_wells
        
        # Mock features for each well
        mock_features = [
            MagicMock(
                cell_id=1, area=100.0, perimeter=40.0, circularity=0.8,
                eccentricity=0.5, solidity=0.9, centroid_x=50.0, centroid_y=60.0,
                intensity_features={}, custom_features={}
            ),
            MagicMock(
                cell_id=2, area=120.0, perimeter=45.0, circularity=0.75,
                eccentricity=0.6, solidity=0.85, centroid_x=80.0, centroid_y=90.0,
                intensity_features={}, custom_features={}
            )
        ]
        
        # Mock query for features
        feature_query = MagicMock()
        feature_query.join.return_value.join.return_value.filter.return_value.all.return_value = mock_features
        
        # Mock separate queries for plate and features
        self.mock_db.query.side_effect = [
            MagicMock(filter=MagicMock(return_value=MagicMock(all=MagicMock(return_value=mock_wells)))),
            feature_query,  # First well features
            feature_query   # Second well features
        ]
        
        df = aggregator.aggregate_plate_features(self.plate_id, db=self.mock_db)
        
        assert len(df) == 2
        assert "well_position" in df.columns
        assert "area_mean" in df.columns

    def test_aggregate_well_features(self):
        """Test aggregating features for a single well."""
        aggregator = PlateAggregator()
        
        # Mock cell features
        mock_features = [
            MagicMock(
                cell_id=1, area=100.0, perimeter=40.0, circularity=0.8,
                eccentricity=0.5, solidity=0.9, centroid_x=50.0, centroid_y=60.0,
                intensity_features={}, custom_features={}
            ),
            MagicMock(
                cell_id=2, area=120.0, perimeter=45.0, circularity=0.75,
                eccentricity=0.6, solidity=0.85, centroid_x=80.0, centroid_y=90.0,
                intensity_features={}, custom_features={}
            )
        ]
        
        query_mock = MagicMock()
        query_mock.join.return_value.join.return_value.filter.return_value.all.return_value = mock_features
        self.mock_db.query.return_value = query_mock
        
        features = aggregator.aggregate_well_features(self.well_id, db=self.mock_db)
        
        assert "area_mean" in features
        assert "total_cell_count" in features
        assert features["total_cell_count"] == 2
        assert features["area_mean"] == 110.0  # (100 + 120) / 2

    def test_plate_heatmap_data(self):
        """Test generating heatmap data for plate visualization."""
        aggregator = PlateAggregator()
        
        # Mock aggregated DataFrame
        import pandas as pd
        mock_df = pd.DataFrame({
            'well_position': ['A01', 'A02', 'B01'],
            'area_mean': [100.0, 110.0, 95.0],
            'well_id': [str(uuid.uuid4()) for _ in range(3)]
        })
        
        with patch.object(aggregator, 'aggregate_plate_features', return_value=mock_df):
            heatmap_data = aggregator.get_plate_heatmap_data(
                self.plate_id, "area", "mean", self.mock_db
            )
        
        assert len(heatmap_data) == 3
        assert heatmap_data["A01"] == 100.0
        assert heatmap_data["A02"] == 110.0
        assert heatmap_data["B01"] == 95.0

    def test_calculate_zprime(self):
        """Test Z' factor calculation for assay quality."""
        qc_metrics = QCMetrics()
        
        # Mock aggregated DataFrame with control wells
        import pandas as pd
        mock_df = pd.DataFrame({
            'well_position': ['A01', 'A02', 'H01', 'H02', 'D05', 'D06'],
            'area_mean': [50.0, 52.0, 150.0, 148.0, 100.0, 105.0],  # neg, neg, pos, pos, sample, sample
            'well_id': [str(uuid.uuid4()) for _ in range(6)]
        })
        
        with patch.object(qc_metrics.aggregator, 'aggregate_plate_features', return_value=mock_df):
            result = qc_metrics.calculate_zprime(
                self.plate_id, "area", ["H01", "H02"], ["A01", "A02"], db=self.mock_db
            )
        
        assert "zprime_factor" in result
        assert "assay_quality" in result
        assert result["positive_controls"]["mean"] == 149.0
        assert result["negative_controls"]["mean"] == 51.0

    def test_calculate_zscore(self):
        """Test Z-score calculation for outlier detection."""
        qc_metrics = QCMetrics()
        
        # Mock aggregated DataFrame
        import pandas as pd
        mock_df = pd.DataFrame({
            'well_position': ['A01', 'A02', 'A03', 'A04'],
            'area_mean': [100.0, 102.0, 98.0, 200.0],  # Last one is outlier
            'well_id': [str(uuid.uuid4()) for _ in range(4)]
        })
        
        with patch.object(qc_metrics.aggregator, 'aggregate_plate_features', return_value=mock_df):
            result = qc_metrics.calculate_zscore(
                self.plate_id, "area", db=self.mock_db
            )
        
        assert "zscore_data" in result
        assert "outliers_2sigma" in result
        assert len(result["zscore_data"]) == 4
        
        # Check that outliers are detected (A04 with value 200 should be an outlier)
        # Calculate expected Z-score: (200 - 125) / std ≈ 3.9 (should be > 2)
        # Mean of [100, 102, 98, 200] = 125, std ≈ 49.24
        assert len(result["outliers_2sigma"]) >= 1
        assert "A04" in result["outliers_2sigma"] or len(result["outliers_3sigma"]) >= 1

    def test_calculate_cv(self):
        """Test coefficient of variation calculation for a well."""
        qc_metrics = QCMetrics()
        
        # Mock cell features with varying areas
        mock_features = [
            MagicMock(area=100.0),
            MagicMock(area=110.0),
            MagicMock(area=90.0),
            MagicMock(area=105.0),
            MagicMock(area=95.0)
        ]
        
        query_mock = MagicMock()
        query_mock.join.return_value.join.return_value.filter.return_value.all.return_value = mock_features
        self.mock_db.query.return_value = query_mock
        
        result = qc_metrics.calculate_cv(self.well_id, "area", self.mock_db)
        
        assert "cv_percent" in result
        assert "mean" in result
        assert "std" in result
        assert result["cell_count"] == 5
        assert result["mean"] == 100.0  # (100+110+90+105+95)/5
        
        # CV should be reasonable for this data
        assert 0 < result["cv_percent"] < 50
