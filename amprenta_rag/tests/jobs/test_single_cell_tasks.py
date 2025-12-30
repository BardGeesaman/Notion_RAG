"""Tests for single-cell Celery tasks."""

import os
from datetime import datetime, timezone
from unittest.mock import MagicMock, patch
from uuid import uuid4

import pytest

from amprenta_rag.jobs.tasks.single_cell import process_single_cell


class TestSingleCellTasks:
    """Test single-cell Celery tasks."""
    
    def setup_method(self):
        os.environ["CELERY_TASK_ALWAYS_EAGER"] = "true"
    
    def teardown_method(self):
        os.environ.pop("CELERY_TASK_ALWAYS_EAGER", None)
    
    @patch('amprenta_rag.database.session.db_session')
    @patch('amprenta_rag.single_cell.h5ad_parser.load_h5ad')
    @patch('amprenta_rag.single_cell.h5ad_parser.validate_h5ad')
    @patch('amprenta_rag.single_cell.h5ad_parser.extract_metadata')
    @patch('amprenta_rag.single_cell.scanpy_pipeline.run_preprocessing')
    @patch('amprenta_rag.single_cell.marker_discovery.find_markers')
    @patch('amprenta_rag.single_cell.gene_mapper.map_genes_to_features')
    def test_single_cell_success(self, mock_gene_mapper, mock_find_markers, mock_preprocessing, 
                                 mock_extract_metadata, mock_validate_h5ad, mock_load_h5ad, mock_db_session):
        """Test successful single-cell processing."""
        dataset_id = uuid4()
        
        # Mock h5ad data
        mock_adata = MagicMock()
        mock_adata.n_obs = 1000
        mock_adata.n_vars = 500
        mock_load_h5ad.return_value = mock_adata
        
        # Mock metadata
        mock_extract_metadata.return_value = {
            "n_cells": 1000,
            "n_genes": 500
        }
        
        # Mock processed data
        mock_processed = MagicMock()
        mock_processed.n_obs = 1000
        mock_processed.n_vars = 500
        mock_processed.obs = MagicMock()
        mock_processed.obs.index = [f"cell_{i}" for i in range(10)]  # Sample cells
        mock_processed.obs.iloc = [MagicMock() for _ in range(10)]
        for i, cell in enumerate(mock_processed.obs.iloc):
            cell.get.return_value = i % 3  # Mock leiden cluster
        mock_processed.obs.__getitem__.return_value.unique.return_value = [0, 1, 2]  # 3 clusters
        mock_processed.obs.__getitem__.return_value.__eq__.return_value.sum.return_value = 3  # cells per cluster
        # Mock UMAP as numpy-like array that supports tuple indexing
        import numpy as np
        mock_processed.obsm = {"X_umap": np.array([[i, i+1] for i in range(10)])}  # Mock UMAP coords
        mock_preprocessing.return_value = mock_processed
        
        # Mock markers
        mock_markers = MagicMock()
        mock_markers.empty = False
        mock_markers.iterrows.return_value = [
            (0, {"gene_symbol": "ACTB", "cluster_id": 0, "log2_fold_change": 1.5, "pval": 0.01, "pval_adj": 0.05, "pct_in_cluster": 0.8, "pct_out_cluster": 0.3})
        ]
        mock_markers.get.return_value.tolist.return_value = ["ACTB"]
        mock_find_markers.return_value = mock_markers
        
        # Mock gene mapping
        mock_gene_mapper.return_value = {"ACTB": uuid4()}
        
        # Mock dataset
        mock_dataset = MagicMock()
        mock_dataset.id = dataset_id
        mock_dataset.h5ad_path = "/path/to/data.h5ad"
        mock_dataset.clustering_resolution = 0.5
        mock_dataset.n_cells = 1000
        mock_dataset.n_genes = 500
        
        # Mock database session
        mock_db = MagicMock()
        mock_db_session.return_value.__enter__.return_value = mock_db
        mock_db.query.return_value.filter_by.return_value.first.return_value = mock_dataset
        
        # Execute task
        result = process_single_cell(str(dataset_id))
        
        # Verify result
        assert result["status"] == "completed"
        assert result["dataset_id"] == str(dataset_id)
        assert result["n_cells"] == 1000
        assert result["n_genes"] == 500
        assert result["n_clusters"] == 3
        
        # Verify pipeline was called
        mock_load_h5ad.assert_called_once_with("/path/to/data.h5ad")
        mock_validate_h5ad.assert_called_once_with(mock_adata)
        mock_preprocessing.assert_called_once_with(mock_adata, resolution=0.5)
        mock_find_markers.assert_called_once_with(mock_processed, groupby="leiden", top_n=50)
        
        # Verify dataset status was updated
        assert mock_dataset.processing_status == "completed"
        assert mock_dataset.processed_at is not None
    
    @patch('amprenta_rag.database.session.db_session')
    def test_single_cell_dataset_not_found(self, mock_db_session):
        """Test single-cell task when dataset doesn't exist."""
        dataset_id = uuid4()
        
        # Mock database session with no dataset found
        mock_db = MagicMock()
        mock_db_session.return_value.__enter__.return_value = mock_db
        mock_db.query.return_value.filter_by.return_value.first.return_value = None
        
        # Execute task
        result = process_single_cell(str(dataset_id))
        
        # Verify result
        assert result["status"] == "failed"
        assert result["error"] == "Dataset not found"
        assert result["dataset_id"] == str(dataset_id)
    
    @patch('amprenta_rag.database.session.db_session')
    @patch('amprenta_rag.single_cell.h5ad_parser.load_h5ad')
    def test_single_cell_h5ad_load_error(self, mock_load_h5ad, mock_db_session):
        """Test single-cell task when h5ad file loading fails."""
        dataset_id = uuid4()
        
        # Mock dataset
        mock_dataset = MagicMock()
        mock_dataset.id = dataset_id
        mock_dataset.h5ad_path = "/path/to/missing.h5ad"
        
        # Mock database session
        mock_db = MagicMock()
        mock_db_session.return_value.__enter__.return_value = mock_db
        mock_db.query.return_value.filter_by.return_value.first.return_value = mock_dataset
        
        # Mock h5ad loading to raise exception
        mock_load_h5ad.side_effect = FileNotFoundError("H5AD file not found")
        
        # Execute task - should raise exception for retry
        with pytest.raises(Exception) as exc_info:
            process_single_cell(str(dataset_id))
        
        # Verify correct error
        assert "H5AD file not found" in str(exc_info.value)
        
        # Verify dataset status was updated to failed
        assert mock_dataset.processing_status == "failed"
        assert "H5AD file not found" in mock_dataset.processing_log
    
    @patch('amprenta_rag.database.session.db_session')
    @patch('amprenta_rag.single_cell.h5ad_parser.load_h5ad')
    @patch('amprenta_rag.single_cell.h5ad_parser.validate_h5ad')
    def test_single_cell_validation_error(self, mock_validate_h5ad, mock_load_h5ad, mock_db_session):
        """Test single-cell task when h5ad validation fails."""
        dataset_id = uuid4()
        
        # Mock h5ad data
        mock_adata = MagicMock()
        mock_load_h5ad.return_value = mock_adata
        
        # Mock validation to raise exception
        mock_validate_h5ad.side_effect = ValueError("Invalid h5ad format")
        
        # Mock dataset
        mock_dataset = MagicMock()
        mock_dataset.id = dataset_id
        mock_dataset.h5ad_path = "/path/to/invalid.h5ad"
        
        # Mock database session
        mock_db = MagicMock()
        mock_db_session.return_value.__enter__.return_value = mock_db
        mock_db.query.return_value.filter_by.return_value.first.return_value = mock_dataset
        
        # Execute task - should raise exception for retry
        with pytest.raises(Exception) as exc_info:
            process_single_cell(str(dataset_id))
        
        # Verify correct error
        assert "Invalid h5ad format" in str(exc_info.value)
    
    @patch('amprenta_rag.database.session.db_session')
    @patch('amprenta_rag.single_cell.h5ad_parser.load_h5ad')
    @patch('amprenta_rag.single_cell.h5ad_parser.validate_h5ad')
    @patch('amprenta_rag.single_cell.h5ad_parser.extract_metadata')
    @patch('amprenta_rag.single_cell.scanpy_pipeline.run_preprocessing')
    def test_single_cell_preprocessing_error(self, mock_preprocessing, mock_extract_metadata, 
                                           mock_validate_h5ad, mock_load_h5ad, mock_db_session):
        """Test single-cell task when scanpy pipeline fails."""
        dataset_id = uuid4()
        
        # Mock h5ad data
        mock_adata = MagicMock()
        mock_load_h5ad.return_value = mock_adata
        mock_extract_metadata.return_value = {"n_cells": 1000, "n_genes": 500}
        
        # Mock preprocessing to raise exception
        mock_preprocessing.side_effect = RuntimeError("Scanpy pipeline failed")
        
        # Mock dataset
        mock_dataset = MagicMock()
        mock_dataset.id = dataset_id
        mock_dataset.h5ad_path = "/path/to/data.h5ad"
        mock_dataset.clustering_resolution = 0.5
        
        # Mock database session
        mock_db = MagicMock()
        mock_db_session.return_value.__enter__.return_value = mock_db
        mock_db.query.return_value.filter_by.return_value.first.return_value = mock_dataset
        
        # Execute task - should raise exception for retry
        with pytest.raises(Exception) as exc_info:
            process_single_cell(str(dataset_id))
        
        # Verify correct error
        assert "Scanpy pipeline failed" in str(exc_info.value)
    
    @patch('amprenta_rag.database.session.db_session')
    @patch('amprenta_rag.single_cell.h5ad_parser.load_h5ad')
    @patch('amprenta_rag.single_cell.h5ad_parser.validate_h5ad')
    @patch('amprenta_rag.single_cell.h5ad_parser.extract_metadata')
    @patch('amprenta_rag.single_cell.scanpy_pipeline.run_preprocessing')
    @patch('amprenta_rag.single_cell.marker_discovery.find_markers')
    def test_single_cell_marker_discovery(self, mock_find_markers, mock_preprocessing, 
                                         mock_extract_metadata, mock_validate_h5ad, mock_load_h5ad, mock_db_session):
        """Test single-cell task validates marker discovery integration."""
        dataset_id = uuid4()
        
        # Mock h5ad data
        mock_adata = MagicMock()
        mock_load_h5ad.return_value = mock_adata
        mock_extract_metadata.return_value = {"n_cells": 100, "n_genes": 200}
        
        # Mock processed data
        mock_processed = MagicMock()
        mock_processed.n_obs = 100
        mock_processed.n_vars = 200
        mock_processed.obs = MagicMock()
        mock_processed.obs.index = ["cell_1", "cell_2"]
        mock_processed.obs.iloc = [MagicMock(), MagicMock()]
        mock_processed.obs.iloc[0].get.return_value = 0
        mock_processed.obs.iloc[1].get.return_value = 1
        mock_processed.obs.__getitem__.return_value.unique.return_value = [0, 1]
        mock_processed.obs.__getitem__.return_value.__eq__.return_value.sum.return_value = 1
        import numpy as np
        mock_processed.obsm = {"X_umap": np.array([[0, 1], [2, 3]])}
        mock_preprocessing.return_value = mock_processed
        
        # Mock markers with specific genes
        mock_markers = MagicMock()
        mock_markers.empty = False
        mock_markers.iterrows.return_value = [
            (0, {"gene_symbol": "CD3E", "cluster_id": 0, "log2_fold_change": 2.0, "pval": 0.001, "pval_adj": 0.01, "pct_in_cluster": 0.9, "pct_out_cluster": 0.1}),
            (1, {"gene_symbol": "CD19", "cluster_id": 1, "log2_fold_change": 1.8, "pval": 0.002, "pval_adj": 0.02, "pct_in_cluster": 0.85, "pct_out_cluster": 0.15})
        ]
        mock_markers.get.return_value.tolist.return_value = ["CD3E", "CD19"]
        mock_find_markers.return_value = mock_markers
        
        # Mock dataset
        mock_dataset = MagicMock()
        mock_dataset.id = dataset_id
        mock_dataset.h5ad_path = "/path/to/data.h5ad"
        mock_dataset.clustering_resolution = 0.5
        
        # Mock database session
        mock_db = MagicMock()
        mock_db_session.return_value.__enter__.return_value = mock_db
        mock_db.query.return_value.filter_by.return_value.first.return_value = mock_dataset
        
        # Mock gene mapping
        with patch('amprenta_rag.single_cell.gene_mapper.map_genes_to_features') as mock_gene_mapper:
            mock_gene_mapper.return_value = {"CD3E": uuid4(), "CD19": uuid4()}
            
            # Execute task
            result = process_single_cell(str(dataset_id))
        
        # Verify marker discovery was called correctly
        mock_find_markers.assert_called_once_with(mock_processed, groupby="leiden", top_n=50)
        
        # Verify result includes marker data
        assert result["status"] == "completed"
        assert result["n_clusters"] == 2
    
    @patch('amprenta_rag.database.session.db_session')
    @patch('amprenta_rag.single_cell.h5ad_parser.load_h5ad')
    @patch('amprenta_rag.single_cell.h5ad_parser.validate_h5ad')
    @patch('amprenta_rag.single_cell.h5ad_parser.extract_metadata')
    @patch('amprenta_rag.single_cell.scanpy_pipeline.run_preprocessing')
    @patch('amprenta_rag.single_cell.marker_discovery.find_markers')
    @patch('amprenta_rag.single_cell.gene_mapper.map_genes_to_features')
    def test_single_cell_status_transitions(self, mock_gene_mapper, mock_find_markers, mock_preprocessing,
                                           mock_extract_metadata, mock_validate_h5ad, mock_load_h5ad, mock_db_session):
        """Test single-cell task status transitions."""
        dataset_id = uuid4()
        
        # Mock h5ad data
        mock_adata = MagicMock()
        mock_load_h5ad.return_value = mock_adata
        mock_extract_metadata.return_value = {"n_cells": 50, "n_genes": 100}
        
        # Mock processed data (minimal)
        mock_processed = MagicMock()
        mock_processed.n_obs = 50
        mock_processed.n_vars = 100
        mock_processed.obs = MagicMock()
        mock_processed.obs.index = ["cell_1"]
        mock_processed.obs.iloc = [MagicMock()]
        mock_processed.obs.iloc[0].get.return_value = 0
        mock_processed.obs.__getitem__.return_value.unique.return_value = [0]
        mock_processed.obs.__getitem__.return_value.__eq__.return_value.sum.return_value = 1
        import numpy as np
        mock_processed.obsm = {"X_umap": np.array([[0, 1]])}
        mock_preprocessing.return_value = mock_processed
        
        # Mock empty markers
        mock_markers = MagicMock()
        mock_markers.empty = True
        mock_markers.iterrows.return_value = []
        mock_find_markers.return_value = mock_markers
        mock_gene_mapper.return_value = {}
        
        # Mock dataset
        mock_dataset = MagicMock()
        mock_dataset.id = dataset_id
        mock_dataset.h5ad_path = "/path/to/data.h5ad"
        mock_dataset.clustering_resolution = None  # Test default resolution
        
        # Mock database session
        mock_db = MagicMock()
        mock_db_session.return_value.__enter__.return_value = mock_db
        mock_db.query.return_value.filter_by.return_value.first.return_value = mock_dataset
        
        # Execute task
        result = process_single_cell(str(dataset_id))
        
        # Verify status transitions
        # First call sets to "running", final call sets to "completed"
        assert mock_dataset.processing_status == "completed"
        assert mock_dataset.processed_at is not None
        assert result["status"] == "completed"
    
    @patch('amprenta_rag.database.session.db_session')
    @patch('amprenta_rag.single_cell.h5ad_parser.load_h5ad')
    def test_single_cell_retry_on_failure(self, mock_load_h5ad, mock_db_session):
        """Test single-cell task handles failure and re-raises for retry."""
        dataset_id = uuid4()
        
        # Mock dataset
        mock_dataset = MagicMock()
        mock_dataset.id = dataset_id
        mock_dataset.h5ad_path = "/path/to/data.h5ad"
        
        # Mock database session
        mock_db = MagicMock()
        mock_db_session.return_value.__enter__.return_value = mock_db
        mock_db.query.return_value.filter_by.return_value.first.return_value = mock_dataset
        
        # Mock h5ad loading to raise exception
        mock_load_h5ad.side_effect = Exception("Processing error")
        
        # Execute task - should raise exception for retry
        with pytest.raises(Exception) as exc_info:
            process_single_cell(str(dataset_id))
        
        # Verify correct error
        assert "Processing error" in str(exc_info.value)
    
    @patch('amprenta_rag.database.session.db_session')
    @patch('amprenta_rag.single_cell.h5ad_parser.load_h5ad')
    def test_single_cell_max_retries_exceeded(self, mock_load_h5ad, mock_db_session):
        """Test single-cell task returns failure dict when max retries exceeded."""
        dataset_id = uuid4()
        
        # Mock dataset
        mock_dataset = MagicMock()
        mock_dataset.id = dataset_id
        mock_dataset.h5ad_path = "/path/to/data.h5ad"
        
        # Mock database session
        mock_db = MagicMock()
        mock_db_session.return_value.__enter__.return_value = mock_db
        mock_db.query.return_value.filter_by.return_value.first.return_value = mock_dataset
        
        # Mock h5ad loading to raise exception
        mock_load_h5ad.side_effect = Exception("Persistent failure")
        
        # Execute task - should raise exception for retry
        with pytest.raises(Exception) as exc_info:
            process_single_cell(str(dataset_id))
        
        # Verify correct error
        assert "Persistent failure" in str(exc_info.value)
    
    @patch('amprenta_rag.database.session.db_session')
    @patch('amprenta_rag.single_cell.h5ad_parser.load_h5ad')
    @patch('amprenta_rag.single_cell.h5ad_parser.validate_h5ad')
    @patch('amprenta_rag.single_cell.h5ad_parser.extract_metadata')
    @patch('amprenta_rag.single_cell.scanpy_pipeline.run_preprocessing')
    @patch('amprenta_rag.single_cell.marker_discovery.find_markers')
    @patch('amprenta_rag.single_cell.gene_mapper.map_genes_to_features')
    def test_single_cell_cluster_creation(self, mock_gene_mapper, mock_find_markers, mock_preprocessing,
                                         mock_extract_metadata, mock_validate_h5ad, mock_load_h5ad, mock_db_session):
        """Test single-cell task validates CellCluster insertion."""
        dataset_id = uuid4()
        
        # Mock h5ad data
        mock_adata = MagicMock()
        mock_load_h5ad.return_value = mock_adata
        mock_extract_metadata.return_value = {"n_cells": 300, "n_genes": 150}
        
        # Mock processed data with 3 clusters
        mock_processed = MagicMock()
        mock_processed.n_obs = 300
        mock_processed.n_vars = 150
        mock_processed.obs = MagicMock()
        mock_processed.obs.index = [f"cell_{i}" for i in range(3)]
        mock_processed.obs.iloc = [MagicMock() for _ in range(3)]
        for i, cell in enumerate(mock_processed.obs.iloc):
            cell.get.return_value = i  # Each cell in different cluster
        mock_processed.obs.__getitem__.return_value.unique.return_value = [0, 1, 2]
        mock_processed.obs.__getitem__.return_value.__eq__.return_value.sum.return_value = 100  # 100 cells per cluster
        import numpy as np
        mock_processed.obsm = {"X_umap": np.array([[i, i+1] for i in range(3)])}
        mock_preprocessing.return_value = mock_processed
        
        # Mock empty markers
        mock_markers = MagicMock()
        mock_markers.empty = True
        mock_markers.iterrows.return_value = []
        mock_find_markers.return_value = mock_markers
        mock_gene_mapper.return_value = {}
        
        # Mock dataset
        mock_dataset = MagicMock()
        mock_dataset.id = dataset_id
        mock_dataset.h5ad_path = "/path/to/data.h5ad"
        mock_dataset.clustering_resolution = 0.5
        
        # Mock database session
        mock_db = MagicMock()
        mock_db_session.return_value.__enter__.return_value = mock_db
        mock_db.query.return_value.filter_by.return_value.first.return_value = mock_dataset
        
        # Execute task
        result = process_single_cell(str(dataset_id))
        
        # Verify clusters were created
        assert result["status"] == "completed"
        assert result["n_clusters"] == 3
        
        # Verify CellCluster objects were added (3 clusters)
        # Note: mock_db.add should be called for each cluster + annotations + dataset
        assert mock_db.add.call_count >= 3  # At least 3 cluster records
    
    @patch('amprenta_rag.database.session.db_session')
    @patch('amprenta_rag.single_cell.h5ad_parser.load_h5ad')
    @patch('amprenta_rag.single_cell.h5ad_parser.validate_h5ad')
    @patch('amprenta_rag.single_cell.h5ad_parser.extract_metadata')
    @patch('amprenta_rag.single_cell.scanpy_pipeline.run_preprocessing')
    @patch('amprenta_rag.single_cell.marker_discovery.find_markers')
    @patch('amprenta_rag.single_cell.gene_mapper.map_genes_to_features')
    def test_single_cell_annotation_creation(self, mock_gene_mapper, mock_find_markers, mock_preprocessing,
                                            mock_extract_metadata, mock_validate_h5ad, mock_load_h5ad, mock_db_session):
        """Test single-cell task validates CellAnnotation insertion."""
        dataset_id = uuid4()
        
        # Mock h5ad data
        mock_adata = MagicMock()
        mock_load_h5ad.return_value = mock_adata
        mock_extract_metadata.return_value = {"n_cells": 5, "n_genes": 100}
        
        # Mock processed data with specific cell annotations
        mock_processed = MagicMock()
        mock_processed.n_obs = 5
        mock_processed.n_vars = 100
        mock_processed.obs = MagicMock()
        mock_processed.obs.index = ["AAACCTGAGAAGGCCT", "AAACCTGAGACAGACC", "AAACCTGAGATCTGAA", "AAACCTGAGCACGACG", "AAACCTGAGCAGCGTA"]
        mock_processed.obs.iloc = [MagicMock() for _ in range(5)]
        for i, cell in enumerate(mock_processed.obs.iloc):
            cell.get.return_value = i % 2  # Alternate between clusters 0 and 1
        mock_processed.obs.__getitem__.return_value.unique.return_value = [0, 1]
        mock_processed.obs.__getitem__.return_value.__eq__.return_value.sum.return_value = 2
        import numpy as np
        mock_processed.obsm = {"X_umap": np.array([[i*0.5, i*0.3] for i in range(5)])}  # Specific UMAP coordinates
        mock_preprocessing.return_value = mock_processed
        
        # Mock empty markers
        mock_markers = MagicMock()
        mock_markers.empty = True
        mock_markers.iterrows.return_value = []
        mock_find_markers.return_value = mock_markers
        mock_gene_mapper.return_value = {}
        
        # Mock dataset
        mock_dataset = MagicMock()
        mock_dataset.id = dataset_id
        mock_dataset.h5ad_path = "/path/to/data.h5ad"
        mock_dataset.clustering_resolution = 0.8
        
        # Mock database session
        mock_db = MagicMock()
        mock_db_session.return_value.__enter__.return_value = mock_db
        mock_db.query.return_value.filter_by.return_value.first.return_value = mock_dataset
        
        # Execute task
        result = process_single_cell(str(dataset_id))
        
        # Verify annotations were created
        assert result["status"] == "completed"
        assert result["n_cells"] == 5
        
        # Verify CellAnnotation objects were added (5 cells + 2 clusters + dataset)
        assert mock_db.add.call_count >= 5  # At least 5 annotation records
    
    @patch('amprenta_rag.database.session.db_session')
    @patch('amprenta_rag.single_cell.h5ad_parser.load_h5ad')
    def test_single_cell_db_update_on_failure(self, mock_load_h5ad, mock_db_session):
        """Test single-cell task updates database processing_status on failure."""
        dataset_id = uuid4()
        
        # Mock dataset
        mock_dataset = MagicMock()
        mock_dataset.id = dataset_id
        mock_dataset.h5ad_path = "/path/to/data.h5ad"
        mock_dataset.processing_status = "pending"
        mock_dataset.processing_log = None
        mock_dataset.processed_at = None
        
        # Mock database session
        mock_db = MagicMock()
        mock_db_session.return_value.__enter__.return_value = mock_db
        mock_db.query.return_value.filter_by.return_value.first.return_value = mock_dataset
        
        # Mock h5ad loading to raise exception
        mock_load_h5ad.side_effect = Exception("File processing error")
        
        # Execute task - should raise exception for retry
        with pytest.raises(Exception) as exc_info:
            process_single_cell(str(dataset_id))
        
        # Verify correct error
        assert "File processing error" in str(exc_info.value)
        
        # Verify dataset status was updated to failed
        assert mock_dataset.processing_status == "failed"
        assert "File processing error" in mock_dataset.processing_log
        assert mock_dataset.processed_at is not None
        mock_db.add.assert_called_with(mock_dataset)
        mock_db.commit.assert_called()
