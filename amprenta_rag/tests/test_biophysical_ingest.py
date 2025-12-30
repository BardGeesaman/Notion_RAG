"""
Unit tests for biophysical ingest service.

Tests orchestration functions for SPR, MST, and DSC file ingestion,
background processing, and entity linking.
"""

import tempfile
import threading
import time
from pathlib import Path
from unittest.mock import MagicMock, patch, call
from uuid import UUID, uuid4

import numpy as np
import pytest

from amprenta_rag.biophysical.ingest_service import (
    ingest_spr_file, ingest_mst_file, ingest_dsc_file,
    link_to_compound, link_to_target, get_compound_biophysical_profile,
    get_processing_status, reprocess_experiment,
    _process_spr_async, _process_mst_async, _process_dsc_async
)
from amprenta_rag.database.models_biophysical import (
    SPRExperiment, MSTExperiment, DSCExperiment,
    SPRExperimentType
)


class TestBiophysicalIngest:
    """Test biophysical file ingestion functions."""
    
    def create_mock_spr_file(self) -> str:
        """Create a mock SPR CSV file for testing."""
        content = """Instrument: Biacore T200
Chip: CM5
Temperature: 25.0
Flow Rate: 30
Buffer: PBS + 0.05% Tween-20
Ligand: EGFR-His6
Analyte: Compound-001

Time	Response	Concentration	Cycle
0.0	0.0	100.0	1
1.0	10.0	100.0	1
2.0	25.0	100.0	1
3.0	35.0	100.0	1
4.0	30.0	100.0	1
5.0	20.0	100.0	1
6.0	15.0	100.0	1
7.0	10.0	100.0	1
8.0	8.0	100.0	1
9.0	5.0	100.0	1
10.0	2.0	100.0	1"""
        
        with tempfile.NamedTemporaryFile(mode='w', suffix='.csv', delete=False) as f:
            f.write(content)
            return f.name
    
    def create_mock_mst_file(self) -> str:
        """Create a mock MST CSV file for testing."""
        content = """Concentration,Fnorm,Error,Cold,Hot
1e-12,1.0,0.1,1000,1000
1e-11,1.5,0.1,1000,1015
1e-10,2.0,0.1,1000,1020
1e-9,3.0,0.1,1000,1030
1e-8,5.0,0.1,1000,1050
1e-7,8.0,0.1,1000,1080
1e-6,10.0,0.1,1000,1100"""
        
        with tempfile.NamedTemporaryFile(mode='w', suffix='.csv', delete=False) as f:
            f.write(content)
            return f.name
    
    def create_mock_dsc_file(self) -> str:
        """Create a mock DSC CSV file for testing."""
        content = """Temperature,Cp
20.0,2.0
25.0,2.1
30.0,2.2
35.0,2.3
40.0,2.5
45.0,2.8
50.0,3.2
55.0,3.8
60.0,4.5
65.0,5.2
70.0,4.8
75.0,4.2
80.0,3.5
85.0,2.8
90.0,2.5"""
        
        with tempfile.NamedTemporaryFile(mode='w', suffix='.csv', delete=False) as f:
            f.write(content)
            return f.name
    
    @patch('amprenta_rag.biophysical.ingest_service.db_session')
    @patch('amprenta_rag.biophysical.ingest_service.threading.Thread')
    @patch('amprenta_rag.biophysical.ingest_service.parse_biacore_csv')
    def test_ingest_spr_file_success(self, mock_parse, mock_thread, mock_db_session):
        """Test successful SPR file ingestion."""
        # Setup mocks
        mock_db = MagicMock()
        mock_db_session.return_value = mock_db
        
        mock_spr_data = MagicMock()
        mock_spr_data.sensorgrams = [MagicMock(), MagicMock()]  # Multiple sensorgrams
        mock_spr_data.ligand_name = "Test Ligand"
        mock_spr_data.temperature = 25.0
        mock_spr_data.buffer = "PBS"
        mock_spr_data.chip_type = "CM5"
        mock_spr_data.flow_rate = 30.0
        mock_spr_data.metadata = {"test": "data"}
        mock_parse.return_value = mock_spr_data
        
        # Create mock file
        file_path = self.create_mock_spr_file()
        
        try:
            # Test ingestion
            compound_id = uuid4()
            target_id = uuid4()
            user_id = uuid4()
            
            result = ingest_spr_file(
                file_path=file_path,
                compound_id=compound_id,
                target_id=target_id,
                target_name="Test Target",
                user_id=user_id
            )
            
            # Verify parsing was called
            mock_parse.assert_called_once_with(file_path)
            
            # Verify database operations
            mock_db.add.assert_called_once()
            mock_db.commit.assert_called_once()
            mock_db.expunge.assert_called_once()
            
            # Verify experiment creation
            added_experiment = mock_db.add.call_args[0][0]
            assert isinstance(added_experiment, SPRExperiment)
            assert added_experiment.compound_id == compound_id
            assert added_experiment.target_id == target_id
            assert added_experiment.target_name == "Test Target"
            assert added_experiment.experiment_type == SPRExperimentType.KINETICS
            assert added_experiment.processing_status == "pending"
            assert added_experiment.created_by_id == user_id
            
            # Verify background thread was started
            mock_thread.assert_called_once()
            thread_args = mock_thread.call_args
            assert thread_args[1]['target'] == _process_spr_async
            assert thread_args[1]['daemon'] == True
            
        finally:
            Path(file_path).unlink()  # Cleanup
    
    @patch('amprenta_rag.biophysical.ingest_service.db_session')
    @patch('amprenta_rag.biophysical.ingest_service.threading.Thread')
    @patch('amprenta_rag.biophysical.ingest_service.parse_mst_csv')
    def test_ingest_mst_file_success(self, mock_parse, mock_thread, mock_db_session):
        """Test successful MST file ingestion."""
        # Setup mocks
        mock_db = MagicMock()
        mock_db_session.return_value = mock_db
        
        mock_mst_data = MagicMock()
        mock_mst_data.temperature = 25.0
        mock_mst_data.buffer = "PBS"
        mock_mst_data.capillary_type = "Standard"
        mock_mst_data.mst_power = 20.0
        mock_mst_data.led_power = 50.0
        mock_mst_data.metadata = {"test": "data"}
        mock_parse.return_value = mock_mst_data
        
        # Create mock file
        file_path = self.create_mock_mst_file()
        
        try:
            # Test ingestion
            compound_id = uuid4()
            result = ingest_mst_file(
                file_path=file_path,
                compound_id=compound_id,
                target_name="Test Protein"
            )
            
            # Verify parsing was called
            mock_parse.assert_called_once_with(file_path)
            
            # Verify database operations
            mock_db.add.assert_called_once()
            mock_db.commit.assert_called_once()
            
            # Verify experiment creation
            added_experiment = mock_db.add.call_args[0][0]
            assert isinstance(added_experiment, MSTExperiment)
            assert added_experiment.compound_id == compound_id
            assert added_experiment.target_name == "Test Protein"
            assert added_experiment.processing_status == "pending"
            
            # Verify background thread was started
            mock_thread.assert_called_once()
            
        finally:
            Path(file_path).unlink()  # Cleanup
    
    @patch('amprenta_rag.biophysical.ingest_service.db_session')
    @patch('amprenta_rag.biophysical.ingest_service.threading.Thread')
    @patch('amprenta_rag.biophysical.ingest_service.parse_microcal_csv')
    def test_ingest_dsc_file_success(self, mock_parse, mock_thread, mock_db_session):
        """Test successful DSC file ingestion."""
        # Setup mocks
        mock_db = MagicMock()
        mock_db_session.return_value = mock_db
        
        mock_dsc_data = MagicMock()
        mock_dsc_data.protein_name = "Test Protein"
        mock_dsc_data.scan_rate = 1.0
        mock_dsc_data.protein_concentration = 1.0
        mock_dsc_data.buffer = "PBS"
        mock_dsc_data.ph = 7.4
        mock_dsc_data.metadata = {"test": "data"}
        mock_parse.return_value = mock_dsc_data
        
        # Create mock file
        file_path = self.create_mock_dsc_file()
        
        try:
            # Test ingestion
            protein_id = uuid4()
            result = ingest_dsc_file(
                file_path=file_path,
                protein_id=protein_id,
                protein_name="Test Protein"
            )
            
            # Verify parsing was called
            mock_parse.assert_called_once_with(file_path)
            
            # Verify database operations
            mock_db.add.assert_called_once()
            mock_db.commit.assert_called_once()
            
            # Verify experiment creation
            added_experiment = mock_db.add.call_args[0][0]
            assert isinstance(added_experiment, DSCExperiment)
            assert added_experiment.target_id == protein_id
            assert added_experiment.target_name == "Test Protein"
            assert added_experiment.processing_status == "pending"
            
            # Verify background thread was started
            mock_thread.assert_called_once()
            
        finally:
            Path(file_path).unlink()  # Cleanup
    
    def test_ingest_spr_file_not_found(self):
        """Test error handling for missing SPR file."""
        with pytest.raises(FileNotFoundError, match="SPR file not found"):
            ingest_spr_file("/nonexistent/file.csv")
    
    def test_ingest_mst_file_not_found(self):
        """Test error handling for missing MST file."""
        with pytest.raises(FileNotFoundError, match="MST file not found"):
            ingest_mst_file("/nonexistent/file.xlsx")
    
    def test_ingest_dsc_file_not_found(self):
        """Test error handling for missing DSC file."""
        with pytest.raises(FileNotFoundError, match="DSC file not found"):
            ingest_dsc_file("/nonexistent/file.csv")
    
    @patch('amprenta_rag.biophysical.ingest_service.parse_biacore_csv')
    def test_ingest_spr_file_invalid_format(self, mock_parse):
        """Test error handling for invalid SPR file format."""
        mock_parse.side_effect = Exception("Invalid format")
        
        file_path = self.create_mock_spr_file()
        
        try:
            with pytest.raises(ValueError, match="Invalid SPR file format"):
                ingest_spr_file(file_path)
        finally:
            Path(file_path).unlink()  # Cleanup


class TestEntityLinking:
    """Test entity linking functions."""
    
    @patch('amprenta_rag.biophysical.ingest_service.db_session')
    def test_link_to_compound_spr(self, mock_db_session):
        """Test linking SPR experiment to compound."""
        # Setup mocks
        mock_db = MagicMock()
        mock_db_session.return_value.__enter__.return_value = mock_db
        
        mock_experiment = MagicMock()
        mock_db.query.return_value.filter.return_value.first.return_value = mock_experiment
        
        # Test linking
        experiment_id = uuid4()
        compound_id = uuid4()
        
        link_to_compound(experiment_id, compound_id, "spr")
        
        # Verify database operations
        mock_db.query.assert_called_once()
        mock_db.commit.assert_called_once()
        assert mock_experiment.compound_id == compound_id
    
    @patch('amprenta_rag.biophysical.ingest_service.db_session')
    def test_link_to_target_mst(self, mock_db_session):
        """Test linking MST experiment to target."""
        # Setup mocks
        mock_db = MagicMock()
        mock_db_session.return_value.__enter__.return_value = mock_db
        
        mock_experiment = MagicMock()
        mock_db.query.return_value.filter.return_value.first.return_value = mock_experiment
        
        # Test linking
        experiment_id = uuid4()
        target_id = uuid4()
        
        link_to_target(experiment_id, target_id, "mst")
        
        # Verify database operations
        mock_db.query.assert_called_once()
        mock_db.commit.assert_called_once()
        assert mock_experiment.target_id == target_id
    
    def test_link_to_compound_invalid_assay_type(self):
        """Test error handling for invalid assay type."""
        with pytest.raises(ValueError, match="Invalid assay type"):
            link_to_compound(uuid4(), uuid4(), "invalid")
    
    def test_link_to_target_invalid_assay_type(self):
        """Test error handling for invalid assay type."""
        with pytest.raises(ValueError, match="Invalid assay type"):
            link_to_target(uuid4(), uuid4(), "invalid")


class TestCompoundProfile:
    """Test compound biophysical profile retrieval."""
    
    @patch('amprenta_rag.biophysical.ingest_service.db_session')
    def test_get_compound_biophysical_profile(self, mock_db_session):
        """Test retrieving complete biophysical profile for a compound."""
        # Setup mocks
        mock_db = MagicMock()
        mock_db_session.return_value.__enter__.return_value = mock_db
        
        # Mock experiments
        mock_spr_exp = MagicMock()
        mock_mst_exp = MagicMock()
        mock_dsc_exp = MagicMock()
        
        # Configure query returns
        def mock_query_side_effect(model):
            if model == SPRExperiment:
                return MagicMock(filter=MagicMock(return_value=MagicMock(all=MagicMock(return_value=[mock_spr_exp]))))
            elif model == MSTExperiment:
                return MagicMock(filter=MagicMock(return_value=MagicMock(all=MagicMock(return_value=[mock_mst_exp]))))
            elif model == DSCExperiment:
                return MagicMock(filter=MagicMock(return_value=MagicMock(all=MagicMock(return_value=[mock_dsc_exp]))))
            return MagicMock()
        
        mock_db.query.side_effect = mock_query_side_effect
        
        # Test profile retrieval
        compound_id = uuid4()
        profile = get_compound_biophysical_profile(compound_id)
        
        # Verify results
        assert 'spr' in profile
        assert 'mst' in profile
        assert 'dsc' in profile
        assert len(profile['spr']) == 1
        assert len(profile['mst']) == 1
        assert len(profile['dsc']) == 1
        assert profile['spr'][0] == mock_spr_exp
        assert profile['mst'][0] == mock_mst_exp
        assert profile['dsc'][0] == mock_dsc_exp
        
        # Verify expunge was called for each experiment
        assert mock_db.expunge.call_count == 3


class TestStatusManagement:
    """Test processing status management functions."""
    
    @patch('amprenta_rag.biophysical.ingest_service.db_session')
    def test_get_processing_status(self, mock_db_session):
        """Test getting processing status of an experiment."""
        # Setup mocks
        mock_db = MagicMock()
        mock_db_session.return_value.__enter__.return_value = mock_db
        
        mock_experiment = MagicMock()
        mock_experiment.processing_status = "completed"
        mock_db.query.return_value.filter.return_value.first.return_value = mock_experiment
        
        # Test status retrieval
        experiment_id = uuid4()
        status = get_processing_status(experiment_id, "spr")
        
        assert status == "completed"
        mock_db.query.assert_called_once()
    
    @patch('amprenta_rag.biophysical.ingest_service.db_session')
    def test_get_processing_status_not_found(self, mock_db_session):
        """Test error handling for experiment not found."""
        # Setup mocks
        mock_db = MagicMock()
        mock_db_session.return_value.__enter__.return_value = mock_db
        mock_db.query.return_value.filter.return_value.first.return_value = None
        
        # Test error case
        with pytest.raises(ValueError, match="SPR experiment not found"):
            get_processing_status(uuid4(), "spr")
    
    @patch('amprenta_rag.biophysical.ingest_service.db_session')
    @patch('amprenta_rag.biophysical.ingest_service.threading.Thread')
    def test_reprocess_experiment(self, mock_thread, mock_db_session):
        """Test reprocessing an existing experiment."""
        # Setup mocks
        mock_db = MagicMock()
        mock_db_session.return_value.__enter__.return_value = mock_db
        
        mock_experiment = MagicMock()
        mock_db.query.return_value.filter.return_value.first.return_value = mock_experiment
        
        # Test reprocessing
        experiment_id = uuid4()
        reprocess_experiment(experiment_id, "dsc")
        
        # Verify status reset
        assert mock_experiment.processing_status == "pending"
        assert mock_experiment.error_message is None
        mock_db.commit.assert_called_once()
        
        # Verify background thread was started
        mock_thread.assert_called_once()
        thread_args = mock_thread.call_args
        assert thread_args[1]['target'] == _process_dsc_async
        assert thread_args[1]['args'][0] == experiment_id
    
    def test_get_processing_status_invalid_assay_type(self):
        """Test error handling for invalid assay type."""
        with pytest.raises(ValueError, match="Invalid assay type"):
            get_processing_status(uuid4(), "invalid")
    
    def test_reprocess_experiment_invalid_assay_type(self):
        """Test error handling for invalid assay type."""
        with pytest.raises(ValueError, match="Invalid assay type"):
            reprocess_experiment(uuid4(), "invalid")


class TestBackgroundProcessing:
    """Test background processing functions (basic functionality)."""
    
    @patch('amprenta_rag.biophysical.ingest_service.db_session')
    @patch('amprenta_rag.biophysical.ingest_service.parse_biacore_csv')
    @patch('amprenta_rag.biophysical.ingest_service.fit_1_to_1_langmuir')
    def test_process_spr_async_success(self, mock_fit, mock_parse, mock_db_session):
        """Test successful SPR background processing."""
        # Setup mocks
        mock_db = MagicMock()
        mock_db_session.return_value.__enter__.return_value = mock_db
        
        mock_experiment = MagicMock()
        mock_experiment.raw_data_path = "test.csv"
        mock_db.query.return_value.filter.return_value.first.return_value = mock_experiment
        
        mock_sensorgram = MagicMock()
        mock_sensorgram.cycle = 1
        mock_sensorgram.concentration = 100e-9
        mock_sensorgram.time = np.array([0, 1, 2, 3, 4])
        mock_sensorgram.response = np.array([0, 10, 20, 15, 10])
        mock_sensorgram.association_start = 0
        mock_sensorgram.dissociation_start = 2
        mock_sensorgram.flow_cell = 1
        mock_sensorgram.reference_subtracted = True
        
        mock_spr_data = MagicMock()
        mock_spr_data.sensorgrams = [mock_sensorgram]
        mock_parse.return_value = mock_spr_data
        
        mock_kinetic_fit = MagicMock()
        mock_kinetic_fit.ka = 1e5
        mock_kinetic_fit.kd = 1e-3
        mock_kinetic_fit.kd_affinity = 1e-8
        mock_kinetic_fit.rmax = 100
        mock_kinetic_fit.chi_squared = 5.0
        mock_fit.return_value = mock_kinetic_fit
        
        # Test background processing
        experiment_id = uuid4()
        _process_spr_async(experiment_id)
        
        # Verify parsing and analysis were called
        mock_parse.assert_called_once_with("test.csv")
        mock_fit.assert_called_once()
        
        # Verify database updates
        assert mock_experiment.processing_status == "completed"
        assert mock_experiment.ka == 1e5
        assert mock_experiment.kd == 1e-3
        assert mock_experiment.kd_equilibrium == 1e-8
        mock_db.add.assert_called()  # Sensorgram added
        mock_db.commit.assert_called()
    
    @patch('amprenta_rag.biophysical.ingest_service.db_session')
    def test_process_spr_async_experiment_not_found(self, mock_db_session):
        """Test error handling when experiment is not found."""
        # Setup mocks
        mock_db = MagicMock()
        mock_db_session.return_value.__enter__.return_value = mock_db
        mock_db.query.return_value.filter.return_value.first.return_value = None
        
        # Test with non-existent experiment
        experiment_id = uuid4()
        
        # Should not raise exception (logs error and returns)
        _process_spr_async(experiment_id)
        
        # Verify query was attempted
        mock_db.query.assert_called_once()


class TestIntegrationWorkflow:
    """Integration tests for complete workflows."""
    
    @patch('amprenta_rag.biophysical.ingest_service.db_session')
    @patch('amprenta_rag.biophysical.ingest_service.threading.Thread')
    @patch('amprenta_rag.biophysical.ingest_service.parse_biacore_csv')
    def test_complete_spr_workflow(self, mock_parse, mock_thread, mock_db_session):
        """Test complete SPR workflow from ingestion to linking."""
        # Setup mocks for ingestion
        mock_db = MagicMock()
        mock_db_session.return_value = mock_db
        mock_db_session.return_value.__enter__.return_value = mock_db
        
        mock_spr_data = MagicMock()
        mock_spr_data.sensorgrams = [MagicMock()]
        mock_spr_data.ligand_name = "Test Ligand"
        mock_spr_data.temperature = 25.0
        mock_spr_data.buffer = "PBS"
        mock_spr_data.chip_type = "CM5"
        mock_spr_data.flow_rate = 30.0
        mock_spr_data.metadata = {}
        mock_parse.return_value = mock_spr_data
        
        # Create mock file
        file_path = self.create_mock_spr_file()
        
        try:
            # Step 1: Ingest file
            compound_id = uuid4()
            spr_experiment = ingest_spr_file(
                file_path=file_path,
                compound_id=compound_id
            )
            
            # Verify ingestion
            mock_db.add.assert_called_once()
            mock_thread.assert_called_once()
            
            # Step 2: Link to compound (already linked during ingestion)
            experiment_id = spr_experiment.id
            
            # Mock experiment for linking test
            mock_experiment = MagicMock()
            mock_db.query.return_value.filter.return_value.first.return_value = mock_experiment
            
            # Test additional linking
            target_id = uuid4()
            link_to_target(experiment_id, target_id, "spr")
            
            # Verify linking
            assert mock_experiment.target_id == target_id
            
            # Step 3: Get compound profile
            profile = get_compound_biophysical_profile(compound_id)
            
            # Verify profile structure
            assert 'spr' in profile
            assert 'mst' in profile
            assert 'dsc' in profile
            
        finally:
            Path(file_path).unlink()  # Cleanup
    
    def create_mock_spr_file(self) -> str:
        """Create a mock SPR CSV file for testing."""
        content = """Instrument: Biacore T200
Time	Response	Concentration	Cycle
0.0	0.0	100.0	1
1.0	10.0	100.0	1
2.0	25.0	100.0	1"""
        
        with tempfile.NamedTemporaryFile(mode='w', suffix='.csv', delete=False) as f:
            f.write(content)
            return f.name
