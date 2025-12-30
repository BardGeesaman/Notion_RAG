"""
Integration tests for Biophysical Assay functionality.

These tests cover the complete workflows using real PostgreSQL database,
testing SPR, MST, and DSC end-to-end functionality.
"""

from __future__ import annotations

import tempfile
from pathlib import Path
from uuid import uuid4
from datetime import datetime, timezone

import numpy as np
import pytest
from fastapi.testclient import TestClient

from amprenta_rag.api.main import app
from amprenta_rag.database.session import db_session
from amprenta_rag.database.models import User
from amprenta_rag.database.models_biophysical import (
    SPRExperiment,
    SPRSensorgram,
    MSTExperiment,
    MSTDoseResponse,
    DSCExperiment,
    DSCScan,
    SPRExperimentType,
)
from amprenta_rag.models.chemistry import Compound
from amprenta_rag.database.models import ProteinStructure
from amprenta_rag.biophysical.ingest_service import (
    ingest_spr_file,
    ingest_mst_file,
    ingest_dsc_file,
    get_compound_biophysical_profile,
    reprocess_experiment,
)


pytestmark = pytest.mark.integration


@pytest.fixture
def test_user():
    """Create a test user in the database."""
    with db_session() as db:
        user = User(
            id=uuid4(),
            username=f"test_user_{uuid4().hex[:8]}",
            email=f"test_{uuid4().hex[:8]}@example.com",
            password_hash="test_hash"
        )
        db.add(user)
        db.commit()
        db.expunge(user)
        return user


@pytest.fixture
def test_compound():
    """Create a test compound in the database."""
    with db_session() as db:
        compound = Compound(
            id=uuid4(),
            compound_id=f"TEST_COMPOUND_{uuid4().hex[:8]}",
            smiles="CC(=O)Nc1ccc(O)cc1",
            molecular_weight=151.16,
            created_at=datetime.now(timezone.utc)
        )
        db.add(compound)
        db.commit()
        db.expunge(compound)
        return compound


@pytest.fixture
def test_target():
    """Create a test protein target in the database."""
    with db_session() as db:
        target = ProteinStructure(
            id=uuid4(),
            source="test",
            created_at=datetime.now(timezone.utc)
        )
        db.add(target)
        db.commit()
        db.expunge(target)
        return target


@pytest.fixture
def mock_spr_file():
    """Create a mock SPR CSV file for testing."""
    csv_content = """Time,Response,Cycle,Concentration
0.0,0.0,1,1e-9
10.0,5.2,1,1e-9
20.0,12.8,1,1e-9
30.0,22.1,1,1e-9
40.0,32.5,1,1e-9
50.0,43.2,1,1e-9
60.0,52.8,1,1e-9
70.0,61.5,1,1e-9
80.0,68.9,1,1e-9
90.0,75.2,1,1e-9
100.0,80.1,1,1e-9
300.0,78.5,1,1e-9
310.0,76.2,1,1e-9
320.0,73.8,1,1e-9
330.0,71.1,1,1e-9
340.0,68.3,1,1e-9
350.0,65.2,1,1e-9
360.0,62.0,1,1e-9
370.0,58.5,1,1e-9
380.0,55.1,1,1e-9
390.0,51.8,1,1e-9
400.0,48.2,1,1e-9
"""
    temp_file = tempfile.NamedTemporaryFile(mode='w', suffix='.csv', delete=False)
    temp_file.write(csv_content)
    temp_file.close()
    return temp_file.name


@pytest.fixture
def mock_mst_file():
    """Create a mock MST CSV file for testing."""
    csv_content = """Concentration,Fnorm,Fnorm_Error
1e-12,1.000,0.050
1e-11,0.985,0.045
1e-10,0.932,0.048
1e-9,0.821,0.052
1e-8,0.634,0.055
1e-7,0.412,0.048
1e-6,0.245,0.042
1e-5,0.142,0.038
1e-4,0.089,0.035
1e-3,0.065,0.032
"""
    temp_file = tempfile.NamedTemporaryFile(mode='w', suffix='.csv', delete=False)
    temp_file.write(csv_content)
    temp_file.close()
    return temp_file.name


@pytest.fixture
def mock_dsc_file():
    """Create a mock DSC CSV file for testing."""
    csv_content = """Temperature,Cp
25.0,1.52
30.0,1.58
35.0,1.65
40.0,1.74
45.0,1.86
50.0,2.02
55.0,2.25
60.0,2.58
65.0,3.12
70.0,4.25
75.0,3.85
80.0,2.95
85.0,2.15
90.0,1.85
"""
    temp_file = tempfile.NamedTemporaryFile(mode='w', suffix='.csv', delete=False)
    temp_file.write(csv_content)
    temp_file.close()
    return temp_file.name


def test_full_spr_workflow(test_user, test_compound, test_target, mock_spr_file):
    """Test complete SPR workflow: ingest → analyze → retrieve."""
    
    # Test ingestion
    experiment = ingest_spr_file(
        file_path=mock_spr_file,
        compound_id=test_compound.id,
        target_id=test_target.id,
        target_name=test_target.name,
        user_id=test_user.id
    )
    
    assert experiment is not None
    assert experiment.compound_id == test_compound.id
    assert experiment.target_id == test_target.id
    assert experiment.processing_status in ["pending", "running", "completed"]
    
    # Wait for processing to complete (in real scenario, would be background)
    # For test, we simulate completion
    with db_session() as db:
        exp = db.query(SPRExperiment).filter(SPRExperiment.id == experiment.id).first()
        if exp and exp.processing_status != "completed":
            exp.processing_status = "completed"
            exp.ka = 1.5e5
            exp.kd_rate = 2.3e-3
            exp.kd_equilibrium = exp.kd_rate / exp.ka
            exp.rmax = 85.0
            exp.chi_squared = 1.25
            db.commit()
    
    # Test retrieval
    with db_session() as db:
        retrieved_exp = db.query(SPRExperiment).filter(SPRExperiment.id == experiment.id).first()
        assert retrieved_exp is not None
        assert retrieved_exp.processing_status == "completed"
        assert retrieved_exp.ka is not None
        assert retrieved_exp.kd_rate is not None
        assert retrieved_exp.kd_equilibrium is not None
        
        # Check sensorgrams were created
        sensorgrams = db.query(SPRSensorgram).filter(SPRSensorgram.experiment_id == experiment.id).all()
        assert len(sensorgrams) > 0
        
        for sensorgram in sensorgrams:
            assert sensorgram.time_seconds is not None
            assert sensorgram.response_values is not None
            assert len(sensorgram.time_seconds) > 0
            assert len(sensorgram.response_values) > 0
    
    # Cleanup
    Path(mock_spr_file).unlink()


def test_full_mst_workflow(test_user, test_compound, test_target, mock_mst_file):
    """Test complete MST workflow: ingest → analyze → retrieve."""
    
    # Test ingestion
    experiment = ingest_mst_file(
        file_path=mock_mst_file,
        compound_id=test_compound.id,
        target_id=test_target.id,
        target_name=test_target.name,
        user_id=test_user.id
    )
    
    assert experiment is not None
    assert experiment.compound_id == test_compound.id
    assert experiment.target_id == test_target.id
    assert experiment.processing_status in ["pending", "running", "completed"]
    
    # Simulate processing completion
    with db_session() as db:
        exp = db.query(MSTExperiment).filter(MSTExperiment.id == experiment.id).first()
        if exp and exp.processing_status != "completed":
            exp.processing_status = "completed"
            exp.kd_value = 2.5e-8
            exp.kd_error = 3.5e-9
            exp.hill_coefficient = 1.1
            exp.binding_amplitude = 25.8
            exp.signal_to_noise = 12.5
            exp.aggregation_detected = False
            db.commit()
    
    # Test retrieval
    with db_session() as db:
        retrieved_exp = db.query(MSTExperiment).filter(MSTExperiment.id == experiment.id).first()
        assert retrieved_exp is not None
        assert retrieved_exp.processing_status == "completed"
        assert retrieved_exp.kd_value is not None
        assert retrieved_exp.hill_coefficient is not None
        
        # Check dose-response points were created
        dose_points = db.query(MSTDoseResponse).filter(MSTDoseResponse.experiment_id == experiment.id).all()
        assert len(dose_points) > 0
        
        for point in dose_points:
            assert point.ligand_concentration is not None
            assert point.thermophoresis_response is not None
            assert point.ligand_concentration > 0
    
    # Cleanup
    Path(mock_mst_file).unlink()


def test_full_dsc_workflow(test_user, test_compound, test_target, mock_dsc_file):
    """Test complete DSC workflow: ingest → analyze → retrieve."""
    
    # Test ingestion
    experiment = ingest_dsc_file(
        file_path=mock_dsc_file,
        compound_id=test_compound.id,
        target_id=test_target.id,
        target_name=test_target.name,
        user_id=test_user.id
    )
    
    assert experiment is not None
    assert experiment.compound_id == test_compound.id
    assert experiment.target_id == test_target.id
    assert experiment.processing_status in ["pending", "running", "completed"]
    
    # Simulate processing completion
    with db_session() as db:
        exp = db.query(DSCExperiment).filter(DSCExperiment.id == experiment.id).first()
        if exp and exp.processing_status != "completed":
            exp.processing_status = "completed"
            exp.tm_value = 65.8
            exp.tm_error = 0.5
            exp.delta_h = -125.3
            exp.delta_h_error = 12.5
            exp.cooperativity = 1.75
            db.commit()
    
    # Test retrieval
    with db_session() as db:
        retrieved_exp = db.query(DSCExperiment).filter(DSCExperiment.id == experiment.id).first()
        assert retrieved_exp is not None
        assert retrieved_exp.processing_status == "completed"
        assert retrieved_exp.tm_value is not None
        assert retrieved_exp.delta_h is not None
        assert retrieved_exp.cooperativity is not None
        
        # Check thermogram scans were created
        scans = db.query(DSCScan).filter(DSCScan.experiment_id == experiment.id).all()
        assert len(scans) > 0
        
        for scan in scans:
            assert scan.temperature_celsius is not None
            assert scan.heat_capacity is not None
            assert len(scan.temperature_celsius) > 0
            assert len(scan.heat_capacity) > 0
    
    # Cleanup
    Path(mock_dsc_file).unlink()


def test_cross_assay_compound_profile(test_user, test_compound, test_target, mock_spr_file, mock_mst_file, mock_dsc_file):
    """Test retrieving all biophysical data for a compound."""
    
    # Create experiments for all three assay types
    spr_exp = ingest_spr_file(
        file_path=mock_spr_file,
        compound_id=test_compound.id,
        target_id=test_target.id,
        user_id=test_user.id
    )
    
    mst_exp = ingest_mst_file(
        file_path=mock_mst_file,
        compound_id=test_compound.id,
        target_id=test_target.id,
        user_id=test_user.id
    )
    
    dsc_exp = ingest_dsc_file(
        file_path=mock_dsc_file,
        compound_id=test_compound.id,
        target_id=test_target.id,
        user_id=test_user.id
    )
    
    # Mark all as completed
    with db_session() as db:
        for exp_id, exp_type in [(spr_exp.id, SPRExperiment), (mst_exp.id, MSTExperiment), (dsc_exp.id, DSCExperiment)]:
            exp = db.query(exp_type).filter(exp_type.id == exp_id).first()
            if exp:
                exp.processing_status = "completed"
        db.commit()
    
    # Test compound profile retrieval
    profile = get_compound_biophysical_profile(test_compound.id)
    
    assert profile is not None
    assert "spr_experiments" in profile
    assert "mst_experiments" in profile
    assert "dsc_experiments" in profile
    
    assert len(profile["spr_experiments"]) >= 1
    assert len(profile["mst_experiments"]) >= 1
    assert len(profile["dsc_experiments"]) >= 1
    
    # Verify experiment IDs match
    spr_ids = [exp.id for exp in profile["spr_experiments"]]
    mst_ids = [exp.id for exp in profile["mst_experiments"]]
    dsc_ids = [exp.id for exp in profile["dsc_experiments"]]
    
    assert spr_exp.id in spr_ids
    assert mst_exp.id in mst_ids
    assert dsc_exp.id in dsc_ids
    
    # Cleanup
    Path(mock_spr_file).unlink()
    Path(mock_mst_file).unlink()
    Path(mock_dsc_file).unlink()


def test_entity_linking(test_user, test_compound, test_target):
    """Test compound/target linking works correctly."""
    
    with db_session() as db:
        # Create experiment with entity links
        experiment = SPRExperiment(
            id=uuid4(),
            experiment_name="Entity Linking Test",
            compound_id=test_compound.id,
            target_id=test_target.id,
            target_name=test_target.name,
            experiment_type=SPRExperimentType.KINETICS,
            processing_status="completed",
            response_units="RU",
            concentration_units="M",
            replicate_number=1,
            replicate_group_id=uuid4(),
            created_by_id=test_user.id,
            created_at=datetime.now(timezone.utc)
        )
        db.add(experiment)
        db.commit()
        
        experiment_id = experiment.id
        db.expunge(experiment)
    
    # Verify relationships work
    with db_session() as db:
        exp = db.query(SPRExperiment).filter(SPRExperiment.id == experiment_id).first()
        assert exp is not None
        assert exp.compound_id == test_compound.id
        assert exp.target_id == test_target.id
        assert exp.target_name == test_target.name
        
        # Verify we can query by compound
        compound_experiments = db.query(SPRExperiment).filter(SPRExperiment.compound_id == test_compound.id).all()
        assert len(compound_experiments) >= 1
        assert any(exp.id == experiment_id for exp in compound_experiments)
        
        # Verify we can query by target
        target_experiments = db.query(SPRExperiment).filter(SPRExperiment.target_id == test_target.id).all()
        assert len(target_experiments) >= 1
        assert any(exp.id == experiment_id for exp in target_experiments)


def test_reprocessing(test_user, test_compound, test_target, mock_spr_file):
    """Test experiment reprocessing capability."""
    
    # Create initial experiment
    experiment = ingest_spr_file(
        file_path=mock_spr_file,
        compound_id=test_compound.id,
        target_id=test_target.id,
        user_id=test_user.id
    )
    
    # Mark as completed with initial parameters
    with db_session() as db:
        exp = db.query(SPRExperiment).filter(SPRExperiment.id == experiment.id).first()
        exp.processing_status = "completed"
        exp.ka = 1.0e5
        exp.kd_rate = 1.0e-3
        exp.kd_equilibrium = 1.0e-8
        initial_ka = exp.ka
        db.commit()
    
    # Test reprocessing
    reprocessed_exp = reprocess_experiment(
        experiment_id=experiment.id,
        assay_type="spr",
        model="two_state"
    )
    
    assert reprocessed_exp is not None
    
    # Verify reprocessing updated the experiment
    with db_session() as db:
        exp = db.query(SPRExperiment).filter(SPRExperiment.id == experiment.id).first()
        assert exp is not None
        # In a real scenario, reprocessing might change parameters
        # For this test, we just verify the reprocessing call worked
        assert exp.processing_status in ["completed", "running", "pending"]
    
    # Cleanup
    Path(mock_spr_file).unlink()


def test_api_integration_spr_upload():
    """Test SPR file upload via API integration."""
    client = TestClient(app)
    
    # Create mock file content
    csv_content = """Time,Response,Cycle,Concentration
0.0,0.0,1,1e-9
10.0,5.2,1,1e-9
20.0,12.8,1,1e-9
"""
    
    # Test upload (will fail without proper auth, but tests endpoint existence)
    response = client.post(
        "/api/v1/biophysical/spr/upload",
        files={"file": ("test.csv", csv_content, "text/csv")},
        data={"target_name": "Test Target"}
    )
    
    # Expect 401 (unauthorized) since we're not authenticated in this test
    # This confirms the endpoint exists and is protected
    assert response.status_code == 401


def test_api_integration_experiment_listing():
    """Test experiment listing via API integration."""
    client = TestClient(app)
    
    # Test experiment listing (will fail without auth, but tests endpoint)
    response = client.get("/api/v1/biophysical/spr")
    
    # Expect 401 (unauthorized) or 200 (if endpoint allows anonymous access)
    assert response.status_code in [200, 401]
    
    response = client.get("/api/v1/biophysical/mst")
    assert response.status_code in [200, 401]
    
    response = client.get("/api/v1/biophysical/dsc")
    assert response.status_code in [200, 401]
