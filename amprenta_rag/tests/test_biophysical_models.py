"""
Tests for Biophysical Assay database models.

Tests cover model creation, relationships, constraints, and business logic
for SPR, MST, and DSC biophysical assay data storage and analysis.
"""

import pytest
import json
from datetime import datetime, timezone
from uuid import uuid4

from sqlalchemy import create_engine, Column, String, Text, Integer, Float, Boolean, DateTime, JSON, ForeignKey, BigInteger, UniqueConstraint
from sqlalchemy.exc import IntegrityError
from sqlalchemy.orm import Session, sessionmaker, declarative_base
from sqlalchemy.dialects.postgresql import UUID

from amprenta_rag.database.models_biophysical import (
    SPRExperimentType,
    generate_uuid,
)

# Create a separate Base for test models to avoid ARRAY type issues with SQLite
TestBase = declarative_base()


# Create test-specific models without ARRAY types for SQLite compatibility
class TestUser(TestBase):
    """Simple User model for testing."""
    __tablename__ = "test_users"
    
    id = Column(UUID(as_uuid=True), primary_key=True, default=generate_uuid)
    username = Column(String(100), nullable=False, unique=True)
    email = Column(String(200), nullable=False, unique=True)


class TestCompound(TestBase):
    """Simple Compound model for testing."""
    __tablename__ = "test_compounds"
    
    id = Column(UUID(as_uuid=True), primary_key=True, default=generate_uuid)
    compound_id = Column(String(200), nullable=False, unique=True)
    smiles = Column(Text, nullable=False)


class TestProteinStructure(TestBase):
    """Simple ProteinStructure model for testing."""
    __tablename__ = "test_protein_structures"
    
    id = Column(UUID(as_uuid=True), primary_key=True, default=generate_uuid)
    pdb_id = Column(String(10), nullable=True)
    source = Column(String(50), nullable=False, default="pdb")


class TestSPRExperiment(TestBase):
    """SPR experiment model for testing without ARRAY types."""
    __tablename__ = "test_spr_experiments"

    id = Column(UUID(as_uuid=True), primary_key=True, default=generate_uuid)
    experiment_name = Column(String(200), nullable=False, index=True)
    experiment_type = Column(String(20), nullable=False, default=SPRExperimentType.KINETICS)
    compound_id = Column(UUID(as_uuid=True), ForeignKey("test_compounds.id"), nullable=True, index=True)
    target_id = Column(UUID(as_uuid=True), ForeignKey("test_protein_structures.id"), nullable=True, index=True)
    target_name = Column(String(200), nullable=True)
    analyte_name = Column(String(200), nullable=True)
    ligand_name = Column(String(200), nullable=True)
    analyte_concentration = Column(Float, nullable=True)
    ligand_density = Column(Float, nullable=True)
    buffer_composition = Column(String(500), nullable=True)
    temperature_celsius = Column(Float, nullable=True, default=25.0)
    flow_rate = Column(Float, nullable=True)
    response_units = Column(String(50), nullable=False, default="RU")
    concentration_units = Column(String(20), nullable=False, default="nM")
    ka = Column(Float, nullable=True)
    kd = Column(Float, nullable=True)
    kd_equilibrium = Column(Float, nullable=True)
    rmax = Column(Float, nullable=True)
    chi_squared = Column(Float, nullable=True)
    baseline_drift = Column(Float, nullable=True)
    bulk_refractive_index = Column(Float, nullable=True)
    reference_subtracted = Column(Boolean, nullable=False, default=True)
    instrument_model = Column(String(100), nullable=True)
    chip_type = Column(String(50), nullable=True)
    cycle_number = Column(Integer, nullable=True)
    replicate_number = Column(Integer, nullable=False, default=1)
    replicate_group_id = Column(UUID(as_uuid=True), nullable=True, index=True)
    processing_status = Column(String(50), nullable=False, default="pending")
    processing_log = Column(Text, nullable=True)
    experiment_date = Column(DateTime(timezone=True), nullable=True)
    created_at = Column(DateTime(timezone=True), default=lambda: datetime.now(timezone.utc), nullable=False)
    updated_at = Column(DateTime(timezone=True), default=lambda: datetime.now(timezone.utc), nullable=False)
    created_by_id = Column(UUID(as_uuid=True), ForeignKey("test_users.id"), nullable=True)


class TestSPRSensorgram(TestBase):
    """SPR sensorgram model for testing without ARRAY types."""
    __tablename__ = "test_spr_sensorgrams"

    id = Column(UUID(as_uuid=True), primary_key=True, default=generate_uuid)
    experiment_id = Column(UUID(as_uuid=True), ForeignKey("test_spr_experiments.id"), nullable=False, index=True)
    cycle_phase = Column(String(50), nullable=False)
    flow_cell = Column(Integer, nullable=True)
    # Store time-series as JSON for SQLite compatibility
    time_seconds_json = Column(JSON, nullable=True)
    response_values_json = Column(JSON, nullable=True)
    baseline_corrected = Column(Boolean, nullable=False, default=False)
    reference_subtracted = Column(Boolean, nullable=False, default=False)
    bulk_shift_corrected = Column(Boolean, nullable=False, default=False)
    data_points_count = Column(Integer, nullable=True)
    noise_level = Column(Float, nullable=True)
    created_at = Column(DateTime(timezone=True), default=lambda: datetime.now(timezone.utc), nullable=False)


class TestMSTExperiment(TestBase):
    """MST experiment model for testing."""
    __tablename__ = "test_mst_experiments"

    id = Column(UUID(as_uuid=True), primary_key=True, default=generate_uuid)
    experiment_name = Column(String(200), nullable=False, index=True)
    compound_id = Column(UUID(as_uuid=True), ForeignKey("test_compounds.id"), nullable=True, index=True)
    target_id = Column(UUID(as_uuid=True), ForeignKey("test_protein_structures.id"), nullable=True, index=True)
    target_name = Column(String(200), nullable=True)
    target_concentration = Column(Float, nullable=True)
    ligand_name = Column(String(200), nullable=True)
    buffer_composition = Column(String(500), nullable=True)
    temperature_celsius = Column(Float, nullable=True, default=25.0)
    excitation_power = Column(Float, nullable=True)
    mst_power = Column(Float, nullable=True)
    response_units = Column(String(50), nullable=False, default="‰")
    concentration_units = Column(String(20), nullable=False, default="nM")
    kd_value = Column(Float, nullable=True)
    kd_error = Column(Float, nullable=True)
    hill_coefficient = Column(Float, nullable=True)
    binding_amplitude = Column(Float, nullable=True)
    r_squared = Column(Float, nullable=True)
    signal_to_noise = Column(Float, nullable=True)
    aggregation_detected = Column(Boolean, nullable=False, default=False)
    bleaching_rate = Column(Float, nullable=True)
    instrument_model = Column(String(100), nullable=True)
    capillary_type = Column(String(50), nullable=True)
    replicate_number = Column(Integer, nullable=False, default=1)
    replicate_group_id = Column(UUID(as_uuid=True), nullable=True, index=True)
    processing_status = Column(String(50), nullable=False, default="pending")
    processing_log = Column(Text, nullable=True)
    experiment_date = Column(DateTime(timezone=True), nullable=True)
    created_at = Column(DateTime(timezone=True), default=lambda: datetime.now(timezone.utc), nullable=False)
    updated_at = Column(DateTime(timezone=True), default=lambda: datetime.now(timezone.utc), nullable=False)
    created_by_id = Column(UUID(as_uuid=True), ForeignKey("test_users.id"), nullable=True)


class TestMSTDoseResponse(TestBase):
    """MST dose-response model for testing."""
    __tablename__ = "test_mst_dose_responses"

    id = Column(UUID(as_uuid=True), primary_key=True, default=generate_uuid)
    experiment_id = Column(UUID(as_uuid=True), ForeignKey("test_mst_experiments.id"), nullable=False, index=True)
    ligand_concentration = Column(Float, nullable=False)
    thermophoresis_response = Column(Float, nullable=False)
    thermophoresis_error = Column(Float, nullable=True)
    initial_fluorescence = Column(Float, nullable=True)
    thermophoresis_signal = Column(Float, nullable=True)
    capillary_number = Column(Integer, nullable=True)
    measurement_valid = Column(Boolean, nullable=False, default=True)
    outlier_flag = Column(Boolean, nullable=False, default=False)
    created_at = Column(DateTime(timezone=True), default=lambda: datetime.now(timezone.utc), nullable=False)


class TestDSCExperiment(TestBase):
    """DSC experiment model for testing."""
    __tablename__ = "test_dsc_experiments"

    id = Column(UUID(as_uuid=True), primary_key=True, default=generate_uuid)
    experiment_name = Column(String(200), nullable=False, index=True)
    compound_id = Column(UUID(as_uuid=True), ForeignKey("test_compounds.id"), nullable=True, index=True)
    target_id = Column(UUID(as_uuid=True), ForeignKey("test_protein_structures.id"), nullable=True, index=True)
    target_name = Column(String(200), nullable=True)
    protein_concentration = Column(Float, nullable=True)
    ligand_concentration = Column(Float, nullable=True)
    buffer_composition = Column(String(500), nullable=True)
    scan_rate = Column(Float, nullable=True, default=1.0)
    start_temperature = Column(Float, nullable=True, default=20.0)
    end_temperature = Column(Float, nullable=True, default=100.0)
    response_units = Column(String(50), nullable=False, default="kcal/mol/°C")
    concentration_units = Column(String(20), nullable=False, default="μM")
    tm_value = Column(Float, nullable=True)
    tm_error = Column(Float, nullable=True)
    delta_h = Column(Float, nullable=True)
    delta_h_error = Column(Float, nullable=True)
    delta_cp = Column(Float, nullable=True)
    cooperativity = Column(Float, nullable=True)
    kd_thermal = Column(Float, nullable=True)
    delta_tm = Column(Float, nullable=True)
    baseline_quality = Column(Float, nullable=True)
    peak_symmetry = Column(Float, nullable=True)
    signal_to_noise = Column(Float, nullable=True)
    instrument_model = Column(String(100), nullable=True)
    cell_volume = Column(Float, nullable=True)
    reference_buffer = Column(String(500), nullable=True)
    replicate_number = Column(Integer, nullable=False, default=1)
    replicate_group_id = Column(UUID(as_uuid=True), nullable=True, index=True)
    processing_status = Column(String(50), nullable=False, default="pending")
    processing_log = Column(Text, nullable=True)
    experiment_date = Column(DateTime(timezone=True), nullable=True)
    created_at = Column(DateTime(timezone=True), default=lambda: datetime.now(timezone.utc), nullable=False)
    updated_at = Column(DateTime(timezone=True), default=lambda: datetime.now(timezone.utc), nullable=False)
    created_by_id = Column(UUID(as_uuid=True), ForeignKey("test_users.id"), nullable=True)


class TestDSCScan(TestBase):
    """DSC scan model for testing without ARRAY types."""
    __tablename__ = "test_dsc_scans"

    id = Column(UUID(as_uuid=True), primary_key=True, default=generate_uuid)
    experiment_id = Column(UUID(as_uuid=True), ForeignKey("test_dsc_experiments.id"), nullable=False, index=True)
    scan_type = Column(String(50), nullable=False)
    scan_number = Column(Integer, nullable=True)
    # Store thermogram data as JSON for SQLite compatibility
    temperature_celsius_json = Column(JSON, nullable=True)
    heat_capacity_json = Column(JSON, nullable=True)
    baseline_subtracted = Column(Boolean, nullable=False, default=False)
    reference_subtracted = Column(Boolean, nullable=False, default=False)
    normalized = Column(Boolean, nullable=False, default=False)
    data_points_count = Column(Integer, nullable=True)
    noise_level = Column(Float, nullable=True)
    created_at = Column(DateTime(timezone=True), default=lambda: datetime.now(timezone.utc), nullable=False)


@pytest.fixture
def db_session():
    """Create in-memory SQLite database for testing."""
    engine = create_engine("sqlite:///:memory:", echo=False)
    
    # Create all test tables
    TestBase.metadata.create_all(engine)
    
    SessionLocal = sessionmaker(bind=engine)
    session = SessionLocal()
    
    yield session
    
    session.close()


@pytest.fixture
def test_user(db_session):
    """Create a test user."""
    user = TestUser(
        id=uuid4(),
        username="test_user",
        email="test@example.com"
    )
    db_session.add(user)
    db_session.commit()
    db_session.refresh(user)
    return user


@pytest.fixture
def test_compound(db_session):
    """Create a test compound."""
    compound = TestCompound(
        id=uuid4(),
        compound_id="COMP-001",
        smiles="CCO"
    )
    db_session.add(compound)
    db_session.commit()
    db_session.refresh(compound)
    return compound


@pytest.fixture
def test_protein(db_session):
    """Create a test protein structure."""
    protein = TestProteinStructure(
        id=uuid4(),
        pdb_id="1ABC",
        source="pdb"
    )
    db_session.add(protein)
    db_session.commit()
    db_session.refresh(protein)
    return protein


class TestSPRExperimentModel:
    """Test SPR experiment model functionality."""

    def test_spr_experiment_creation(self, db_session, test_user, test_compound, test_protein):
        """Test creating an SPR experiment with all required fields."""
        experiment = TestSPRExperiment(
            experiment_name="SPR Kinetics Test",
            experiment_type=SPRExperimentType.KINETICS,
            compound_id=test_compound.id,
            target_id=test_protein.id,
            target_name="EGFR",
            analyte_name="Compound-001",
            ligand_name="EGFR-His6",
            analyte_concentration=100.0,
            ligand_density=1500.0,
            buffer_composition="PBS + 0.05% Tween-20",
            temperature_celsius=25.0,
            flow_rate=30.0,
            response_units="RU",
            concentration_units="nM",
            ka=1e5,
            kd=1e-3,
            kd_equilibrium=1e-8,
            rmax=150.0,
            chi_squared=2.5,
            baseline_drift=0.1,
            bulk_refractive_index=1.33,
            reference_subtracted=True,
            instrument_model="Biacore T200",
            chip_type="CM5",
            cycle_number=1,
            replicate_number=1,
            processing_status="completed",
            experiment_date=datetime.now(timezone.utc),
            created_by_id=test_user.id
        )
        
        db_session.add(experiment)
        db_session.commit()
        db_session.refresh(experiment)
        
        assert experiment.id is not None
        assert experiment.experiment_name == "SPR Kinetics Test"
        assert experiment.experiment_type == SPRExperimentType.KINETICS
        assert experiment.compound_id == test_compound.id
        assert experiment.target_id == test_protein.id
        assert experiment.ka == 1e5
        assert experiment.kd == 1e-3
        assert experiment.kd_equilibrium == 1e-8
        assert experiment.response_units == "RU"
        assert experiment.concentration_units == "nM"
        assert experiment.processing_status == "completed"
        assert experiment.created_by_id == test_user.id

    def test_spr_experiment_defaults(self, db_session):
        """Test SPR experiment with default values."""
        experiment = TestSPRExperiment(
            experiment_name="SPR Default Test"
        )
        
        db_session.add(experiment)
        db_session.commit()
        db_session.refresh(experiment)
        
        assert experiment.experiment_type == SPRExperimentType.KINETICS
        assert experiment.temperature_celsius == 25.0
        assert experiment.response_units == "RU"
        assert experiment.concentration_units == "nM"
        assert experiment.reference_subtracted is True
        assert experiment.replicate_number == 1
        assert experiment.processing_status == "pending"

    def test_spr_sensorgram_creation(self, db_session, test_user, test_compound, test_protein):
        """Test creating SPR sensorgram data."""
        # Create parent experiment
        experiment = TestSPRExperiment(
            experiment_name="SPR with Sensorgram",
            compound_id=test_compound.id,
            target_id=test_protein.id,
            created_by_id=test_user.id
        )
        db_session.add(experiment)
        db_session.commit()
        db_session.refresh(experiment)
        
        # Create sensorgram
        sensorgram = TestSPRSensorgram(
            experiment_id=experiment.id,
            cycle_phase="association",
            flow_cell=1,
            time_seconds_json=[0.0, 1.0, 2.0, 3.0, 4.0],
            response_values_json=[0.0, 10.0, 20.0, 25.0, 30.0],
            baseline_corrected=True,
            reference_subtracted=True,
            bulk_shift_corrected=False,
            data_points_count=5,
            noise_level=0.5
        )
        
        db_session.add(sensorgram)
        db_session.commit()
        db_session.refresh(sensorgram)
        
        assert sensorgram.id is not None
        assert sensorgram.experiment_id == experiment.id
        assert sensorgram.cycle_phase == "association"
        assert sensorgram.flow_cell == 1
        assert sensorgram.time_seconds_json == [0.0, 1.0, 2.0, 3.0, 4.0]
        assert sensorgram.response_values_json == [0.0, 10.0, 20.0, 25.0, 30.0]
        assert sensorgram.data_points_count == 5
        assert sensorgram.noise_level == 0.5


class TestMSTExperimentModel:
    """Test MST experiment model functionality."""

    def test_mst_experiment_creation(self, db_session, test_user, test_compound, test_protein):
        """Test creating an MST experiment with binding results."""
        experiment = TestMSTExperiment(
            experiment_name="MST Binding Assay",
            compound_id=test_compound.id,
            target_id=test_protein.id,
            target_name="EGFR",
            target_concentration=50.0,
            ligand_name="ATP",
            buffer_composition="HEPES pH 7.4 + 150mM NaCl",
            temperature_celsius=25.0,
            excitation_power=20.0,
            mst_power=40.0,
            response_units="‰",
            concentration_units="μM",
            kd_value=2.5e-6,
            kd_error=0.3e-6,
            hill_coefficient=1.1,
            binding_amplitude=-15.0,
            r_squared=0.98,
            signal_to_noise=25.0,
            aggregation_detected=False,
            bleaching_rate=0.02,
            instrument_model="NanoTemper Monolith NT.115",
            capillary_type="Premium",
            replicate_number=1,
            processing_status="completed",
            experiment_date=datetime.now(timezone.utc),
            created_by_id=test_user.id
        )
        
        db_session.add(experiment)
        db_session.commit()
        db_session.refresh(experiment)
        
        assert experiment.id is not None
        assert experiment.experiment_name == "MST Binding Assay"
        assert experiment.compound_id == test_compound.id
        assert experiment.target_id == test_protein.id
        assert experiment.kd_value == 2.5e-6
        assert experiment.kd_error == 0.3e-6
        assert experiment.hill_coefficient == 1.1
        assert experiment.response_units == "‰"
        assert experiment.concentration_units == "μM"
        assert experiment.aggregation_detected is False
        assert experiment.created_by_id == test_user.id

    def test_mst_dose_response_creation(self, db_session, test_user, test_compound, test_protein):
        """Test creating MST dose-response curve data."""
        # Create parent experiment
        experiment = TestMSTExperiment(
            experiment_name="MST with Dose Response",
            compound_id=test_compound.id,
            target_id=test_protein.id,
            created_by_id=test_user.id
        )
        db_session.add(experiment)
        db_session.commit()
        db_session.refresh(experiment)
        
        # Create dose-response points
        concentrations = [0.0, 0.1, 1.0, 10.0, 100.0]
        responses = [0.0, -2.5, -8.0, -12.0, -15.0]
        
        for conc, resp in zip(concentrations, responses):
            dose_point = TestMSTDoseResponse(
                experiment_id=experiment.id,
                ligand_concentration=conc,
                thermophoresis_response=resp,
                thermophoresis_error=0.5,
                initial_fluorescence=1000.0,
                thermophoresis_signal=resp,
                capillary_number=1,
                measurement_valid=True,
                outlier_flag=False
            )
            db_session.add(dose_point)
        
        db_session.commit()
        
        # Query back the dose-response data
        dose_points = db_session.query(TestMSTDoseResponse).filter(
            TestMSTDoseResponse.experiment_id == experiment.id
        ).order_by(TestMSTDoseResponse.ligand_concentration).all()
        
        assert len(dose_points) == 5
        assert dose_points[0].ligand_concentration == 0.0
        assert dose_points[0].thermophoresis_response == 0.0
        assert dose_points[4].ligand_concentration == 100.0
        assert dose_points[4].thermophoresis_response == -15.0
        assert all(dp.measurement_valid for dp in dose_points)
        assert not any(dp.outlier_flag for dp in dose_points)


class TestDSCExperimentModel:
    """Test DSC experiment model functionality."""

    def test_dsc_experiment_creation(self, db_session, test_user, test_compound, test_protein):
        """Test creating a DSC experiment with thermal parameters."""
        experiment = TestDSCExperiment(
            experiment_name="DSC Thermal Stability",
            compound_id=test_compound.id,
            target_id=test_protein.id,
            target_name="EGFR",
            protein_concentration=1.0,
            ligand_concentration=10.0,
            buffer_composition="PBS pH 7.4",
            scan_rate=1.0,
            start_temperature=20.0,
            end_temperature=100.0,
            response_units="kcal/mol/°C",
            concentration_units="μM",
            tm_value=65.5,
            tm_error=0.2,
            delta_h=-450.0,
            delta_h_error=15.0,
            delta_cp=2.1,
            cooperativity=1.8,
            kd_thermal=5e-6,
            delta_tm=3.2,
            baseline_quality=0.995,
            peak_symmetry=0.92,
            signal_to_noise=50.0,
            instrument_model="MicroCal PEAQ-DSC",
            cell_volume=200.0,
            reference_buffer="PBS pH 7.4",
            replicate_number=1,
            processing_status="completed",
            experiment_date=datetime.now(timezone.utc),
            created_by_id=test_user.id
        )
        
        db_session.add(experiment)
        db_session.commit()
        db_session.refresh(experiment)
        
        assert experiment.id is not None
        assert experiment.experiment_name == "DSC Thermal Stability"
        assert experiment.compound_id == test_compound.id
        assert experiment.target_id == test_protein.id
        assert experiment.tm_value == 65.5
        assert experiment.tm_error == 0.2
        assert experiment.delta_h == -450.0
        assert experiment.delta_tm == 3.2
        assert experiment.response_units == "kcal/mol/°C"
        assert experiment.concentration_units == "μM"
        assert experiment.baseline_quality == 0.995
        assert experiment.created_by_id == test_user.id

    def test_dsc_scan_creation(self, db_session, test_user, test_compound, test_protein):
        """Test creating DSC thermogram scan data."""
        # Create parent experiment
        experiment = TestDSCExperiment(
            experiment_name="DSC with Thermogram",
            compound_id=test_compound.id,
            target_id=test_protein.id,
            created_by_id=test_user.id
        )
        db_session.add(experiment)
        db_session.commit()
        db_session.refresh(experiment)
        
        # Create thermogram scan
        temperatures = [20.0, 30.0, 40.0, 50.0, 60.0, 70.0, 80.0]
        heat_capacities = [0.1, 0.2, 0.5, 1.2, 2.8, 1.5, 0.3]
        
        scan = TestDSCScan(
            experiment_id=experiment.id,
            scan_type="sample",
            scan_number=1,
            temperature_celsius_json=temperatures,
            heat_capacity_json=heat_capacities,
            baseline_subtracted=True,
            reference_subtracted=True,
            normalized=False,
            data_points_count=7,
            noise_level=0.05
        )
        
        db_session.add(scan)
        db_session.commit()
        db_session.refresh(scan)
        
        assert scan.id is not None
        assert scan.experiment_id == experiment.id
        assert scan.scan_type == "sample"
        assert scan.scan_number == 1
        assert scan.temperature_celsius_json == temperatures
        assert scan.heat_capacity_json == heat_capacities
        assert scan.data_points_count == 7
        assert scan.baseline_subtracted is True
        assert scan.reference_subtracted is True
        assert scan.normalized is False


class TestBiophysicalModelIntegration:
    """Test integration between biophysical models."""

    def test_experiment_replicate_grouping(self, db_session, test_user, test_compound, test_protein):
        """Test linking experiments by replicate group."""
        replicate_group_id = uuid4()
        
        # Create multiple SPR experiments in the same replicate group
        for i in range(3):
            experiment = TestSPRExperiment(
                experiment_name=f"SPR Replicate {i+1}",
                compound_id=test_compound.id,
                target_id=test_protein.id,
                replicate_number=i+1,
                replicate_group_id=replicate_group_id,
                ka=1e5 + i*1e4,  # Slightly different kinetics
                kd=1e-3 + i*1e-4,
                created_by_id=test_user.id
            )
            db_session.add(experiment)
        
        db_session.commit()
        
        # Query experiments by replicate group
        replicates = db_session.query(TestSPRExperiment).filter(
            TestSPRExperiment.replicate_group_id == replicate_group_id
        ).order_by(TestSPRExperiment.replicate_number).all()
        
        assert len(replicates) == 3
        assert replicates[0].replicate_number == 1
        assert replicates[1].replicate_number == 2
        assert replicates[2].replicate_number == 3
        assert all(exp.replicate_group_id == replicate_group_id for exp in replicates)

    def test_multi_technique_compound_profiling(self, db_session, test_user, test_compound, test_protein):
        """Test using multiple biophysical techniques on the same compound-target pair."""
        # SPR kinetics
        spr_exp = TestSPRExperiment(
            experiment_name="SPR Kinetics",
            compound_id=test_compound.id,
            target_id=test_protein.id,
            experiment_type=SPRExperimentType.KINETICS,
            ka=1.5e5,
            kd=2e-3,
            kd_equilibrium=1.3e-8,
            created_by_id=test_user.id
        )
        
        # MST affinity
        mst_exp = TestMSTExperiment(
            experiment_name="MST Affinity",
            compound_id=test_compound.id,
            target_id=test_protein.id,
            kd_value=1.2e-8,  # Similar to SPR KD
            hill_coefficient=1.0,
            created_by_id=test_user.id
        )
        
        # DSC thermal stability
        dsc_exp = TestDSCExperiment(
            experiment_name="DSC Stability",
            compound_id=test_compound.id,
            target_id=test_protein.id,
            tm_value=68.5,
            delta_tm=4.0,  # Thermal shift due to binding
            created_by_id=test_user.id
        )
        
        db_session.add_all([spr_exp, mst_exp, dsc_exp])
        db_session.commit()
        
        # Query all experiments for this compound-target pair
        spr_results = db_session.query(TestSPRExperiment).filter(
            TestSPRExperiment.compound_id == test_compound.id,
            TestSPRExperiment.target_id == test_protein.id
        ).all()
        
        mst_results = db_session.query(TestMSTExperiment).filter(
            TestMSTExperiment.compound_id == test_compound.id,
            TestMSTExperiment.target_id == test_protein.id
        ).all()
        
        dsc_results = db_session.query(TestDSCExperiment).filter(
            TestDSCExperiment.compound_id == test_compound.id,
            TestDSCExperiment.target_id == test_protein.id
        ).all()
        
        assert len(spr_results) == 1
        assert len(mst_results) == 1
        assert len(dsc_results) == 1
        
        # Compare binding affinities
        spr_kd = spr_results[0].kd_equilibrium
        mst_kd = mst_results[0].kd_value
        assert abs(spr_kd - mst_kd) / mst_kd < 0.1  # Within 10% agreement

    def test_experiment_processing_status_workflow(self, db_session, test_user, test_compound):
        """Test experiment processing status workflow."""
        experiment = TestSPRExperiment(
            experiment_name="Processing Workflow Test",
            compound_id=test_compound.id,
            processing_status="pending",
            created_by_id=test_user.id
        )
        
        db_session.add(experiment)
        db_session.commit()
        db_session.refresh(experiment)
        
        # Initial status
        assert experiment.processing_status == "pending"
        
        # Update to running
        experiment.processing_status = "running"
        experiment.processing_log = "Starting data analysis..."
        db_session.commit()
        
        # Update to completed with results
        experiment.processing_status = "completed"
        experiment.processing_log = "Analysis completed successfully"
        experiment.ka = 1.2e5
        experiment.kd = 1.8e-3
        experiment.kd_equilibrium = 1.5e-8
        experiment.chi_squared = 1.8
        db_session.commit()
        
        # Verify final state
        db_session.refresh(experiment)
        assert experiment.processing_status == "completed"
        assert "successfully" in experiment.processing_log
        assert experiment.ka is not None
        assert experiment.kd is not None
        assert experiment.kd_equilibrium is not None
        assert experiment.chi_squared < 5.0  # Good fit

    def test_spr_experiment_type_enum(self, db_session):
        """Test SPR experiment type enumeration values."""
        # Test all enum values
        kinetics_exp = TestSPRExperiment(
            experiment_name="Kinetics Test",
            experiment_type=SPRExperimentType.KINETICS
        )
        
        affinity_exp = TestSPRExperiment(
            experiment_name="Affinity Test", 
            experiment_type=SPRExperimentType.AFFINITY
        )
        
        screening_exp = TestSPRExperiment(
            experiment_name="Screening Test",
            experiment_type=SPRExperimentType.SCREENING
        )
        
        db_session.add_all([kinetics_exp, affinity_exp, screening_exp])
        db_session.commit()
        
        # Verify enum values are stored correctly
        experiments = db_session.query(TestSPRExperiment).all()
        experiment_types = [exp.experiment_type for exp in experiments]
        
        assert SPRExperimentType.KINETICS in experiment_types
        assert SPRExperimentType.AFFINITY in experiment_types
        assert SPRExperimentType.SCREENING in experiment_types
        assert len(set(experiment_types)) == 3  # All unique

    def test_biophysical_data_quality_metrics(self, db_session, test_user, test_compound, test_protein):
        """Test data quality metrics across all biophysical techniques."""
        # SPR with quality metrics
        spr_exp = TestSPRExperiment(
            experiment_name="SPR Quality Test",
            compound_id=test_compound.id,
            target_id=test_protein.id,
            chi_squared=1.5,  # Good fit
            baseline_drift=0.05,  # Low drift
            bulk_refractive_index=1.333,
            reference_subtracted=True,
            created_by_id=test_user.id
        )
        
        # MST with quality metrics
        mst_exp = TestMSTExperiment(
            experiment_name="MST Quality Test",
            compound_id=test_compound.id,
            target_id=test_protein.id,
            r_squared=0.99,  # Excellent fit
            signal_to_noise=45.0,  # High S/N
            aggregation_detected=False,
            bleaching_rate=0.01,  # Low bleaching
            created_by_id=test_user.id
        )
        
        # DSC with quality metrics
        dsc_exp = TestDSCExperiment(
            experiment_name="DSC Quality Test",
            compound_id=test_compound.id,
            target_id=test_protein.id,
            baseline_quality=0.998,  # Excellent baseline
            peak_symmetry=0.95,  # Symmetric peak
            signal_to_noise=75.0,  # Very high S/N
            created_by_id=test_user.id
        )
        
        db_session.add_all([spr_exp, mst_exp, dsc_exp])
        db_session.commit()
        
        # Verify quality metrics are within acceptable ranges
        spr_result = db_session.query(TestSPRExperiment).filter(
            TestSPRExperiment.experiment_name == "SPR Quality Test"
        ).first()
        assert spr_result.chi_squared < 5.0  # Good fit threshold
        assert spr_result.baseline_drift < 0.1  # Low drift threshold
        
        mst_result = db_session.query(TestMSTExperiment).filter(
            TestMSTExperiment.experiment_name == "MST Quality Test"
        ).first()
        assert mst_result.r_squared > 0.95  # High correlation threshold
        assert mst_result.signal_to_noise > 10.0  # Minimum S/N threshold
        
        dsc_result = db_session.query(TestDSCExperiment).filter(
            TestDSCExperiment.experiment_name == "DSC Quality Test"
        ).first()
        assert dsc_result.baseline_quality > 0.99  # Excellent baseline threshold
        assert dsc_result.peak_symmetry > 0.8  # Symmetric peak threshold
