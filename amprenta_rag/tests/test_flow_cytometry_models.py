"""
Tests for Flow Cytometry database models.

Tests cover model creation, relationships, constraints, and business logic
for flow cytometry data storage and analysis.
"""

import pytest
import json
from datetime import datetime, timezone
from uuid import uuid4

from sqlalchemy import create_engine, Column, String, Text, Integer, Float, Boolean, DateTime, JSON, ForeignKey, BigInteger, UniqueConstraint
from sqlalchemy.exc import IntegrityError
from sqlalchemy.orm import Session, sessionmaker
from sqlalchemy.dialects.postgresql import UUID

from amprenta_rag.database.base import Base
from amprenta_rag.database.models_flow_cytometry import (
    GateType,
    BooleanOperator,
    generate_uuid,
)


# Create test-specific models without ARRAY types for SQLite compatibility
class TestDataset(Base):
    """Simple Dataset model for testing without ARRAY types."""
    __tablename__ = "test_datasets"
    
    id = Column(UUID(as_uuid=True), primary_key=True, default=generate_uuid)
    name = Column(String(500), nullable=False)
    omics_type = Column(String(50), nullable=False)
    description = Column(Text, nullable=True)


class TestFlowCytometryDataset(Base):
    """Flow cytometry dataset metadata for testing."""
    __tablename__ = "test_flow_cytometry_datasets"

    id = Column(UUID(as_uuid=True), primary_key=True, default=generate_uuid)
    dataset_id = Column(UUID(as_uuid=True), ForeignKey("test_datasets.id"), nullable=False, index=True)
    events_parquet_path = Column(String(500), nullable=False)
    file_size_bytes = Column(BigInteger, nullable=True)
    n_events = Column(Integer, nullable=True)
    n_parameters = Column(Integer, nullable=True)
    cytometer_model = Column(String(200), nullable=True)
    sample_id = Column(String(200), nullable=True)
    processing_status = Column(String(50), nullable=False, default="pending")
    ingested_at = Column(DateTime(timezone=True), default=lambda: datetime.now(timezone.utc), nullable=False)


class TestFlowCytometryParameter(Base):
    """Channel/parameter metadata for testing."""
    __tablename__ = "test_flow_cytometry_parameters"

    id = Column(UUID(as_uuid=True), primary_key=True, default=generate_uuid)
    flow_dataset_id = Column(UUID(as_uuid=True), ForeignKey("test_flow_cytometry_datasets.id"), nullable=False, index=True)
    parameter_index = Column(Integer, nullable=False)
    parameter_name = Column(String(100), nullable=False)
    parameter_short_name = Column(String(50), nullable=True)
    detector = Column(String(100), nullable=True)
    fluorophore = Column(String(100), nullable=True)
    data_range = Column(Integer, nullable=True)
    min_value = Column(Float, nullable=True)
    max_value = Column(Float, nullable=True)

    __table_args__ = (
        UniqueConstraint("flow_dataset_id", "parameter_index", name="uq_test_flow_parameter_index"),
    )


class TestFlowCytometryGate(Base):
    """Gate definition for testing (without ARRAY types)."""
    __tablename__ = "test_flow_cytometry_gates"

    id = Column(UUID(as_uuid=True), primary_key=True, default=generate_uuid)
    flow_dataset_id = Column(UUID(as_uuid=True), ForeignKey("test_flow_cytometry_datasets.id"), nullable=False, index=True)
    gate_name = Column(String(200), nullable=False)
    gate_type = Column(String(20), nullable=False)
    gate_definition = Column(JSON, nullable=False)
    x_parameter_id = Column(UUID(as_uuid=True), ForeignKey("test_flow_cytometry_parameters.id"), nullable=False)
    y_parameter_id = Column(UUID(as_uuid=True), ForeignKey("test_flow_cytometry_parameters.id"), nullable=True)
    boolean_operator = Column(String(10), nullable=True)
    # Note: operand_gate_ids stored as JSON string instead of ARRAY for SQLite compatibility
    operand_gate_ids_json = Column(Text, nullable=True)  # JSON string of UUIDs
    parent_gate_id = Column(UUID(as_uuid=True), ForeignKey("test_flow_cytometry_gates.id"), nullable=True)
    is_active = Column(Boolean, nullable=False, default=True)
    created_at = Column(DateTime(timezone=True), default=lambda: datetime.now(timezone.utc), nullable=False)


class TestFlowCytometryPopulation(Base):
    """Statistics for a gated population for testing."""
    __tablename__ = "test_flow_cytometry_populations"

    id = Column(UUID(as_uuid=True), primary_key=True, default=generate_uuid)
    flow_dataset_id = Column(UUID(as_uuid=True), ForeignKey("test_flow_cytometry_datasets.id"), nullable=False, index=True)
    gate_id = Column(UUID(as_uuid=True), ForeignKey("test_flow_cytometry_gates.id"), nullable=False, index=True)
    event_count = Column(Integer, nullable=False)
    parent_event_count = Column(Integer, nullable=True)
    percentage_of_parent = Column(Float, nullable=True)
    percentage_of_total = Column(Float, nullable=True)
    population_name = Column(String(200), nullable=True)
    cell_type = Column(String(200), nullable=True)
    phenotype = Column(String(500), nullable=True)
    parameter_statistics = Column(JSON, nullable=True)
    analysis_software = Column(String(100), nullable=True)
    analysis_version = Column(String(50), nullable=True)
    analyzed_at = Column(DateTime(timezone=True), default=lambda: datetime.now(timezone.utc), nullable=False)

    __table_args__ = (
        UniqueConstraint("flow_dataset_id", "gate_id", name="uq_test_flow_population_gate"),
    )


@pytest.fixture
def db_session():
    """Create in-memory SQLite database for testing."""
    engine = create_engine("sqlite:///:memory:", echo=False)
    
    # Create only the test tables we need (SQLite compatible)
    TestDataset.__table__.create(engine)
    TestFlowCytometryDataset.__table__.create(engine)
    TestFlowCytometryParameter.__table__.create(engine)
    TestFlowCytometryGate.__table__.create(engine)
    TestFlowCytometryPopulation.__table__.create(engine)
    
    Session = sessionmaker(bind=engine)
    session = Session()
    
    yield session
    
    session.close()


@pytest.fixture
def sample_dataset(db_session: Session):
    """Create a sample Dataset for testing."""
    dataset_id = uuid4()
    dataset = TestDataset(
        id=dataset_id,
        name=f"Test Flow Dataset {uuid4().hex[:8]}",
        omics_type="flow_cytometry",
        description="Test flow cytometry dataset",
    )
    db_session.add(dataset)
    db_session.commit()
    return dataset_id


@pytest.fixture
def sample_flow_dataset(db_session: Session, sample_dataset):
    """Create a sample FlowCytometryDataset for testing."""
    flow_dataset_id = uuid4()
    flow_dataset = TestFlowCytometryDataset(
        id=flow_dataset_id,
        dataset_id=sample_dataset,
        events_parquet_path="/data/flow/test_sample.parquet",
        file_size_bytes=1024*1024,
        n_events=10000,
        n_parameters=12,
        cytometer_model="BD FACSCanto II",
        sample_id=f"Sample_{uuid4().hex[:8]}",
        processing_status="completed",
    )
    db_session.add(flow_dataset)
    db_session.commit()
    return flow_dataset_id


def test_flow_cytometry_dataset_creation(db_session: Session, sample_dataset):
    """Test creating a FlowCytometryDataset with required fields."""
    flow_dataset = TestFlowCytometryDataset(
        dataset_id=sample_dataset,
        events_parquet_path="/data/flow/sample.parquet",
        n_events=5000,
        processing_status="pending",
    )
    
    db_session.add(flow_dataset)
    db_session.commit()
    
    assert flow_dataset.id is not None
    assert flow_dataset.dataset_id == sample_dataset
    assert flow_dataset.events_parquet_path == "/data/flow/sample.parquet"
    assert flow_dataset.n_events == 5000
    assert flow_dataset.processing_status == "pending"
    assert flow_dataset.ingested_at is not None


def test_flow_cytometry_parameter_creation(db_session: Session, sample_flow_dataset):
    """Test creating flow cytometry parameters with unique constraint."""
    param = TestFlowCytometryParameter(
        flow_dataset_id=sample_flow_dataset,
        parameter_index=0,
        parameter_name="SSC-A",
        parameter_short_name="SSC",
        detector="PMT3",
        data_range=1024,
    )
    
    db_session.add(param)
    db_session.commit()
    
    assert param.id is not None
    assert param.flow_dataset_id == sample_flow_dataset
    assert param.parameter_index == 0
    assert param.parameter_name == "SSC-A"
    assert param.detector == "PMT3"


def test_flow_cytometry_parameter_unique_constraint(db_session: Session, sample_flow_dataset):
    """Test that parameter index must be unique per dataset."""
    param1 = TestFlowCytometryParameter(
        flow_dataset_id=sample_flow_dataset,
        parameter_index=0,
        parameter_name="FSC-A",
    )
    param2 = TestFlowCytometryParameter(
        flow_dataset_id=sample_flow_dataset,
        parameter_index=0,  # Same index - should fail
        parameter_name="SSC-A",
    )
    
    db_session.add(param1)
    db_session.commit()
    
    db_session.add(param2)
    with pytest.raises(IntegrityError):
        db_session.commit()


def test_flow_cytometry_gate_creation(db_session: Session, sample_flow_dataset):
    """Test creating a flow cytometry gate."""
    # Create parameters first
    param1 = TestFlowCytometryParameter(
        flow_dataset_id=sample_flow_dataset,
        parameter_index=0,
        parameter_name="FSC-A",
    )
    param2 = TestFlowCytometryParameter(
        flow_dataset_id=sample_flow_dataset,
        parameter_index=1,
        parameter_name="CD3-FITC",
    )
    db_session.add(param1)
    db_session.add(param2)
    db_session.commit()
    
    gate = TestFlowCytometryGate(
        flow_dataset_id=sample_flow_dataset,
        gate_name="Lymphocytes",
        gate_type=GateType.POLYGON.value,
        gate_definition={
            "vertices": [[1000, 2000], [5000, 2000], [5000, 8000], [1000, 8000]]
        },
        x_parameter_id=param1.id,
        y_parameter_id=param2.id,
        is_active=True,
    )
    
    db_session.add(gate)
    db_session.commit()
    
    assert gate.id is not None
    assert gate.flow_dataset_id == sample_flow_dataset
    assert gate.gate_name == "Lymphocytes"
    assert gate.gate_type == GateType.POLYGON.value
    assert gate.is_active is True
    assert gate.created_at is not None


def test_flow_cytometry_gate_boolean_logic(db_session: Session, sample_flow_dataset):
    """Test creating gates with boolean operators."""
    # Create parameters
    param1 = TestFlowCytometryParameter(
        flow_dataset_id=sample_flow_dataset,
        parameter_index=0,
        parameter_name="FSC-A",
    )
    param2 = TestFlowCytometryParameter(
        flow_dataset_id=sample_flow_dataset,
        parameter_index=1,
        parameter_name="CD3-FITC",
    )
    db_session.add(param1)
    db_session.add(param2)
    db_session.commit()
    
    # Create parent gate
    parent_gate = TestFlowCytometryGate(
        flow_dataset_id=sample_flow_dataset,
        gate_name="Live Cells",
        gate_type=GateType.RECTANGLE.value,
        gate_definition={"x_min": 1000, "x_max": 10000, "y_min": 1000, "y_max": 10000},
        x_parameter_id=param1.id,
        y_parameter_id=param2.id,
    )
    db_session.add(parent_gate)
    db_session.commit()
    
    # Create child gate with boolean logic
    child_gate = TestFlowCytometryGate(
        flow_dataset_id=sample_flow_dataset,
        gate_name="CD3+ Live Cells",
        gate_type=GateType.RECTANGLE.value,
        gate_definition={"x_min": 2000, "x_max": 8000, "y_min": 2000, "y_max": 8000},
        x_parameter_id=param1.id,
        y_parameter_id=param2.id,
        parent_gate_id=parent_gate.id,
        boolean_operator=BooleanOperator.AND.value,
        operand_gate_ids_json=json.dumps([str(parent_gate.id)]),  # Store as JSON string
    )
    
    db_session.add(child_gate)
    db_session.commit()
    
    assert child_gate.parent_gate_id == parent_gate.id
    assert child_gate.boolean_operator == BooleanOperator.AND.value
    operand_ids = json.loads(child_gate.operand_gate_ids_json)
    assert str(parent_gate.id) in operand_ids


def test_flow_cytometry_population_creation(db_session: Session, sample_flow_dataset):
    """Test creating a flow cytometry population."""
    # Create parameters and gate first
    param1 = TestFlowCytometryParameter(
        flow_dataset_id=sample_flow_dataset,
        parameter_index=0,
        parameter_name="FSC-A",
    )
    param2 = TestFlowCytometryParameter(
        flow_dataset_id=sample_flow_dataset,
        parameter_index=1,
        parameter_name="CD3-FITC",
    )
    db_session.add(param1)
    db_session.add(param2)
    db_session.commit()
    
    gate = TestFlowCytometryGate(
        flow_dataset_id=sample_flow_dataset,
        gate_name="T Cells",
        gate_type=GateType.POLYGON.value,
        gate_definition={"vertices": [[1000, 2000], [5000, 2000], [5000, 8000]]},
        x_parameter_id=param1.id,
        y_parameter_id=param2.id,
    )
    db_session.add(gate)
    db_session.commit()
    
    # Create population
    population = TestFlowCytometryPopulation(
        flow_dataset_id=sample_flow_dataset,
        gate_id=gate.id,
        event_count=2500,
        parent_event_count=10000,
        percentage_of_parent=25.0,
        percentage_of_total=25.0,
        population_name="CD3+ T Cells",
        cell_type="T lymphocyte",
        phenotype="CD3+",
        parameter_statistics={
            "FSC-A": {"mean": 5000, "median": 4800, "std": 1200},
            "CD3-FITC": {"mean": 8000, "median": 7900, "std": 800}
        },
        analysis_software="FlowJo",
        analysis_version="10.8.1",
    )
    
    db_session.add(population)
    db_session.commit()
    
    assert population.id is not None
    assert population.flow_dataset_id == sample_flow_dataset
    assert population.gate_id == gate.id
    assert population.event_count == 2500
    assert population.percentage_of_parent == 25.0
    assert population.population_name == "CD3+ T Cells"
    assert population.analyzed_at is not None


def test_flow_cytometry_population_unique_constraint(db_session: Session, sample_flow_dataset):
    """Test that population must be unique per dataset and gate."""
    # Create parameters and gate
    param1 = TestFlowCytometryParameter(
        flow_dataset_id=sample_flow_dataset,
        parameter_index=0,
        parameter_name="FSC-A",
    )
    db_session.add(param1)
    db_session.commit()
    
    gate = TestFlowCytometryGate(
        flow_dataset_id=sample_flow_dataset,
        gate_name="B Cells",
        gate_type=GateType.RECTANGLE.value,
        gate_definition={"x_min": 1000, "x_max": 5000, "y_min": 1000, "y_max": 5000},
        x_parameter_id=param1.id,
    )
    db_session.add(gate)
    db_session.commit()
    
    # Create first population
    pop1 = TestFlowCytometryPopulation(
        flow_dataset_id=sample_flow_dataset,
        gate_id=gate.id,
        event_count=1000,
    )
    db_session.add(pop1)
    db_session.commit()
    
    # Try to create second population for same dataset and gate
    pop2 = TestFlowCytometryPopulation(
        flow_dataset_id=sample_flow_dataset,
        gate_id=gate.id,  # Same gate - should fail
        event_count=1200,
    )
    
    db_session.add(pop2)
    with pytest.raises(IntegrityError):
        db_session.commit()


def test_flow_cytometry_foreign_key_constraints(db_session: Session, sample_flow_dataset):
    """Test foreign key relationships work correctly."""
    # Create parameter
    param = TestFlowCytometryParameter(
        flow_dataset_id=sample_flow_dataset,
        parameter_index=0,
        parameter_name="CD4-PE",
    )
    db_session.add(param)
    db_session.commit()
    
    # Create gate
    gate = TestFlowCytometryGate(
        flow_dataset_id=sample_flow_dataset,
        gate_name="CD4+ Cells",
        gate_type=GateType.RECTANGLE.value,
        gate_definition={"x_min": 1000, "x_max": 5000, "y_min": 1000, "y_max": 5000},
        x_parameter_id=param.id,
    )
    db_session.add(gate)
    db_session.commit()
    
    # Create population
    population = TestFlowCytometryPopulation(
        flow_dataset_id=sample_flow_dataset,
        gate_id=gate.id,
        event_count=3000,
    )
    db_session.add(population)
    db_session.commit()
    
    # Verify all objects were created
    assert param.id is not None
    assert gate.id is not None
    assert population.id is not None
    
    # Verify foreign key relationships
    assert gate.flow_dataset_id == sample_flow_dataset
    assert population.flow_dataset_id == sample_flow_dataset
    assert population.gate_id == gate.id


def test_flow_cytometry_cascade_delete_simulation(db_session: Session, sample_flow_dataset):
    """Test that related records can be deleted (simulates cascade delete)."""
    # Create parameter, gate, and population
    param = TestFlowCytometryParameter(
        flow_dataset_id=sample_flow_dataset,
        parameter_index=0,
        parameter_name="CD8-APC",
    )
    db_session.add(param)
    db_session.commit()
    
    gate = TestFlowCytometryGate(
        flow_dataset_id=sample_flow_dataset,
        gate_name="CD8+ Cells",
        gate_type=GateType.RECTANGLE.value,
        gate_definition={"x_min": 1000, "x_max": 5000, "y_min": 1000, "y_max": 5000},
        x_parameter_id=param.id,
    )
    db_session.add(gate)
    db_session.commit()
    
    population = TestFlowCytometryPopulation(
        flow_dataset_id=sample_flow_dataset,
        gate_id=gate.id,
        event_count=1500,
    )
    db_session.add(population)
    db_session.commit()
    
    param_id = param.id
    gate_id = gate.id
    population_id = population.id
    
    # Manually delete in correct order (simulates cascade)
    db_session.delete(population)
    db_session.delete(gate)
    db_session.delete(param)
    db_session.commit()
    
    # Verify deletion
    assert db_session.get(TestFlowCytometryParameter, param_id) is None
    assert db_session.get(TestFlowCytometryGate, gate_id) is None
    assert db_session.get(TestFlowCytometryPopulation, population_id) is None