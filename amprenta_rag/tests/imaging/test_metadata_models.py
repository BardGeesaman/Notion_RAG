"""Unit tests for imaging metadata models."""

import pytest
from datetime import datetime, timezone
from uuid import uuid4

from sqlalchemy import create_engine
from sqlalchemy.orm import Session, sessionmaker

from amprenta_rag.imaging.models_metadata import (
    Microscope,
    Objective,
    LightSource,
    FilterSet,
    ChannelConfig,
    AcquisitionSettings,
    ImageFileSet,
)


@pytest.fixture
def db():
    """Create in-memory SQLite database for testing."""
    engine = create_engine("sqlite:///:memory:", echo=False)
    
    # Create all metadata tables
    Microscope.__table__.create(engine)
    Objective.__table__.create(engine)
    LightSource.__table__.create(engine)
    FilterSet.__table__.create(engine)
    ChannelConfig.__table__.create(engine)
    AcquisitionSettings.__table__.create(engine)
    ImageFileSet.__table__.create(engine)
    
    Session = sessionmaker(bind=engine)
    session = Session()
    yield session
    session.close()


class TestMicroscope:
    """Test Microscope model."""

    def test_create_microscope(self, db: Session):
        """Test creating a microscope with basic fields."""
        microscope = Microscope(
            id=uuid4(),
            name="Nikon Ti2-E",
            manufacturer="Nikon",
            model="Ti2-E",
            serial_number="NK123456",
            facility_location="Building A, Room 101",
            is_active=True
        )
        db.add(microscope)
        db.commit()
        
        assert microscope.id is not None
        assert microscope.name == "Nikon Ti2-E"
        assert microscope.manufacturer == "Nikon"
        assert microscope.model == "Ti2-E"
        assert microscope.serial_number == "NK123456"
        assert microscope.facility_location == "Building A, Room 101"
        assert microscope.is_active is True
        assert microscope.created_at is not None

    def test_microscope_with_relationships(self, db: Session):
        """Test microscope with objectives and light sources."""
        microscope = Microscope(
            id=uuid4(),
            name="Zeiss Axio Observer",
            manufacturer="Zeiss",
            model="Axio Observer",
            is_active=True
        )
        db.add(microscope)
        db.commit()
        
        # Add objective
        objective = Objective(
            id=uuid4(),
            microscope_id=microscope.id,
            name="Plan Apo 20x",
            magnification=20.0,
            numerical_aperture=0.75,
            immersion="air"
        )
        db.add(objective)
        
        # Add light source
        light_source = LightSource(
            id=uuid4(),
            microscope_id=microscope.id,
            name="488nm Laser",
            source_type="laser",
            wavelength_nm=488,
            max_power_mw=50.0
        )
        db.add(light_source)
        db.commit()
        
        # Verify relationships
        db.refresh(microscope)
        assert len(microscope.objectives) == 1
        assert len(microscope.light_sources) == 1
        assert microscope.objectives[0].name == "Plan Apo 20x"
        assert microscope.light_sources[0].wavelength_nm == 488


class TestObjective:
    """Test Objective model."""

    def test_create_objective(self, db: Session):
        """Test creating an objective lens."""
        objective = Objective(
            id=uuid4(),
            name="Plan Fluor 40x",
            magnification=40.0,
            numerical_aperture=1.3,
            immersion="oil",
            working_distance_mm=0.21,
            correction="Plan Fluor",
            is_active=True
        )
        db.add(objective)
        db.commit()
        
        assert objective.id is not None
        assert objective.name == "Plan Fluor 40x"
        assert objective.magnification == 40.0
        assert objective.numerical_aperture == 1.3
        assert objective.immersion == "oil"
        assert objective.working_distance_mm == 0.21
        assert objective.correction == "Plan Fluor"
        assert objective.is_active is True

    def test_objective_without_microscope(self, db: Session):
        """Test creating a shared objective not tied to specific microscope."""
        objective = Objective(
            id=uuid4(),
            microscope_id=None,  # Shared objective
            name="Universal 10x",
            magnification=10.0,
            numerical_aperture=0.25,
            immersion="air"
        )
        db.add(objective)
        db.commit()
        
        assert objective.microscope_id is None
        assert objective.microscope is None


class TestLightSource:
    """Test LightSource model."""

    def test_create_light_source(self, db: Session):
        """Test creating a light source."""
        microscope = Microscope(
            id=uuid4(),
            name="Test Microscope",
            manufacturer="Test",
            model="Test Model"
        )
        db.add(microscope)
        db.commit()
        
        light_source = LightSource(
            id=uuid4(),
            microscope_id=microscope.id,
            name="561nm DPSS Laser",
            source_type="laser",
            wavelength_nm=561,
            max_power_mw=100.0,
            is_active=True
        )
        db.add(light_source)
        db.commit()
        
        assert light_source.id is not None
        assert light_source.microscope_id == microscope.id
        assert light_source.name == "561nm DPSS Laser"
        assert light_source.source_type == "laser"
        assert light_source.wavelength_nm == 561
        assert light_source.max_power_mw == 100.0
        assert light_source.is_active is True

    def test_led_light_source(self, db: Session):
        """Test creating an LED light source."""
        microscope = Microscope(
            id=uuid4(),
            name="Test Microscope",
            manufacturer="Test",
            model="Test Model"
        )
        db.add(microscope)
        db.commit()
        
        led_source = LightSource(
            id=uuid4(),
            microscope_id=microscope.id,
            name="White LED",
            source_type="led",
            wavelength_nm=None,  # Broadband
            max_power_mw=None,
            is_active=True
        )
        db.add(led_source)
        db.commit()
        
        assert led_source.source_type == "led"
        assert led_source.wavelength_nm is None


class TestFilterSet:
    """Test FilterSet model."""

    def test_create_filter_set(self, db: Session):
        """Test creating a filter set."""
        filter_set = FilterSet(
            id=uuid4(),
            name="FITC",
            excitation_center_nm=485,
            excitation_bandwidth_nm=20,
            emission_center_nm=528,
            emission_bandwidth_nm=38,
            dichroic_cutoff_nm=505
        )
        db.add(filter_set)
        db.commit()
        
        assert filter_set.id is not None
        assert filter_set.name == "FITC"
        assert filter_set.excitation_center_nm == 485
        assert filter_set.excitation_bandwidth_nm == 20
        assert filter_set.emission_center_nm == 528
        assert filter_set.emission_bandwidth_nm == 38
        assert filter_set.dichroic_cutoff_nm == 505

    def test_filter_set_without_dichroic(self, db: Session):
        """Test creating filter set without dichroic mirror."""
        filter_set = FilterSet(
            id=uuid4(),
            name="DAPI",
            excitation_center_nm=358,
            excitation_bandwidth_nm=15,
            emission_center_nm=461,
            emission_bandwidth_nm=25,
            dichroic_cutoff_nm=None
        )
        db.add(filter_set)
        db.commit()
        
        assert filter_set.dichroic_cutoff_nm is None


class TestChannelConfig:
    """Test ChannelConfig model."""

    def test_create_channel_config(self, db: Session):
        """Test creating a channel configuration."""
        microscope = Microscope(
            id=uuid4(),
            name="Test Microscope",
            manufacturer="Test",
            model="Test Model"
        )
        db.add(microscope)
        
        light_source = LightSource(
            id=uuid4(),
            microscope_id=microscope.id,
            name="488nm Laser",
            source_type="laser",
            wavelength_nm=488
        )
        db.add(light_source)
        
        filter_set = FilterSet(
            id=uuid4(),
            name="GFP",
            excitation_center_nm=485,
            excitation_bandwidth_nm=20,
            emission_center_nm=528,
            emission_bandwidth_nm=38
        )
        db.add(filter_set)
        db.commit()
        
        channel_config = ChannelConfig(
            id=uuid4(),
            microscope_id=microscope.id,
            channel_name="GFP",
            fluorophore="Enhanced GFP",
            light_source_id=light_source.id,
            filter_set_id=filter_set.id,
            default_exposure_ms=100.0,
            default_gain=1.5
        )
        db.add(channel_config)
        db.commit()
        
        assert channel_config.id is not None
        assert channel_config.microscope_id == microscope.id
        assert channel_config.channel_name == "GFP"
        assert channel_config.fluorophore == "Enhanced GFP"
        assert channel_config.light_source_id == light_source.id
        assert channel_config.filter_set_id == filter_set.id
        assert channel_config.default_exposure_ms == 100.0
        assert channel_config.default_gain == 1.5

    def test_channel_config_minimal(self, db: Session):
        """Test creating minimal channel config."""
        microscope = Microscope(
            id=uuid4(),
            name="Test Microscope",
            manufacturer="Test",
            model="Test Model"
        )
        db.add(microscope)
        db.commit()
        
        channel_config = ChannelConfig(
            id=uuid4(),
            microscope_id=microscope.id,
            channel_name="Brightfield",
            fluorophore=None,
            light_source_id=None,
            filter_set_id=None,
            default_exposure_ms=None,
            default_gain=None
        )
        db.add(channel_config)
        db.commit()
        
        assert channel_config.channel_name == "Brightfield"
        assert channel_config.fluorophore is None


class TestAcquisitionSettings:
    """Test AcquisitionSettings model."""

    def test_create_acquisition_settings(self, db: Session):
        """Test creating acquisition settings."""
        settings = AcquisitionSettings(
            id=uuid4(),
            exposure_ms=50.0,
            gain=2.0,
            laser_power_percent=10.0,
            binning=1,
            z_position_um=0.0,
            autofocus_offset_um=0.5,
            temperature_celsius=37.0,
            raw_settings={"vendor_param1": "value1", "vendor_param2": 42}
        )
        db.add(settings)
        db.commit()
        
        assert settings.id is not None
        assert settings.exposure_ms == 50.0
        assert settings.gain == 2.0
        assert settings.laser_power_percent == 10.0
        assert settings.binning == 1
        assert settings.z_position_um == 0.0
        assert settings.autofocus_offset_um == 0.5
        assert settings.temperature_celsius == 37.0
        assert settings.raw_settings == {"vendor_param1": "value1", "vendor_param2": 42}

    def test_acquisition_settings_minimal(self, db: Session):
        """Test creating minimal acquisition settings."""
        settings = AcquisitionSettings(
            id=uuid4(),
            exposure_ms=100.0,  # Required field
            gain=None,
            laser_power_percent=None,
            binning=2,  # Default 1, but testing override
            z_position_um=None,
            autofocus_offset_um=None,
            temperature_celsius=None,
            raw_settings=None
        )
        db.add(settings)
        db.commit()
        
        assert settings.exposure_ms == 100.0
        assert settings.binning == 2
        assert settings.gain is None


class TestImageFileSet:
    """Test ImageFileSet model."""

    def test_create_image_file_set(self, db: Session):
        """Test creating an image file set."""
        file_set = ImageFileSet(
            id=uuid4(),
            plate_id=None,  # Standalone import
            vendor="opera",
            import_path="/data/opera/experiment_001",
            file_count=96,
            image_count=384,
            import_status="completed",
            error_message=None,
            completed_at=datetime.now(timezone.utc)
        )
        db.add(file_set)
        db.commit()
        
        assert file_set.id is not None
        assert file_set.plate_id is None
        assert file_set.vendor == "opera"
        assert file_set.import_path == "/data/opera/experiment_001"
        assert file_set.file_count == 96
        assert file_set.image_count == 384
        assert file_set.import_status == "completed"
        assert file_set.error_message is None
        assert file_set.completed_at is not None

    def test_image_file_set_with_error(self, db: Session):
        """Test creating file set with error status."""
        file_set = ImageFileSet(
            id=uuid4(),
            vendor="imagexpress",
            import_path="/data/imagexpress/failed_run",
            file_count=48,
            image_count=0,
            import_status="failed",
            error_message="File format not recognized",
            completed_at=None
        )
        db.add(file_set)
        db.commit()
        
        assert file_set.import_status == "failed"
        assert file_set.error_message == "File format not recognized"
        assert file_set.completed_at is None

    def test_image_file_set_defaults(self, db: Session):
        """Test file set with default values."""
        file_set = ImageFileSet(
            id=uuid4(),
            vendor="cellvoyager",
            import_path="/data/cellvoyager/new_experiment",
            file_count=192
            # image_count defaults to 0
            # import_status defaults to "pending"
        )
        db.add(file_set)
        db.commit()
        
        assert file_set.image_count == 0
        assert file_set.import_status == "pending"
        assert file_set.error_message is None
        assert file_set.completed_at is None
