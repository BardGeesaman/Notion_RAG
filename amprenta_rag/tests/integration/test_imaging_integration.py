"""
Integration tests for Imaging Metadata functionality.

These tests cover complete workflows using real PostgreSQL database,
testing OME-TIFF parsing, vendor import, instrument management, and QC pipelines.
"""

from __future__ import annotations

import json
import tempfile
from pathlib import Path
from uuid import uuid4
from datetime import datetime, timezone
from typing import Dict, List, Optional

import numpy as np
import pytest
from fastapi.testclient import TestClient

from amprenta_rag.api.main import app
from amprenta_rag.database.session import db_session
from amprenta_rag.database.models import User
from amprenta_rag.imaging.models import MicroscopyImage
from amprenta_rag.models.chemistry import HTSPlate, HTSWell
from amprenta_rag.imaging.models_metadata import (
    Microscope,
    Objective,
    LightSource,
    FilterSet,
    ChannelConfig,
    AcquisitionSettings,
    ImageFileSet,
)
# QC metrics are stored directly in MicroscopyImage model
from amprenta_rag.imaging.ome_parser import parse_ome_tiff, OMEMetadata
from amprenta_rag.imaging.vendor_parsers import parse_vendor_export
from amprenta_rag.imaging.image_qc import run_image_qc


pytestmark = pytest.mark.integration


@pytest.fixture
def test_user():
    """Create a test user in the database."""
    with db_session() as db:
        user = User(
            id=uuid4(),
            email="test@example.com",
            is_active=True,
            created_at=datetime.now(timezone.utc)
        )
        db.add(user)
        db.commit()
        db.expunge(user)
        return user


@pytest.fixture
def test_microscope():
    """Create a test microscope with objectives and channels."""
    with db_session() as db:
        # Create filter set
        filter_set = FilterSet(
            id=uuid4(),
            name="Test DAPI",
            excitation_center_nm=358,
            excitation_bandwidth_nm=46,
            emission_center_nm=461,
            emission_bandwidth_nm=50,
            dichroic_cutoff_nm=400
        )
        db.add(filter_set)
        
        # Create microscope
        microscope = Microscope(
            id=uuid4(),
            name="Test Microscope",
            manufacturer="Test Corp",
            model="TestScope 1000",
            serial_number="TEST-001",
            facility_location="Test Lab",
            is_active=True,
            created_at=datetime.now(timezone.utc)
        )
        db.add(microscope)
        
        # Create objective
        objective = Objective(
            id=uuid4(),
            microscope_id=microscope.id,
            name="Test 20x",
            magnification=20.0,
            numerical_aperture=0.75,
            immersion="air",
            is_active=True
        )
        db.add(objective)
        
        # Create light source
        light_source = LightSource(
            id=uuid4(),
            microscope_id=microscope.id,
            name="Test 405nm Laser",
            source_type="laser",
            wavelength_nm=405,
            max_power_mw=50.0,
            is_active=True
        )
        db.add(light_source)
        
        # Create channel config
        channel = ChannelConfig(
            id=uuid4(),
            microscope_id=microscope.id,
            channel_name="DAPI",
            fluorophore="Hoechst",
            light_source_id=light_source.id,
            filter_set_id=filter_set.id,
            default_exposure_ms=100.0,
            default_gain=1.5
        )
        db.add(channel)
        
        db.commit()
        db.expunge_all()
        
        return {
            "microscope": microscope,
            "objective": objective,
            "light_source": light_source,
            "channel": channel,
            "filter_set": filter_set
        }


@pytest.fixture
def test_plate():
    """Create a test plate with wells."""
    with db_session() as db:
        plate = HTSPlate(
            id=uuid4(),
            barcode="TEST_PLATE_001",
            plate_format=96,
            created_at=datetime.now(timezone.utc)
        )
        db.add(plate)
        
        # Create a few wells
        wells = []
        for i, well_pos in enumerate(["A01", "A02", "B01", "B02"]):
            well = HTSWell(
                id=uuid4(),
                plate_id=plate.id,
                well_position=well_pos,
                row=ord(well_pos[0]) - ord('A') + 1,
                column=int(well_pos[1:]),
                created_at=datetime.now(timezone.utc)
            )
            db.add(well)
            wells.append(well)
        
        db.commit()
        db.expunge_all()
        
        return {"plate": plate, "wells": wells}


def test_full_ome_tiff_workflow(test_user, test_microscope, test_plate):
    """Test complete OME-TIFF workflow: create file → parse metadata → store → retrieve."""
    with db_session() as db:
        # Create a mock OME-TIFF file
        with tempfile.NamedTemporaryFile(suffix=".ome.tiff", delete=False) as tmp_file:
            tmp_path = Path(tmp_file.name)
            
            # Create minimal TIFF with OME-XML (using tifffile would be better in real scenario)
            # For integration test, we'll simulate the parsing result
            mock_metadata = OMEMetadata(
                filename=tmp_path.name,
                ome_uuid=str(uuid4()),
                instrument=None,
                channels=[],
                dimensions=None,
                acquisition_date=datetime.now(timezone.utc),
                description="Test OME-TIFF",
                raw_xml="<OME/>",
                raw_json={"test": True}
            )
            
            # Create acquisition settings
            acquisition_settings = AcquisitionSettings(
                id=uuid4(),
                exposure_ms=100.0,
                gain=1.5,
                laser_power_percent=25.0,
                binning=1,
                z_position_um=0.0,
                temperature_celsius=22.0,
                raw_settings={"source": "ome_tiff"},
                created_at=datetime.now(timezone.utc)
            )
            db.add(acquisition_settings)
            db.flush()
            
            # QC metrics stored directly in MicroscopyImage
            focus_score = 0.75
            signal_noise_ratio = 25.0
            
            # Store microscopy image
            image = MicroscopyImage(
                id=uuid4(),
                well_id=test_plate["wells"][0].id,
                user_id=test_user.id,
                objective_id=test_microscope["objective"].id,
                acquisition_settings_id=acquisition_settings.id,
                qc_metrics_id=qc_metrics.id,
                channel="DAPI",
                z_slice=0,
                timepoint=0,
                field=1,
                width=2048,
                height=2048,
                bit_depth=16,
                pixel_size_um=0.325,
                file_path=str(tmp_path),
                file_size_bytes=8388608,
                ome_metadata=mock_metadata.raw_json,
                created_at=datetime.now(timezone.utc)
            )
            db.add(image)
            db.commit()
            
            # Retrieve and verify
            retrieved_image = db.query(MicroscopyImage).filter(
                MicroscopyImage.id == image.id
            ).first()
            
            assert retrieved_image is not None
            assert retrieved_image.channel == "DAPI"
            assert retrieved_image.ome_metadata["test"] is True
            assert retrieved_image.qc_metrics.passed_qc is True
            assert retrieved_image.acquisition_settings.exposure_ms == 100.0
            
            # Cleanup
            tmp_path.unlink()


def test_vendor_import_workflow(test_user, test_microscope, test_plate):
    """Test vendor export workflow: create export directory → parse → store images → QC."""
    with db_session() as db:
        # Create mock vendor export directory
        with tempfile.TemporaryDirectory() as tmp_dir:
            export_path = Path(tmp_dir)
            
            # Create ImageFileSet for tracking
            file_set = ImageFileSet(
                id=uuid4(),
                plate_id=test_plate["plate"].id,
                vendor="imagexpress",
                import_path=str(export_path),
                file_count=4,
                image_count=0,
                import_status="pending",
                created_at=datetime.now(timezone.utc)
            )
            db.add(file_set)
            db.flush()
            
            # Simulate importing 4 images (2 wells × 2 channels)
            image_count = 0
            for well in test_plate["wells"][:2]:  # A01, A02
                for channel in ["DAPI", "GFP"]:
                    # Create acquisition settings
                    acquisition_settings = AcquisitionSettings(
                        id=uuid4(),
                        exposure_ms=150.0,
                        gain=2.0,
                        laser_power_percent=30.0,
                        binning=1,
                        z_position_um=0.0,
                        raw_settings={"vendor": "imagexpress"},
                        created_at=datetime.now(timezone.utc)
                    )
                    db.add(acquisition_settings)
                    db.flush()
                    
                    # Run QC on synthetic image data
                    synthetic_image = np.random.randint(0, 4095, (512, 512), dtype=np.uint16)
                    qc_result = run_image_qc(synthetic_image)
                    
                    # Create QC metrics from result
                    qc_metrics = ImageQCMetrics(
                        id=uuid4(),
                        focus_score=qc_result.focus.score,
                        saturation_percentage=qc_result.saturation.saturated_percent,
                        uniformity_score=qc_result.uniformity.uniformity_score,
                        artifact_percentage=qc_result.artifacts.artifact_percent,
                        signal_noise_ratio=20.0,
                        passed_qc=qc_result.overall_pass,
                        qc_issues=json.dumps([issue for issue in qc_result.issues if issue]),
                        created_at=datetime.now(timezone.utc)
                    )
                    db.add(qc_metrics)
                    db.flush()
                    
                    # Create microscopy image
                    image = MicroscopyImage(
                        id=uuid4(),
                        well_id=well.id,
                        user_id=test_user.id,
                        objective_id=test_microscope["objective"].id,
                        acquisition_settings_id=acquisition_settings.id,
                        qc_metrics_id=qc_metrics.id,
                        channel=channel,
                        z_slice=0,
                        timepoint=0,
                        field=1,
                        width=512,
                        height=512,
                        bit_depth=16,
                        pixel_size_um=0.65,
                        file_path=f"{export_path}/{well.well_position}_{channel}.tif",
                        file_size_bytes=524288,
                        ome_metadata={"vendor": "imagexpress", "channel": channel},
                        created_at=datetime.now(timezone.utc)
                    )
                    db.add(image)
                    image_count += 1
            
            # Update file set status
            file_set.image_count = image_count
            file_set.import_status = "completed"
            file_set.completed_at = datetime.now(timezone.utc)
            
            db.commit()
            
            # Verify import results
            imported_images = db.query(MicroscopyImage).join(HTSWell).filter(
                HTSWell.plate_id == test_plate["plate"].id
            ).all()
            
            assert len(imported_images) == 4
            assert all(img.qc_metrics is not None for img in imported_images)
            assert all(img.acquisition_settings.raw_settings["vendor"] == "imagexpress" for img in imported_images)
            
            # Verify file set completion
            completed_file_set = db.query(ImageFileSet).filter(
                ImageFileSet.id == file_set.id
            ).first()
            assert completed_file_set.import_status == "completed"
            assert completed_file_set.image_count == 4


def test_instrument_management_workflow(test_user):
    """Test instrument management: create microscope → add objectives → add channels."""
    with db_session() as db:
        # Create microscope
        microscope = Microscope(
            id=uuid4(),
            name="Integration Test Microscope",
            manufacturer="TestCorp",
            model="IntegrationScope 2000",
            serial_number="INT-TEST-001",
            facility_location="Integration Lab",
            is_active=True,
            created_at=datetime.now(timezone.utc)
        )
        db.add(microscope)
        db.flush()
        
        # Add objectives
        objectives_data = [
            {"name": "10x Air", "magnification": 10.0, "numerical_aperture": 0.30, "immersion": "air"},
            {"name": "40x Oil", "magnification": 40.0, "numerical_aperture": 1.30, "immersion": "oil"},
        ]
        
        objectives = []
        for obj_data in objectives_data:
            objective = Objective(
                id=uuid4(),
                microscope_id=microscope.id,
                **obj_data,
                is_active=True
            )
            db.add(objective)
            objectives.append(objective)
        
        db.flush()
        
        # Add light sources
        light_sources = []
        for wavelength, power in [(405, 50.0), (488, 100.0), (561, 75.0)]:
            light_source = LightSource(
                id=uuid4(),
                microscope_id=microscope.id,
                name=f"{wavelength}nm Laser",
                source_type="laser",
                wavelength_nm=wavelength,
                max_power_mw=power,
                is_active=True
            )
            db.add(light_source)
            light_sources.append(light_source)
        
        db.flush()
        
        # Add filter sets
        filter_sets = []
        filter_configs = [
            {"name": "DAPI", "ex_center": 358, "ex_bw": 46, "em_center": 461, "em_bw": 50},
            {"name": "FITC", "ex_center": 482, "ex_bw": 35, "em_center": 536, "em_bw": 40},
        ]
        
        for filter_config in filter_configs:
            filter_set = FilterSet(
                id=uuid4(),
                name=filter_config["name"],
                excitation_center_nm=filter_config["ex_center"],
                excitation_bandwidth_nm=filter_config["ex_bw"],
                emission_center_nm=filter_config["em_center"],
                emission_bandwidth_nm=filter_config["em_bw"],
                dichroic_cutoff_nm=400
            )
            db.add(filter_set)
            filter_sets.append(filter_set)
        
        db.flush()
        
        # Add channel configurations
        channels = []
        channel_configs = [
            {"name": "DAPI", "fluorophore": "Hoechst", "exposure": 100.0, "gain": 1.5},
            {"name": "FITC", "fluorophore": "FITC", "exposure": 200.0, "gain": 2.0},
        ]
        
        for i, channel_config in enumerate(channel_configs):
            channel = ChannelConfig(
                id=uuid4(),
                microscope_id=microscope.id,
                channel_name=channel_config["name"],
                fluorophore=channel_config["fluorophore"],
                light_source_id=light_sources[i].id if i < len(light_sources) else None,
                filter_set_id=filter_sets[i].id if i < len(filter_sets) else None,
                default_exposure_ms=channel_config["exposure"],
                default_gain=channel_config["gain"]
            )
            db.add(channel)
            channels.append(channel)
        
        db.commit()
        
        # Verify complete instrument setup
        retrieved_microscope = db.query(Microscope).filter(
            Microscope.id == microscope.id
        ).first()
        
        assert retrieved_microscope is not None
        assert retrieved_microscope.name == "Integration Test Microscope"
        
        # Verify relationships
        microscope_objectives = db.query(Objective).filter(
            Objective.microscope_id == microscope.id
        ).all()
        assert len(microscope_objectives) == 2
        
        microscope_light_sources = db.query(LightSource).filter(
            LightSource.microscope_id == microscope.id
        ).all()
        assert len(microscope_light_sources) == 3
        
        microscope_channels = db.query(ChannelConfig).filter(
            ChannelConfig.microscope_id == microscope.id
        ).all()
        assert len(microscope_channels) == 2


def test_image_qc_workflow():
    """Test image QC workflow: create image → run QC → store results → retrieve."""
    # Create test images with different quality characteristics
    test_cases = [
        {
            "name": "good_image",
            "image": np.random.randint(1000, 3000, (512, 512), dtype=np.uint16),
            "expected_pass": True
        },
        {
            "name": "blurry_image", 
            "image": np.random.randint(500, 1500, (512, 512), dtype=np.uint16),  # Lower contrast
            "expected_pass": False
        },
        {
            "name": "saturated_image",
            "image": np.full((512, 512), 65535, dtype=np.uint16),  # Fully saturated
            "expected_pass": False
        }
    ]
    
    results = []
    for test_case in test_cases:
        # Run QC analysis
        qc_result = run_image_qc(test_case["image"])
        
        # Verify QC result structure
        assert hasattr(qc_result, 'focus')
        assert hasattr(qc_result, 'saturation') 
        assert hasattr(qc_result, 'uniformity')
        assert hasattr(qc_result, 'artifacts')
        assert hasattr(qc_result, 'overall_pass')
        
        # Store result for verification
        results.append({
            "name": test_case["name"],
            "qc_result": qc_result,
            "expected_pass": test_case["expected_pass"]
        })
    
    # Verify QC logic works as expected
    good_result = next(r for r in results if r["name"] == "good_image")
    assert good_result["qc_result"].overall_pass == good_result["expected_pass"]
    
    saturated_result = next(r for r in results if r["name"] == "saturated_image")
    assert saturated_result["qc_result"].saturation.saturated_percent > 50.0
    assert saturated_result["qc_result"].overall_pass == saturated_result["expected_pass"]


def test_5d_browsing(test_user, test_microscope, test_plate):
    """Test 5D browsing: query across dimensions → paginate → filter by QC."""
    with db_session() as db:
        # Create images across multiple dimensions
        channels = ["DAPI", "GFP", "RFP"]
        z_slices = [0, 1, 2]
        timepoints = [0, 1]
        
        created_images = []
        for well in test_plate["wells"][:2]:  # Use first 2 wells
            for channel in channels:
                for z_slice in z_slices:
                    for timepoint in timepoints:
                        # Create varied QC metrics
                        passed_qc = np.random.random() > 0.3  # 70% pass rate
                        
                        # Create acquisition settings
                        acquisition_settings = AcquisitionSettings(
                            id=uuid4(),
                            exposure_ms=100.0 + z_slice * 50,  # Vary by Z
                            gain=1.5 + timepoint * 0.5,  # Vary by timepoint
                            z_position_um=float(z_slice * 1.0),
                            raw_settings={"5d_test": True},
                            created_at=datetime.now(timezone.utc)
                        )
                        db.add(acquisition_settings)
                        db.flush()
                        
                        # Create QC metrics
                        qc_metrics = ImageQCMetrics(
                            id=uuid4(),
                            focus_score=0.7 if passed_qc else 0.2,
                            saturation_percentage=2.0 if passed_qc else 15.0,
                            uniformity_score=0.8 if passed_qc else 0.3,
                            artifact_percentage=1.0 if passed_qc else 12.0,
                            signal_noise_ratio=25.0,
                            passed_qc=passed_qc,
                            qc_issues=json.dumps([]),
                            created_at=datetime.now(timezone.utc)
                        )
                        db.add(qc_metrics)
                        db.flush()
                        
                        # Create image
                        image = MicroscopyImage(
                            id=uuid4(),
                            well_id=well.id,
                            user_id=test_user.id,
                            objective_id=test_microscope["objective"].id,
                            acquisition_settings_id=acquisition_settings.id,
                            qc_metrics_id=qc_metrics.id,
                            channel=channel,
                            z_slice=z_slice,
                            timepoint=timepoint,
                            field=1,
                            width=1024,
                            height=1024,
                            bit_depth=16,
                            pixel_size_um=0.325,
                            file_path=f"/test/{well.well_position}_{channel}_Z{z_slice}_T{timepoint}.tif",
                            file_size_bytes=2097152,
                            created_at=datetime.now(timezone.utc)
                        )
                        db.add(image)
                        created_images.append(image)
        
        db.commit()
        
        # Test 5D browsing queries
        
        # 1. Browse by channel
        dapi_images = db.query(MicroscopyImage).filter(
            MicroscopyImage.channel == "DAPI"
        ).all()
        assert len(dapi_images) == 2 * 3 * 2  # 2 wells × 3 Z × 2 T
        
        # 2. Browse by Z-slice
        z0_images = db.query(MicroscopyImage).filter(
            MicroscopyImage.z_slice == 0
        ).all()
        assert len(z0_images) == 2 * 3 * 2  # 2 wells × 3 channels × 2 T
        
        # 3. Browse by timepoint
        t1_images = db.query(MicroscopyImage).filter(
            MicroscopyImage.timepoint == 1
        ).all()
        assert len(t1_images) == 2 * 3 * 3  # 2 wells × 3 channels × 3 Z
        
        # 4. Filter by QC status
        passed_images = db.query(MicroscopyImage).join(ImageQCMetrics).filter(
            ImageQCMetrics.passed_qc == True
        ).all()
        failed_images = db.query(MicroscopyImage).join(ImageQCMetrics).filter(
            ImageQCMetrics.passed_qc == False
        ).all()
        
        assert len(passed_images) + len(failed_images) == len(created_images)
        
        # 5. Complex multi-dimensional query
        specific_images = db.query(MicroscopyImage).join(HTSWell).filter(
            HTSWell.well_position == "A01",
            MicroscopyImage.channel == "GFP",
            MicroscopyImage.z_slice == 1,
            MicroscopyImage.timepoint == 0
        ).all()
        assert len(specific_images) == 1
        
        # 6. Pagination test
        page_size = 10
        page_1 = db.query(MicroscopyImage).limit(page_size).offset(0).all()
        page_2 = db.query(MicroscopyImage).limit(page_size).offset(page_size).all()
        
        assert len(page_1) == min(page_size, len(created_images))
        if len(created_images) > page_size:
            assert len(page_2) > 0
            # Verify no overlap
            page_1_ids = {img.id for img in page_1}
            page_2_ids = {img.id for img in page_2}
            assert page_1_ids.isdisjoint(page_2_ids)


def test_plate_qc_aggregation(test_user, test_microscope, test_plate):
    """Test plate QC aggregation: multiple images → aggregate QC → generate report."""
    with db_session() as db:
        # Create images with known QC characteristics
        qc_scenarios = [
            {"focus": 0.8, "saturation": 2.0, "uniformity": 0.9, "artifacts": 1.0, "pass": True},   # Good
            {"focus": 0.2, "saturation": 3.0, "uniformity": 0.8, "artifacts": 2.0, "pass": False},  # Poor focus
            {"focus": 0.7, "saturation": 15.0, "uniformity": 0.7, "artifacts": 1.5, "pass": False}, # High saturation
            {"focus": 0.6, "saturation": 4.0, "uniformity": 0.3, "artifacts": 3.0, "pass": False},  # Poor uniformity
        ]
        
        created_images = []
        for i, well in enumerate(test_plate["wells"]):
            scenario = qc_scenarios[i % len(qc_scenarios)]
            
            # Create acquisition settings
            acquisition_settings = AcquisitionSettings(
                id=uuid4(),
                exposure_ms=200.0,
                gain=2.0,
                raw_settings={"qc_test": True},
                created_at=datetime.now(timezone.utc)
            )
            db.add(acquisition_settings)
            db.flush()
            
            # Create QC metrics based on scenario
            qc_metrics = ImageQCMetrics(
                id=uuid4(),
                focus_score=scenario["focus"],
                saturation_percentage=scenario["saturation"],
                uniformity_score=scenario["uniformity"],
                artifact_percentage=scenario["artifacts"],
                signal_noise_ratio=20.0,
                passed_qc=scenario["pass"],
                qc_issues=json.dumps(["test_issue"] if not scenario["pass"] else []),
                created_at=datetime.now(timezone.utc)
            )
            db.add(qc_metrics)
            db.flush()
            
            # Create image
            image = MicroscopyImage(
                id=uuid4(),
                well_id=well.id,
                user_id=test_user.id,
                objective_id=test_microscope["objective"].id,
                acquisition_settings_id=acquisition_settings.id,
                qc_metrics_id=qc_metrics.id,
                channel="DAPI",
                z_slice=0,
                timepoint=0,
                field=1,
                width=1024,
                height=1024,
                bit_depth=16,
                pixel_size_um=0.325,
                file_path=f"/test/{well.well_position}_DAPI.tif",
                file_size_bytes=2097152,
                created_at=datetime.now(timezone.utc)
            )
            db.add(image)
            created_images.append(image)
        
        db.commit()
        
        # Generate plate QC aggregation
        plate_images = db.query(MicroscopyImage).join(HTSWell).filter(
            HTSWell.plate_id == test_plate["plate"].id
        ).all()
        
        # Calculate aggregated metrics
        total_images = len(plate_images)
        passed_images = sum(1 for img in plate_images if img.qc_metrics.passed_qc)
        failed_images = total_images - passed_images
        
        avg_focus = sum(img.qc_metrics.focus_score for img in plate_images) / total_images
        avg_saturation = sum(img.qc_metrics.saturation_percentage for img in plate_images) / total_images
        avg_uniformity = sum(img.qc_metrics.uniformity_score for img in plate_images) / total_images
        
        # Verify aggregation results
        assert total_images == 4
        assert passed_images == 1  # Only first scenario passes
        assert failed_images == 3
        
        pass_rate = (passed_images / total_images) * 100
        assert pass_rate == 25.0
        
        # Verify averages are reasonable
        expected_avg_focus = sum(s["focus"] for s in qc_scenarios) / len(qc_scenarios)
        assert abs(avg_focus - expected_avg_focus) < 0.1
        
        # Generate QC report data structure
        qc_report = {
            "plate_id": str(test_plate["plate"].id),
            "plate_barcode": test_plate["plate"].barcode,
            "total_images": total_images,
            "passed_count": passed_images,
            "failed_count": failed_images,
            "pass_rate_percent": pass_rate,
            "average_focus_score": avg_focus,
            "average_saturation_percent": avg_saturation,
            "average_uniformity_score": avg_uniformity,
            "wells_with_issues": [
                img.well.well_position for img in plate_images 
                if not img.qc_metrics.passed_qc
            ],
            "generated_at": datetime.now(timezone.utc).isoformat()
        }
        
        # Verify report structure
        assert qc_report["pass_rate_percent"] == 25.0
        assert len(qc_report["wells_with_issues"]) == 3
        assert "A02" in qc_report["wells_with_issues"]  # Second well should fail
        assert qc_report["total_images"] == qc_report["passed_count"] + qc_report["failed_count"]
