#!/usr/bin/env python3
"""
Seed imaging metadata demo data (idempotent).

Creates:
- Microscope instruments with objectives and channel configurations
- Demo plate with wells and microscopy images
- Image metadata with varied QC scores and acquisition settings
- Sample OME-TIFF metadata for testing

Usage:
  python scripts/seed_imaging_metadata.py
  python scripts/seed_imaging_metadata.py --reset
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path
from datetime import datetime, timezone
from typing import Dict, List, Optional, Tuple
from uuid import uuid4, UUID

import numpy as np

# Ensure repo root is on sys.path when running as `python scripts/seed_imaging_metadata.py`
REPO_ROOT = Path(__file__).resolve().parents[1]
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

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


# Demo microscope configurations
DEMO_MICROSCOPES = [
    {
        "name": "Nikon Ti2-E #1",
        "manufacturer": "Nikon",
        "model": "Ti2-E",
        "serial_number": "NK-TI2E-001",
        "facility_location": "Lab A, Room 101",
        "objectives": [
            {"name": "Plan Apo 10x", "magnification": 10.0, "numerical_aperture": 0.45, "immersion": "air"},
            {"name": "Plan Apo 20x", "magnification": 20.0, "numerical_aperture": 0.75, "immersion": "air"},
            {"name": "Plan Apo 40x Oil", "magnification": 40.0, "numerical_aperture": 1.30, "immersion": "oil"},
            {"name": "Plan Apo 63x Oil", "magnification": 63.0, "numerical_aperture": 1.40, "immersion": "oil"},
        ],
        "light_sources": [
            {"name": "405nm Diode Laser", "source_type": "laser", "wavelength_nm": 405, "max_power_mw": 50.0},
            {"name": "488nm Argon Laser", "source_type": "laser", "wavelength_nm": 488, "max_power_mw": 100.0},
            {"name": "561nm DPSS Laser", "source_type": "laser", "wavelength_nm": 561, "max_power_mw": 75.0},
            {"name": "640nm Diode Laser", "source_type": "laser", "wavelength_nm": 640, "max_power_mw": 80.0},
        ],
        "channels": [
            {"channel_name": "DAPI", "fluorophore": "Hoechst 33342", "default_exposure_ms": 100.0, "default_gain": 1.5},
            {"channel_name": "GFP", "fluorophore": "FITC", "default_exposure_ms": 200.0, "default_gain": 2.0},
            {"channel_name": "RFP", "fluorophore": "Texas Red", "default_exposure_ms": 300.0, "default_gain": 2.5},
            {"channel_name": "Cy5", "fluorophore": "Alexa 647", "default_exposure_ms": 400.0, "default_gain": 3.0},
        ]
    },
    {
        "name": "Zeiss Axio Observer Z1",
        "manufacturer": "Zeiss",
        "model": "Axio Observer Z1",
        "serial_number": "ZS-AO-Z1-002",
        "facility_location": "Lab B, Room 203",
        "objectives": [
            {"name": "EC Plan-Neofluar 10x", "magnification": 10.0, "numerical_aperture": 0.30, "immersion": "air"},
            {"name": "Plan-Apochromat 20x", "magnification": 20.0, "numerical_aperture": 0.80, "immersion": "air"},
            {"name": "Plan-Apochromat 63x Oil", "magnification": 63.0, "numerical_aperture": 1.40, "immersion": "oil"},
        ],
        "light_sources": [
            {"name": "Colibri 7 LED", "source_type": "led", "wavelength_nm": None, "max_power_mw": None},
        ],
        "channels": [
            {"channel_name": "DAPI", "fluorophore": "DAPI", "default_exposure_ms": 150.0, "default_gain": 1.8},
            {"channel_name": "GFP", "fluorophore": "GFP", "default_exposure_ms": 250.0, "default_gain": 2.2},
            {"channel_name": "RFP", "fluorophore": "mCherry", "default_exposure_ms": 350.0, "default_gain": 2.8},
        ]
    },
    {
        "name": "Opera Phenix Plus",
        "manufacturer": "PerkinElmer",
        "model": "Opera Phenix Plus",
        "serial_number": "PE-OPP-003",
        "facility_location": "HCS Core Facility",
        "objectives": [
            {"name": "Air 20x", "magnification": 20.0, "numerical_aperture": 0.70, "immersion": "air"},
            {"name": "Water 40x", "magnification": 40.0, "numerical_aperture": 1.15, "immersion": "water"},
            {"name": "Water 63x", "magnification": 63.0, "numerical_aperture": 1.15, "immersion": "water"},
        ],
        "light_sources": [
            {"name": "405nm Laser", "source_type": "laser", "wavelength_nm": 405, "max_power_mw": 30.0},
            {"name": "488nm Laser", "source_type": "laser", "wavelength_nm": 488, "max_power_mw": 50.0},
            {"name": "561nm Laser", "source_type": "laser", "wavelength_nm": 561, "max_power_mw": 40.0},
            {"name": "640nm Laser", "source_type": "laser", "wavelength_nm": 640, "max_power_mw": 45.0},
        ],
        "channels": [
            {"channel_name": "DAPI", "fluorophore": "Hoechst", "default_exposure_ms": 80.0, "default_gain": 1.2},
            {"channel_name": "FITC", "fluorophore": "FITC", "default_exposure_ms": 120.0, "default_gain": 1.8},
            {"channel_name": "TRITC", "fluorophore": "TRITC", "default_exposure_ms": 180.0, "default_gain": 2.1},
            {"channel_name": "Cy5", "fluorophore": "Cy5", "default_exposure_ms": 220.0, "default_gain": 2.5},
            {"channel_name": "Brightfield", "fluorophore": None, "default_exposure_ms": 50.0, "default_gain": 1.0},
        ]
    }
]

# Filter sets for different fluorophores
DEMO_FILTER_SETS = [
    {"name": "DAPI", "excitation_center_nm": 358, "excitation_bandwidth_nm": 46, 
     "emission_center_nm": 461, "emission_bandwidth_nm": 50, "dichroic_cutoff_nm": 400},
    {"name": "FITC", "excitation_center_nm": 482, "excitation_bandwidth_nm": 35,
     "emission_center_nm": 536, "emission_bandwidth_nm": 40, "dichroic_cutoff_nm": 506},
    {"name": "TRITC", "excitation_center_nm": 543, "excitation_bandwidth_nm": 22,
     "emission_center_nm": 593, "emission_bandwidth_nm": 40, "dichroic_cutoff_nm": 562},
    {"name": "Texas Red", "excitation_center_nm": 562, "excitation_bandwidth_nm": 40,
     "emission_center_nm": 624, "emission_bandwidth_nm": 40, "dichroic_cutoff_nm": 593},
    {"name": "Cy5", "excitation_center_nm": 628, "excitation_bandwidth_nm": 40,
     "emission_center_nm": 692, "emission_bandwidth_nm": 40, "dichroic_cutoff_nm": 660},
]


def seed_filter_sets(db: db_session) -> Dict[str, UUID]:
    """Create demo filter sets."""
    print("Creating filter sets...")
    
    filter_set_ids = {}
    
    for filter_data in DEMO_FILTER_SETS:
        # Check if filter set already exists
        existing = db.query(FilterSet).filter(FilterSet.name == filter_data["name"]).first()
        if existing:
            filter_set_ids[filter_data["name"]] = existing.id
            continue
        
        filter_set = FilterSet(
            id=uuid4(),
            **filter_data
        )
        db.add(filter_set)
        filter_set_ids[filter_data["name"]] = filter_set.id
    
    db.commit()
    print(f"Created {len(filter_set_ids)} filter sets")
    return filter_set_ids


def seed_microscopes(db: db_session, filter_set_ids: Dict[str, UUID]) -> List[UUID]:
    """Create demo microscopes with objectives and channels."""
    print("Creating microscopes...")
    
    microscope_ids = []
    
    for microscope_data in DEMO_MICROSCOPES:
        # Check if microscope already exists
        existing = db.query(Microscope).filter(
            Microscope.name == microscope_data["name"]
        ).first()
        if existing:
            microscope_ids.append(existing.id)
            continue
        
        # Create microscope
        microscope = Microscope(
            id=uuid4(),
            name=microscope_data["name"],
            manufacturer=microscope_data["manufacturer"],
            model=microscope_data["model"],
            serial_number=microscope_data.get("serial_number"),
            facility_location=microscope_data.get("facility_location"),
            is_active=True,
            created_at=datetime.now(timezone.utc)
        )
        db.add(microscope)
        microscope_ids.append(microscope.id)
        
        # Create objectives
        for obj_data in microscope_data["objectives"]:
            objective = Objective(
                id=uuid4(),
                microscope_id=microscope.id,
                name=obj_data["name"],
                magnification=obj_data["magnification"],
                numerical_aperture=obj_data["numerical_aperture"],
                immersion=obj_data["immersion"],
                working_distance_mm=None,
                correction="Plan",
                is_active=True
            )
            db.add(objective)
        
        # Create light sources
        for light_data in microscope_data["light_sources"]:
            light_source = LightSource(
                id=uuid4(),
                microscope_id=microscope.id,
                name=light_data["name"],
                source_type=light_data["source_type"],
                wavelength_nm=light_data.get("wavelength_nm"),
                max_power_mw=light_data.get("max_power_mw"),
                is_active=True
            )
            db.add(light_source)
        
        # Create channel configurations
        for channel_data in microscope_data["channels"]:
            # Find matching filter set
            filter_set_id = None
            channel_name = channel_data["channel_name"]
            if channel_name in filter_set_ids:
                filter_set_id = filter_set_ids[channel_name]
            elif channel_name == "FITC" and "FITC" in filter_set_ids:
                filter_set_id = filter_set_ids["FITC"]
            elif channel_name == "TRITC" and "TRITC" in filter_set_ids:
                filter_set_id = filter_set_ids["TRITC"]
            
            channel = ChannelConfig(
                id=uuid4(),
                microscope_id=microscope.id,
                channel_name=channel_data["channel_name"],
                fluorophore=channel_data.get("fluorophore"),
                light_source_id=None,  # Would link to specific light source
                filter_set_id=filter_set_id,
                default_exposure_ms=channel_data.get("default_exposure_ms"),
                default_gain=channel_data.get("default_gain")
            )
            db.add(channel)
    
    db.commit()
    print(f"Created {len(microscope_ids)} microscopes")
    return microscope_ids


def seed_demo_plate_and_wells(db: db_session) -> Tuple[UUID, List[UUID]]:
    """Create demo plate with wells."""
    print("Creating demo plate and wells...")
    
    # Check if demo plate already exists
    existing_plate = db.query(HTSPlate).filter(HTSPlate.barcode == "IMAGING_DEMO_001").first()
    if existing_plate:
        well_ids = [well.id for well in existing_plate.wells]
        print(f"Found existing plate with {len(well_ids)} wells")
        return existing_plate.id, well_ids
    
    # Create demo plate
    plate = HTSPlate(
        id=uuid4(),
        barcode="IMAGING_DEMO_001",
        plate_format=96,
        created_at=datetime.now(timezone.utc)
    )
    db.add(plate)
    
    # Create wells (A01-H12 for 96-well plate)
    well_ids = []
    for row in range(8):  # A-H
        for col in range(12):  # 1-12
            well_position = f"{chr(ord('A') + row)}{col+1:02d}"
            
            well = HTSWell(
                id=uuid4(),
                plate_id=plate.id,
                well_position=well_position,
                row=row + 1,
                column=col + 1,
                created_at=datetime.now(timezone.utc)
            )
            db.add(well)
            well_ids.append(well.id)
    
    db.commit()
    print(f"Created plate with {len(well_ids)} wells")
    return plate.id, well_ids


def generate_varied_qc_metrics() -> Tuple[float, float]:
    """Generate varied QC metrics for demo images."""
    # Focus score (0-1, higher is better)
    focus_score = np.random.beta(2, 1)  # Skewed toward higher values
    
    # Signal to noise ratio (5-50, higher is better)
    signal_noise_ratio = np.random.uniform(5.0, 50.0)
    
    return focus_score, signal_noise_ratio


def create_sample_ome_metadata(channel: str, z_slice: int, timepoint: int) -> Dict:
    """Create sample OME-TIFF metadata."""
    return {
        "ome_uuid": str(uuid4()),
        "instrument": {
            "microscope_name": "Nikon Ti2-E #1",
            "microscope_model": "Ti2-E",
            "objective_name": "Plan Apo 20x",
            "objective_magnification": 20.0,
            "objective_na": 0.75,
            "objective_immersion": "air"
        },
        "channels": [{
            "name": channel,
            "fluorophore": {"DAPI": "Hoechst 33342", "GFP": "FITC", "RFP": "Texas Red"}.get(channel),
            "excitation_wavelength_nm": {"DAPI": 358, "GFP": 482, "RFP": 562}.get(channel),
            "emission_wavelength_nm": {"DAPI": 461, "GFP": 536, "RFP": 593}.get(channel),
            "exposure_ms": {"DAPI": 100.0, "GFP": 200.0, "RFP": 300.0}.get(channel, 200.0),
            "color": {"DAPI": "#0000FF", "GFP": "#00FF00", "RFP": "#FF0000"}.get(channel, "#FFFFFF")
        }],
        "dimensions": {
            "size_x": 2048,
            "size_y": 2048,
            "size_z": 5,
            "size_c": 1,
            "size_t": 3,
            "pixel_size_x_um": 0.325,
            "pixel_size_y_um": 0.325,
            "pixel_size_z_um": 1.0,
            "pixel_type": "uint16"
        },
        "acquisition_date": datetime.now(timezone.utc).isoformat(),
        "description": f"Demo image {channel} Z{z_slice} T{timepoint}"
    }


def seed_demo_images(db: db_session, microscope_ids: List[UUID], plate_id: UUID, well_ids: List[UUID]) -> None:
    """Create demo images with varied QC scores and metadata."""
    print("Creating demo images...")
    
    # Get first user for ownership
    user = db.query(User).first()
    if not user:
        print("Warning: No users found, creating demo user")
        user = User(
            id=uuid4(),
            email="demo@example.com",
            is_active=True,
            created_at=datetime.now(timezone.utc)
        )
        db.add(user)
        db.commit()
    
    # Get objectives for acquisition settings
    objectives = db.query(Objective).filter(Objective.microscope_id.in_(microscope_ids)).all()
    if not objectives:
        print("Warning: No objectives found")
        return
    
    channels = ["DAPI", "GFP", "RFP"]
    z_slices = [0, 1, 2, 3, 4]
    timepoints = [0, 1, 2]
    
    image_count = 0
    # Create images for first 24 wells (3 rows) to keep demo manageable
    for well_id in well_ids[:24]:
        well = db.query(HTSWell).filter(HTSWell.id == well_id).first()
        if not well:
            continue
            
        for channel in channels:
            for z_slice in z_slices:
                for timepoint in timepoints:
                    # Check if image already exists
                    existing = db.query(MicroscopyImage).filter(
                        MicroscopyImage.well_id == well_id,
                        MicroscopyImage.channel == channel,
                        MicroscopyImage.z_slice == z_slice,
                        MicroscopyImage.timepoint == timepoint
                    ).first()
                    if existing:
                        continue
                    
                    # Create acquisition settings
                    objective = np.random.choice(objectives)
                    acquisition_settings = AcquisitionSettings(
                        id=uuid4(),
                        exposure_ms=float(np.random.uniform(50, 500)),
                        gain=float(np.random.uniform(1.0, 4.0)),
                        laser_power_percent=float(np.random.uniform(5, 50)),
                        binning=1,
                        z_position_um=float(z_slice * 1.0),
                        temperature_celsius=float(np.random.uniform(20, 25)),
                        raw_settings={"demo": True, "field": 1},
                        created_at=datetime.now(timezone.utc)
                    )
                    db.add(acquisition_settings)
                    db.flush()  # Get ID
                    
                    # Generate QC metrics
                    focus_score, signal_noise_ratio = generate_varied_qc_metrics()
                    
                    # Create microscopy image
                    image = MicroscopyImage(
                        id=uuid4(),
                        well_id=well_id,
                        objective_id=objective.id,
                        acquisition_settings_id=acquisition_settings.id,
                        channel=channel,
                        z_slice=z_slice,
                        timepoint=timepoint,
                        width=2048,
                        height=2048,
                        bit_depth=16,
                        pixel_size_um=0.325,
                        image_path=f"/demo/images/{well.well_position}_{channel}_Z{z_slice}_T{timepoint}.tiff",
                        focus_score=focus_score,
                        signal_noise_ratio=signal_noise_ratio,
                        ome_metadata=create_sample_ome_metadata(channel, z_slice, timepoint),
                        created_at=datetime.now(timezone.utc)
                    )
                    db.add(image)
                    image_count += 1
    
    db.commit()
    print(f"Created {image_count} demo images")


def seed_image_file_set(db: db_session, plate_id: UUID) -> UUID:
    """Create demo image file set for batch tracking."""
    print("Creating image file set...")
    
    # Check if file set already exists
    existing = db.query(ImageFileSet).filter(
        ImageFileSet.plate_id == plate_id,
        ImageFileSet.vendor == "demo"
    ).first()
    if existing:
        print("Found existing file set")
        return existing.id
    
    # Count images for this plate
    image_count = db.query(MicroscopyImage).join(HTSWell).filter(
        HTSWell.plate_id == plate_id
    ).count()
    
    file_set = ImageFileSet(
        id=uuid4(),
        plate_id=plate_id,
        vendor="demo",
        import_path="/demo/exports/IMAGING_DEMO_001",
        file_count=image_count,
        image_count=image_count,
        import_status="completed",
        error_message=None,
        created_at=datetime.now(timezone.utc),
        completed_at=datetime.now(timezone.utc)
    )
    db.add(file_set)
    db.commit()
    
    print(f"Created file set with {image_count} images")
    return file_set.id


def reset_imaging_data(db: db_session) -> None:
    """Reset all imaging demo data."""
    print("Resetting imaging demo data...")
    
    # Delete in reverse dependency order
    db.query(MicroscopyImage).filter(
        MicroscopyImage.image_path.like("/demo/images/%")
    ).delete(synchronize_session=False)
    
    db.query(AcquisitionSettings).filter(
        AcquisitionSettings.raw_settings.like("%demo%")
    ).delete(synchronize_session=False)
    
    db.query(ImageFileSet).filter(
        ImageFileSet.vendor == "demo"
    ).delete(synchronize_session=False)
    
    db.query(HTSWell).join(HTSPlate).filter(
        HTSPlate.barcode == "IMAGING_DEMO_001"
    ).delete(synchronize_session=False)
    
    db.query(HTSPlate).filter(
        HTSPlate.barcode == "IMAGING_DEMO_001"
    ).delete(synchronize_session=False)
    
    # Delete microscope-related data
    microscope_names = [m["name"] for m in DEMO_MICROSCOPES]
    microscope_ids = [
        m.id for m in db.query(Microscope).filter(
            Microscope.name.in_(microscope_names)
        ).all()
    ]
    
    if microscope_ids:
        db.query(ChannelConfig).filter(
            ChannelConfig.microscope_id.in_(microscope_ids)
        ).delete(synchronize_session=False)
        
        db.query(LightSource).filter(
            LightSource.microscope_id.in_(microscope_ids)
        ).delete(synchronize_session=False)
        
        db.query(Objective).filter(
            Objective.microscope_id.in_(microscope_ids)
        ).delete(synchronize_session=False)
        
        db.query(Microscope).filter(
            Microscope.id.in_(microscope_ids)
        ).delete(synchronize_session=False)
    
    # Delete filter sets
    filter_names = [f["name"] for f in DEMO_FILTER_SETS]
    db.query(FilterSet).filter(
        FilterSet.name.in_(filter_names)
    ).delete(synchronize_session=False)
    
    db.commit()
    print("Reset complete")


def main() -> None:
    """Main seeding function."""
    parser = argparse.ArgumentParser(description="Seed imaging metadata demo data")
    parser.add_argument(
        "--reset", 
        action="store_true", 
        help="Reset all demo data before seeding"
    )
    args = parser.parse_args()
    
    print("🔬 Imaging Metadata Demo Data Seeder")
    print("=" * 50)
    
    with db_session() as db:
        if args.reset:
            reset_imaging_data(db)
        
        # Create filter sets
        filter_set_ids = seed_filter_sets(db)
        
        # Create microscopes with objectives and channels
        microscope_ids = seed_microscopes(db, filter_set_ids)
        
        # Create demo plate and wells
        plate_id, well_ids = seed_demo_plate_and_wells(db)
        
        # Create demo images with QC metrics
        seed_demo_images(db, microscope_ids, plate_id, well_ids)
        
        # Create image file set for tracking
        file_set_id = seed_image_file_set(db, plate_id)
        
        print("\n✅ Imaging metadata demo data seeded successfully!")
        print(f"   - {len(DEMO_MICROSCOPES)} microscopes")
        print(f"   - {len(DEMO_FILTER_SETS)} filter sets")
        print(f"   - 1 demo plate with {len(well_ids)} wells")
        print(f"   - ~1080 demo images (24 wells × 3 channels × 5 Z × 3 T)")
        print(f"   - File set: {file_set_id}")
        print("\nReady for dashboard testing and API validation!")


if __name__ == "__main__":
    main()
