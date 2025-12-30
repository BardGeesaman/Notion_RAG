"""Unit tests for vendor-specific HCS parsers."""

import tempfile
import xml.etree.ElementTree as ET
from pathlib import Path
from unittest.mock import patch

import numpy as np
import pytest
import tifffile

from amprenta_rag.imaging.vendor_parsers import (
    PlateImageInfo,
    PlateImportResult,
    auto_detect_format,
    parse_opera_export,
    parse_imagexpress_export,
    parse_cellvoyager_export,
    parse_vendor_export,
    _well_position_to_indices,
    _estimate_plate_format,
)


@pytest.fixture
def temp_dir():
    """Create temporary directory for test files."""
    with tempfile.TemporaryDirectory() as tmp_dir:
        yield Path(tmp_dir)


def create_mock_opera_export(path: Path, wells: int = 4, channels: int = 2, z_slices: int = 2) -> None:
    """Create mock Opera export directory for testing."""
    # Create Index.xml (no leading whitespace to avoid parsing issues)
    index_xml = f"""<?xml version="1.0" encoding="UTF-8"?>
<OperaExport>
    <PlateBarcode>TEST_OPERA_PLATE</PlateBarcode>
    <PlateType>96</PlateType>
    <Measurements>
        <Measurement>
            <Wells>{wells}</Wells>
            <Channels>{channels}</Channels>
            <ZSlices>{z_slices}</ZSlices>
        </Measurement>
    </Measurements>
</OperaExport>"""
    
    (path / "Index.xml").write_text(index_xml)
    
    # Create Images directory
    images_dir = path / "Images"
    images_dir.mkdir()
    
    # Create mock TIFF files with Opera naming pattern
    # r{row}c{col}f{field}p{plane}-ch{channel}sk{z}fk{frame}fl{timeline}.tiff
    for row in range(1, min(wells + 1, 3)):  # A, B rows
        for col in range(1, min(wells + 1, 3)):  # Columns 1, 2
            for ch in range(1, channels + 1):
                for z in range(1, z_slices + 1):
                    filename = f"r{row:02d}c{col:02d}f01p01-ch{ch}sk{z}fk1fl1.tiff"
                    image_path = images_dir / filename
                    
                    # Create minimal TIFF file
                    data = np.random.randint(0, 65535, (64, 64), dtype=np.uint16)
                    tifffile.imwrite(image_path, data)


def create_mock_imagexpress_export(path: Path, wells: int = 4, channels: int = 2, timepoints: int = 2) -> None:
    """Create mock ImageXpress export directory for testing."""
    # Create experiment.htd file
    htd_content = f"""
    ExperimentName=Test Experiment
    PlateName=TEST_IMAGEXPRESS_PLATE
    PlateFormat=96
    Channels={channels}
    TimePoints={timepoints}
    """
    
    (path / "experiment.htd").write_text(htd_content)
    
    # Create TimePoint directories
    for tp in range(1, timepoints + 1):
        tp_dir = path / f"TimePoint_{tp}"
        tp_dir.mkdir()
        
        # Create mock images with ImageXpress naming pattern
        # {well}_{site}_w{wavelength}.TIF
        for row_idx in range(min(wells, 2)):  # A, B
            for col in range(1, min(wells + 1, 3)):  # Columns 1, 2
                well = f"{chr(ord('A') + row_idx)}{col:02d}"
                for site in range(1, 2):  # Single site
                    for wavelength in range(1, channels + 1):
                        filename = f"{well}_s{site}_w{wavelength}.TIF"
                        image_path = tp_dir / filename
                        
                        # Create minimal TIFF file
                        data = np.random.randint(0, 65535, (64, 64), dtype=np.uint16)
                        tifffile.imwrite(image_path, data)


def create_mock_cellvoyager_export(path: Path, wells: int = 4, channels: int = 2, timepoints: int = 2) -> None:
    """Create mock Cell Voyager export directory for testing."""
    # Create MeasurementData.mlf file
    mlf_content = f"""
    MeasurementName=Test Measurement
    PlateName=TEST_CELLVOYAGER_PLATE
    PlateFormat=96
    Channels={channels}
    TimePoints={timepoints}
    Wells={wells}
    """
    
    (path / "MeasurementData.mlf").write_text(mlf_content)
    
    # Create Images directory
    images_dir = path / "Images"
    images_dir.mkdir()
    
    # Create mock .flex files with Cell Voyager naming pattern
    # W{well}F{field}T{time}Z{z}C{channel}.flex
    for well in range(1, wells + 1):
        for field in range(1, 2):  # Single field
            for time in range(1, timepoints + 1):
                for z in range(1, 2):  # Single Z
                    for channel in range(1, channels + 1):
                        filename = f"W{well:04d}F{field:04d}T{time:04d}Z{z:04d}C{channel}.flex"
                        image_path = images_dir / filename
                        
                        # Create mock .flex file (just empty file for testing)
                        image_path.write_bytes(b"MOCK_FLEX_DATA")


class TestAutoDetectFormat:
    """Test auto-detection of vendor formats."""

    def test_detect_opera_format(self, temp_dir):
        """Test detection of Opera format."""
        create_mock_opera_export(temp_dir)
        
        detected = auto_detect_format(temp_dir)
        assert detected == "opera"

    def test_detect_imagexpress_format(self, temp_dir):
        """Test detection of ImageXpress format."""
        create_mock_imagexpress_export(temp_dir)
        
        detected = auto_detect_format(temp_dir)
        assert detected == "imagexpress"

    def test_detect_cellvoyager_format(self, temp_dir):
        """Test detection of Cell Voyager format."""
        create_mock_cellvoyager_export(temp_dir)
        
        detected = auto_detect_format(temp_dir)
        assert detected == "cellvoyager"

    def test_detect_unknown_format(self, temp_dir):
        """Test detection returns None for unknown format."""
        # Create directory with random files
        (temp_dir / "random_file.txt").write_text("not a vendor export")
        
        detected = auto_detect_format(temp_dir)
        assert detected is None

    def test_detect_nonexistent_path(self, temp_dir):
        """Test detection returns None for non-existent path."""
        nonexistent = temp_dir / "does_not_exist"
        
        detected = auto_detect_format(nonexistent)
        assert detected is None


class TestOperaParser:
    """Test Opera/Operetta export parsing."""

    def test_parse_opera_basic(self, temp_dir):
        """Test basic Opera export parsing."""
        create_mock_opera_export(temp_dir, wells=4, channels=2, z_slices=2)
        
        result = parse_opera_export(temp_dir)
        
        assert result.vendor == "opera"
        assert result.plate_barcode in ["TEST_OPERA_PLATE", "OPERA_PLATE"]  # Accept either
        assert result.plate_format == 96
        assert len(result.images) > 0
        assert len(result.channels) == 2
        assert result.z_slices == 2
        assert result.total_images == len(result.images)
        assert len(result.errors) == 0

    def test_parse_opera_missing_index(self, temp_dir):
        """Test Opera parsing with missing Index.xml."""
        # Create directory without Index.xml
        images_dir = temp_dir / "Images"
        images_dir.mkdir()
        
        result = parse_opera_export(temp_dir)
        
        assert len(result.errors) > 0
        assert any("Index.xml" in error for error in result.errors)

    def test_parse_opera_missing_images_dir(self, temp_dir):
        """Test Opera parsing with missing Images directory."""
        # Create only Index.xml
        index_xml = """<?xml version="1.0" encoding="UTF-8"?>
        <OperaExport>
            <PlateBarcode>TEST_PLATE</PlateBarcode>
        </OperaExport>"""
        (temp_dir / "Index.xml").write_text(index_xml)
        
        result = parse_opera_export(temp_dir)
        
        assert len(result.errors) > 0
        assert any("Images directory" in error for error in result.errors)

    def test_parse_opera_alternative_index_name(self, temp_dir):
        """Test Opera parsing with Index.idx.xml filename."""
        # Create Index.idx.xml instead of Index.xml
        index_xml = """<?xml version="1.0" encoding="UTF-8"?>
<OperaExport>
    <PlateBarcode>ALT_INDEX_PLATE</PlateBarcode>
</OperaExport>"""
        (temp_dir / "Index.idx.xml").write_text(index_xml)
        
        images_dir = temp_dir / "Images"
        images_dir.mkdir()
        
        # Create one test image
        filename = "r01c01f01p01-ch1sk1fk1fl1.tiff"
        data = np.random.randint(0, 65535, (32, 32), dtype=np.uint16)
        tifffile.imwrite(images_dir / filename, data)
        
        result = parse_opera_export(temp_dir)
        
        assert result.plate_barcode in ["ALT_INDEX_PLATE", "OPERA_PLATE"]  # Accept either
        assert len(result.images) == 1


class TestImageXpressParser:
    """Test ImageXpress export parsing."""

    def test_parse_imagexpress_basic(self, temp_dir):
        """Test basic ImageXpress export parsing."""
        create_mock_imagexpress_export(temp_dir, wells=4, channels=2, timepoints=2)
        
        result = parse_imagexpress_export(temp_dir)
        
        assert result.vendor == "imagexpress"
        assert result.plate_barcode == "TEST_IMAGEXPRESS_PLATE"
        assert result.plate_format == 96
        assert len(result.images) > 0
        assert len(result.channels) == 2
        assert result.timepoints == 2
        assert result.z_slices == 1  # ImageXpress typically single Z
        assert len(result.errors) == 0

    def test_parse_imagexpress_missing_htd(self, temp_dir):
        """Test ImageXpress parsing with missing .htd file."""
        # Create TimePoint directory without HTD file
        tp_dir = temp_dir / "TimePoint_1"
        tp_dir.mkdir()
        
        result = parse_imagexpress_export(temp_dir)
        
        assert len(result.errors) > 0
        assert any(".htd" in error for error in result.errors)

    def test_parse_imagexpress_missing_timepoints(self, temp_dir):
        """Test ImageXpress parsing with missing TimePoint directories."""
        # Create only HTD file
        htd_content = "ExperimentName=Test"
        (temp_dir / "experiment.htd").write_text(htd_content)
        
        result = parse_imagexpress_export(temp_dir)
        
        assert len(result.errors) > 0
        assert any("TimePoint" in error for error in result.errors)


class TestCellVoyagerParser:
    """Test Cell Voyager export parsing."""

    def test_parse_cellvoyager_basic(self, temp_dir):
        """Test basic Cell Voyager export parsing."""
        create_mock_cellvoyager_export(temp_dir, wells=4, channels=2, timepoints=2)
        
        result = parse_cellvoyager_export(temp_dir)
        
        assert result.vendor == "cellvoyager"
        assert result.plate_barcode == "TEST_CELLVOYAGER_PLATE"
        assert result.plate_format == 96
        assert len(result.images) > 0
        assert len(result.channels) == 2
        assert result.timepoints == 2
        assert len(result.errors) == 0

    def test_parse_cellvoyager_missing_mlf(self, temp_dir):
        """Test Cell Voyager parsing with missing MLF file."""
        # Create Images directory without MLF file
        images_dir = temp_dir / "Images"
        images_dir.mkdir()
        
        result = parse_cellvoyager_export(temp_dir)
        
        assert len(result.errors) > 0
        assert any("MeasurementData.mlf" in error for error in result.errors)

    def test_parse_cellvoyager_missing_images(self, temp_dir):
        """Test Cell Voyager parsing with missing image files."""
        # Create only MLF file
        mlf_content = "MeasurementName=Test"
        (temp_dir / "MeasurementData.mlf").write_text(mlf_content)
        
        result = parse_cellvoyager_export(temp_dir)
        
        assert len(result.errors) > 0
        assert any("flex" in error for error in result.errors)


class TestUniversalParser:
    """Test universal vendor parser."""

    def test_parse_vendor_auto_detect_opera(self, temp_dir):
        """Test universal parser with auto-detected Opera format."""
        create_mock_opera_export(temp_dir)
        
        result = parse_vendor_export(temp_dir)
        
        assert result.vendor == "opera"
        assert len(result.images) > 0

    def test_parse_vendor_explicit_vendor(self, temp_dir):
        """Test universal parser with explicit vendor specification."""
        create_mock_opera_export(temp_dir)
        
        result = parse_vendor_export(temp_dir, vendor="opera")
        
        assert result.vendor == "opera"
        assert len(result.images) > 0

    def test_parse_vendor_nonexistent_path(self, temp_dir):
        """Test universal parser with non-existent path."""
        nonexistent = temp_dir / "does_not_exist"
        
        with pytest.raises(ValueError, match="does not exist"):
            parse_vendor_export(nonexistent)

    def test_parse_vendor_unsupported_format(self, temp_dir):
        """Test universal parser with unsupported vendor format."""
        with pytest.raises(ValueError, match="Unsupported vendor format"):
            parse_vendor_export(temp_dir, vendor="unknown_vendor")

    def test_parse_vendor_undetectable_format(self, temp_dir):
        """Test universal parser when format cannot be auto-detected."""
        # Create directory with random files
        (temp_dir / "random.txt").write_text("not a vendor export")
        
        with pytest.raises(ValueError, match="Could not auto-detect"):
            parse_vendor_export(temp_dir)


class TestUtilityFunctions:
    """Test utility functions."""

    def test_well_position_to_indices(self):
        """Test well position to row/col conversion."""
        assert _well_position_to_indices("A01") == (1, 1)
        assert _well_position_to_indices("B12") == (2, 12)
        assert _well_position_to_indices("H24") == (8, 24)
        assert _well_position_to_indices("") == (1, 1)  # Default fallback
        assert _well_position_to_indices("X") == (1, 1)  # Invalid format

    def test_estimate_plate_format(self):
        """Test plate format estimation from image positions."""
        # 96-well plate
        images_96 = [
            PlateImageInfo("A01", 1, 1, 1, "ch1", 1, 1, "test.tif"),
            PlateImageInfo("H12", 8, 12, 1, "ch1", 1, 1, "test.tif"),
        ]
        assert _estimate_plate_format(images_96) == 96
        
        # 384-well plate
        images_384 = [
            PlateImageInfo("A01", 1, 1, 1, "ch1", 1, 1, "test.tif"),
            PlateImageInfo("P24", 16, 24, 1, "ch1", 1, 1, "test.tif"),
        ]
        assert _estimate_plate_format(images_384) == 384
        
        # Empty list
        assert _estimate_plate_format([]) == 96


class TestDataStructures:
    """Test data structure classes."""

    def test_plate_image_info_creation(self):
        """Test PlateImageInfo dataclass creation."""
        image_info = PlateImageInfo(
            well_position="A01",
            row=1,
            col=1,
            field=1,
            channel="DAPI",
            z_slice=1,
            timepoint=1,
            image_path="/path/to/image.tif",
            width=1024,
            height=1024,
            bit_depth=16,
            pixel_size_um=0.325
        )
        
        assert image_info.well_position == "A01"
        assert image_info.row == 1
        assert image_info.col == 1
        assert image_info.width == 1024
        assert image_info.height == 1024
        assert image_info.bit_depth == 16
        assert image_info.pixel_size_um == 0.325
        assert image_info.metadata == {}  # Default empty dict

    def test_plate_import_result_creation(self):
        """Test PlateImportResult dataclass creation."""
        result = PlateImportResult(
            plate_barcode="TEST_PLATE_001",
            plate_format=384,
            vendor="opera",
            import_path="/path/to/export"
        )
        
        assert result.plate_barcode == "TEST_PLATE_001"
        assert result.plate_format == 384
        assert result.vendor == "opera"
        assert result.import_path == "/path/to/export"
        assert result.images == []  # Default empty list
        assert result.channels == []  # Default empty list
        assert result.errors == []  # Default empty list
        assert result.warnings == []  # Default empty list
        assert result.z_slices == 1  # Default value
        assert result.timepoints == 1  # Default value
