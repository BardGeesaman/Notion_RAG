"""Vendor-specific parsers for HCS instrument export formats."""

from __future__ import annotations

import re
import xml.etree.ElementTree as ET
from dataclasses import dataclass, field
from datetime import datetime
from pathlib import Path
from typing import Any, Dict, List, Optional, Union

import tifffile


@dataclass
class PlateImageInfo:
    """Information about a single plate image."""
    well_position: str  # "A01", "B12", etc.
    row: int
    col: int
    field: int  # Field of view index
    channel: str
    z_slice: int
    timepoint: int
    image_path: str
    thumbnail_path: Optional[str] = None
    width: int = 0
    height: int = 0
    bit_depth: int = 16
    pixel_size_um: Optional[float] = None
    exposure_ms: Optional[float] = None
    acquired_at: Optional[datetime] = None
    metadata: Dict[str, Any] = field(default_factory=dict)


@dataclass
class PlateImportResult:
    """Result of parsing a vendor plate export."""
    plate_barcode: str
    plate_format: int  # 96, 384, 1536
    vendor: str  # "opera", "imagexpress", "cellvoyager"
    import_path: str
    images: List[PlateImageInfo] = field(default_factory=list)
    channels: List[str] = field(default_factory=list)
    z_slices: int = 1
    timepoints: int = 1
    wells_with_images: int = 0
    total_images: int = 0
    errors: List[str] = field(default_factory=list)
    warnings: List[str] = field(default_factory=list)


def auto_detect_format(path: Union[str, Path]) -> Optional[str]:
    """
    Auto-detect vendor format from directory structure.
    
    Args:
        path: Path to export directory
        
    Returns:
        Vendor format ("opera", "imagexpress", "cellvoyager") or None
    """
    path = Path(path)
    
    if not path.exists() or not path.is_dir():
        return None
    
    # Check for Opera/Operetta format
    if (path / "Index.xml").exists() or (path / "Index.idx.xml").exists():
        if (path / "Images").exists():
            return "opera"
    
    # Check for ImageXpress format
    htd_files = list(path.glob("*.htd"))
    if htd_files and any((path / f"TimePoint_{i}").exists() for i in range(1, 5)):
        return "imagexpress"
    
    # Check for Cell Voyager format
    if (path / "MeasurementData.mlf").exists():
        if (path / "Images").exists() or list(path.glob("*.flex")):
            return "cellvoyager"
    
    return None


def parse_opera_export(path: Union[str, Path]) -> PlateImportResult:
    """
    Parse PerkinElmer Opera/Operetta export.
    
    Args:
        path: Path to export directory
        
    Returns:
        PlateImportResult with parsed metadata and image list
    """
    path = Path(path)
    result = PlateImportResult(
        plate_barcode="",  # Will be set from XML
        plate_format=96,  # Default, will be updated
        vendor="opera",
        import_path=str(path)
    )
    
    # Look for Index.xml or Index.idx.xml
    index_file = None
    if (path / "Index.xml").exists():
        index_file = path / "Index.xml"
    elif (path / "Index.idx.xml").exists():
        index_file = path / "Index.idx.xml"
    else:
        result.errors.append("No Index.xml or Index.idx.xml file found")
        return result
    
    # Parse index file for plate metadata
    try:
        tree = ET.parse(index_file)
        root = tree.getroot()
        
        # Extract plate barcode if available
        barcode_elem = root.find(".//PlateBarcode") or root.find(".//Barcode")
        if barcode_elem is not None and barcode_elem.text:
            result.plate_barcode = barcode_elem.text.strip()
        
        # If no barcode found, set default
        if not result.plate_barcode:
            result.plate_barcode = "OPERA_PLATE"
        
        # Extract plate format
        plate_type_elem = root.find(".//PlateType") or root.find(".//PlateFormat")
        if plate_type_elem is not None and plate_type_elem.text:
            try:
                result.plate_format = int(plate_type_elem.text.strip())
            except ValueError:
                result.warnings.append(f"Could not parse plate format: {plate_type_elem.text}")
                
    except ET.ParseError as e:
        result.errors.append(f"Failed to parse index file: {e}")
        result.plate_barcode = "OPERA_PLATE"  # Set default on error
        return result
    
    # Look for Images directory
    images_dir = path / "Images"
    if not images_dir.exists():
        result.errors.append("Images directory not found")
        return result
    
    # Parse image files
    # Opera filename pattern: r{row}c{col}f{field}p{plane}-ch{channel}sk{z}fk{frame}fl{timeline}.tiff
    image_pattern = re.compile(
        r'r(\d+)c(\d+)f(\d+)p(\d+)-ch(\d+)sk(\d+)fk(\d+)fl(\d+)\.tiff?',
        re.IGNORECASE
    )
    
    channels_set = set()
    max_z = 0
    max_timepoint = 0
    wells_set = set()
    
    # Get all TIFF files (case insensitive)
    tiff_files = []
    for pattern in ["*.tif", "*.tiff", "*.TIF", "*.TIFF"]:
        tiff_files.extend(images_dir.glob(pattern))
    
    for image_file in tiff_files:
        match = image_pattern.match(image_file.name)
        if not match:
            result.warnings.append(f"Could not parse filename: {image_file.name}")
            continue
        
        row, col, field, plane, channel, z_slice, frame, timeline = map(int, match.groups())
        
        # Convert row/col to well position
        well_position = f"{chr(ord('A') + row - 1)}{col:02d}"
        
        # Extract image dimensions if possible
        width, height, bit_depth = 0, 0, 16
        try:
            with tifffile.TiffFile(image_file) as tif:
                page = tif.pages[0]
                width = page.imagewidth
                height = page.imagelength
                bit_depth = page.bitspersample
        except Exception as e:
            result.warnings.append(f"Could not read image dimensions for {image_file.name}: {e}")
        
        image_info = PlateImageInfo(
            well_position=well_position,
            row=row,
            col=col,
            field=field,
            channel=f"ch{channel}",
            z_slice=z_slice,
            timepoint=timeline,
            image_path=str(image_file),
            width=width,
            height=height,
            bit_depth=bit_depth,
            metadata={
                "plane": plane,
                "frame": frame,
                "original_filename": image_file.name
            }
        )
        
        result.images.append(image_info)
        channels_set.add(f"ch{channel}")
        max_z = max(max_z, z_slice)
        max_timepoint = max(max_timepoint, timeline)
        wells_set.add(well_position)
    
    # Update summary statistics
    result.channels = sorted(channels_set)
    result.z_slices = max_z
    result.timepoints = max_timepoint
    result.wells_with_images = len(wells_set)
    result.total_images = len(result.images)
    
    if result.total_images == 0:
        result.errors.append("No valid images found")
    
    return result


def parse_imagexpress_export(path: Union[str, Path]) -> PlateImportResult:
    """
    Parse Molecular Devices ImageXpress export.
    
    Args:
        path: Path to export directory
        
    Returns:
        PlateImportResult with parsed metadata and image list
    """
    path = Path(path)
    result = PlateImportResult(
        plate_barcode="IMAGEXPRESS_PLATE",
        plate_format=96,
        vendor="imagexpress",
        import_path=str(path)
    )
    
    # Look for .htd experiment file
    htd_files = list(path.glob("*.htd"))
    if not htd_files:
        result.errors.append("No .htd experiment file found")
        return result
    
    htd_file = htd_files[0]
    
    # Parse HTD file for experiment metadata
    try:
        # HTD files are typically text-based with key-value pairs
        with open(htd_file, 'r', encoding='utf-8', errors='ignore') as f:
            content = f.read()
            
        # Extract plate barcode/name
        barcode_match = re.search(r'PlateName[=:]\s*([^\n\r]+)', content, re.IGNORECASE)
        if barcode_match:
            result.plate_barcode = barcode_match.group(1).strip()
        
        # Extract plate format
        format_match = re.search(r'PlateFormat[=:]\s*(\d+)', content, re.IGNORECASE)
        if format_match:
            result.plate_format = int(format_match.group(1))
            
    except Exception as e:
        result.warnings.append(f"Could not parse HTD file: {e}")
    
    # Look for TimePoint directories
    timepoint_dirs = sorted([d for d in path.iterdir() 
                           if d.is_dir() and d.name.startswith("TimePoint_")])
    
    if not timepoint_dirs:
        result.errors.append("No TimePoint directories found")
        return result
    
    # Parse images from TimePoint directories
    # ImageXpress pattern: {well}_{site}_w{wavelength}.TIF
    image_pattern = re.compile(r'([A-Z]\d+)_s(\d+)_w(\d+)\.tiff?', re.IGNORECASE)
    
    channels_set = set()
    max_timepoint = 0
    wells_set = set()
    
    for tp_idx, tp_dir in enumerate(timepoint_dirs, 1):
        timepoint_match = re.search(r'TimePoint_(\d+)', tp_dir.name)
        timepoint = int(timepoint_match.group(1)) if timepoint_match else tp_idx
        max_timepoint = max(max_timepoint, timepoint)
        
        # Get all TIFF files (case insensitive)
        tiff_files = []
        for pattern in ["*.tif", "*.tiff", "*.TIF", "*.TIFF"]:
            tiff_files.extend(tp_dir.glob(pattern))
        
        for image_file in tiff_files:
            match = image_pattern.match(image_file.name)
            if not match:
                result.warnings.append(f"Could not parse filename: {image_file.name}")
                continue
            
            well_position, site, wavelength = match.groups()
            site = int(site)
            wavelength = int(wavelength)
            
            # Convert well position to row/col
            row = ord(well_position[0]) - ord('A') + 1
            col = int(well_position[1:])
            
            # Extract image dimensions
            width, height, bit_depth = 0, 0, 16
            try:
                with tifffile.TiffFile(image_file) as tif:
                    page = tif.pages[0]
                    width = page.imagewidth
                    height = page.imagelength
                    bit_depth = page.bitspersample
            except Exception as e:
                result.warnings.append(f"Could not read image dimensions for {image_file.name}: {e}")
            
            image_info = PlateImageInfo(
                well_position=well_position,
                row=row,
                col=col,
                field=site,
                channel=f"w{wavelength}",
                z_slice=1,  # ImageXpress typically single Z
                timepoint=timepoint,
                image_path=str(image_file),
                width=width,
                height=height,
                bit_depth=bit_depth,
                metadata={
                    "wavelength": wavelength,
                    "original_filename": image_file.name
                }
            )
            
            result.images.append(image_info)
            channels_set.add(f"w{wavelength}")
            wells_set.add(well_position)
    
    # Update summary statistics
    result.channels = sorted(channels_set)
    result.z_slices = 1  # Typically single Z for ImageXpress
    result.timepoints = max_timepoint
    result.wells_with_images = len(wells_set)
    result.total_images = len(result.images)
    
    if result.total_images == 0:
        result.errors.append("No valid images found")
    
    return result


def parse_cellvoyager_export(path: Union[str, Path]) -> PlateImportResult:
    """
    Parse Yokogawa Cell Voyager export.
    
    Args:
        path: Path to export directory
        
    Returns:
        PlateImportResult with parsed metadata and image list
    """
    path = Path(path)
    result = PlateImportResult(
        plate_barcode="CELLVOYAGER_PLATE",
        plate_format=96,
        vendor="cellvoyager",
        import_path=str(path)
    )
    
    # Look for MeasurementData.mlf
    mlf_file = path / "MeasurementData.mlf"
    if not mlf_file.exists():
        result.errors.append("MeasurementData.mlf file not found")
        return result
    
    # Parse MLF file for experiment metadata
    try:
        with open(mlf_file, 'r', encoding='utf-8', errors='ignore') as f:
            content = f.read()
            
        # Extract plate information
        barcode_match = re.search(r'PlateName[=:]\s*([^\n\r]+)', content, re.IGNORECASE)
        if barcode_match:
            result.plate_barcode = barcode_match.group(1).strip()
            
        format_match = re.search(r'PlateFormat[=:]\s*(\d+)', content, re.IGNORECASE)
        if format_match:
            result.plate_format = int(format_match.group(1))
            
    except Exception as e:
        result.warnings.append(f"Could not parse MLF file: {e}")
    
    # Look for image files
    # Cell Voyager pattern: W{well}F{field}T{time}Z{z}C{channel}.flex
    image_pattern = re.compile(r'W(\d+)F(\d+)T(\d+)Z(\d+)C(\d+)\.flex', re.IGNORECASE)
    
    # Check Images directory first, then root directory
    search_dirs = [path / "Images", path]
    image_files = []
    
    for search_dir in search_dirs:
        if search_dir.exists():
            image_files.extend(search_dir.glob("*.flex"))
            break
    
    if not image_files:
        result.errors.append("No .flex image files found")
        return result
    
    channels_set = set()
    max_z = 0
    max_timepoint = 0
    wells_set = set()
    
    for image_file in image_files:
        match = image_pattern.match(image_file.name)
        if not match:
            result.warnings.append(f"Could not parse filename: {image_file.name}")
            continue
        
        well_num, field, timepoint, z_slice, channel = map(int, match.groups())
        
        # Convert well number to position (assuming 96-well plate)
        row = (well_num - 1) // 12 + 1
        col = (well_num - 1) % 12 + 1
        well_position = f"{chr(ord('A') + row - 1)}{col:02d}"
        
        # Note: .flex files are proprietary format, dimensions would need special reader
        width, height, bit_depth = 0, 0, 16
        
        image_info = PlateImageInfo(
            well_position=well_position,
            row=row,
            col=col,
            field=field,
            channel=f"C{channel}",
            z_slice=z_slice,
            timepoint=timepoint,
            image_path=str(image_file),
            width=width,
            height=height,
            bit_depth=bit_depth,
            metadata={
                "well_number": well_num,
                "original_filename": image_file.name
            }
        )
        
        result.images.append(image_info)
        channels_set.add(f"C{channel}")
        max_z = max(max_z, z_slice)
        max_timepoint = max(max_timepoint, timepoint)
        wells_set.add(well_position)
    
    # Update summary statistics
    result.channels = sorted(channels_set)
    result.z_slices = max_z
    result.timepoints = max_timepoint
    result.wells_with_images = len(wells_set)
    result.total_images = len(result.images)
    
    if result.total_images == 0:
        result.errors.append("No valid images found")
    
    return result


def parse_vendor_export(path: Union[str, Path], vendor: Optional[str] = None) -> PlateImportResult:
    """
    Universal parser - auto-detects format or uses specified vendor.
    
    Args:
        path: Path to export directory
        vendor: Optional vendor override ("opera", "imagexpress", "cellvoyager")
        
    Returns:
        PlateImportResult with parsed metadata and image list
        
    Raises:
        ValueError: If vendor format cannot be determined or is unsupported
    """
    path = Path(path)
    
    if not path.exists():
        raise ValueError(f"Export path does not exist: {path}")
    
    if not path.is_dir():
        raise ValueError(f"Export path is not a directory: {path}")
    
    # Auto-detect vendor if not specified
    if vendor is None:
        vendor = auto_detect_format(path)
        if vendor is None:
            raise ValueError(f"Could not auto-detect vendor format for: {path}")
    
    # Route to appropriate parser
    if vendor.lower() == "opera":
        return parse_opera_export(path)
    elif vendor.lower() == "imagexpress":
        return parse_imagexpress_export(path)
    elif vendor.lower() == "cellvoyager":
        return parse_cellvoyager_export(path)
    else:
        raise ValueError(f"Unsupported vendor format: {vendor}")


def _well_position_to_indices(well_position: str) -> tuple[int, int]:
    """Convert well position like 'A01' to row, col indices (1-based)."""
    if not well_position or len(well_position) < 2:
        return 1, 1
    
    row_letter = well_position[0].upper()
    col_str = well_position[1:]
    
    try:
        row = ord(row_letter) - ord('A') + 1
        col = int(col_str)
        return row, col
    except (ValueError, IndexError):
        return 1, 1


def _estimate_plate_format(images: List[PlateImageInfo]) -> int:
    """Estimate plate format (96, 384, 1536) from image well positions."""
    if not images:
        return 96
    
    max_row = 0
    max_col = 0
    
    for image in images:
        row, col = _well_position_to_indices(image.well_position)
        max_row = max(max_row, row)
        max_col = max(max_col, col)
    
    # Determine plate format based on dimensions
    if max_row <= 8 and max_col <= 12:
        return 96
    elif max_row <= 16 and max_col <= 24:
        return 384
    elif max_row <= 32 and max_col <= 48:
        return 1536
    else:
        return 96  # Default fallback
