"""OME-TIFF metadata parser using tifffile and ome-types libraries."""

from __future__ import annotations

import json
import xml.etree.ElementTree as ET
from dataclasses import dataclass, asdict
from datetime import datetime
from pathlib import Path
from typing import Dict, List, Optional, Tuple, Union

import numpy as np
import tifffile
from ome_types import OME, from_xml


@dataclass
class InstrumentInfo:
    """Microscope instrument information from OME metadata."""
    microscope_name: Optional[str] = None
    microscope_model: Optional[str] = None
    microscope_serial: Optional[str] = None
    objective_name: Optional[str] = None
    objective_magnification: Optional[float] = None
    objective_na: Optional[float] = None
    objective_immersion: Optional[str] = None


@dataclass
class ChannelInfo:
    """Channel configuration information."""
    name: str
    fluorophore: Optional[str] = None
    excitation_wavelength_nm: Optional[int] = None
    emission_wavelength_nm: Optional[int] = None
    exposure_ms: Optional[float] = None
    color: Optional[str] = None  # RGB hex


@dataclass
class ImageDimensions:
    """Image dimensions and pixel size information."""
    size_x: int
    size_y: int
    size_z: int = 1
    size_c: int = 1  # Channels
    size_t: int = 1  # Timepoints
    pixel_size_x_um: Optional[float] = None
    pixel_size_y_um: Optional[float] = None
    pixel_size_z_um: Optional[float] = None
    pixel_type: str = "uint16"  # uint8, uint16, float32


@dataclass
class OMEMetadata:
    """Complete OME-TIFF metadata structure."""
    filename: str
    ome_uuid: Optional[str] = None
    instrument: Optional[InstrumentInfo] = None
    channels: List[ChannelInfo] = None
    dimensions: Optional[ImageDimensions] = None
    acquisition_date: Optional[datetime] = None
    description: Optional[str] = None
    raw_xml: str = ""
    raw_json: Dict = None

    def __post_init__(self):
        """Initialize default values for mutable fields."""
        if self.channels is None:
            self.channels = []
        if self.raw_json is None:
            self.raw_json = {}


def parse_ome_tiff(path: Union[str, Path]) -> OMEMetadata:
    """
    Parse OME-TIFF file and extract all metadata.
    
    Args:
        path: Path to OME-TIFF file
        
    Returns:
        OMEMetadata object with parsed information
        
    Raises:
        FileNotFoundError: If file doesn't exist
        ValueError: If file is not a valid OME-TIFF
    """
    path = Path(path)
    
    if not path.exists():
        raise FileNotFoundError(f"File not found: {path}")
    
    try:
        # Read TIFF file and extract OME-XML
        with tifffile.TiffFile(path) as tif:
            if not tif.ome_metadata:
                raise ValueError(f"No OME metadata found in {path}")
            
            ome_xml = tif.ome_metadata
            
            # Parse OME-XML using ome-types
            try:
                ome = from_xml(ome_xml)
            except Exception as e:
                # If ome-types fails, create minimal metadata
                return _create_fallback_metadata(path, ome_xml, tif)
            
            # Extract metadata components
            instrument = extract_instrument_info(ome)
            channels = extract_channels(ome)
            dimensions = extract_dimensions(ome, tif)
            
            # Extract acquisition date
            acquisition_date = None
            if ome.images and ome.images[0].acquisition_date:
                acquisition_date = ome.images[0].acquisition_date
            
            # Extract description
            description = None
            if ome.images and ome.images[0].description:
                description = ome.images[0].description
            
            # Extract UUID
            ome_uuid = None
            if ome.uuid:
                ome_uuid = str(ome.uuid)
            
            return OMEMetadata(
                filename=path.name,
                ome_uuid=ome_uuid,
                instrument=instrument,
                channels=channels,
                dimensions=dimensions,
                acquisition_date=acquisition_date,
                description=description,
                raw_xml=ome_xml,
                raw_json=ome_to_json(ome)
            )
            
    except Exception as e:
        raise ValueError(f"Failed to parse OME-TIFF {path}: {e}")


def extract_instrument_info(ome: OME) -> Optional[InstrumentInfo]:
    """
    Extract microscope and objective information from OME model.
    
    Args:
        ome: Parsed OME model
        
    Returns:
        InstrumentInfo object or None if no instrument data
    """
    if not ome or not ome.instruments:
        return None
    
    instrument = ome.instruments[0]  # Use first instrument
    
    # Extract microscope info
    microscope_name = None
    microscope_model = None
    microscope_serial = None
    
    if instrument.microscope:
        microscope_name = getattr(instrument.microscope, 'manufacturer', None)
        microscope_model = getattr(instrument.microscope, 'model', None)
        microscope_serial = getattr(instrument.microscope, 'serial_number', None)
    
    # Extract objective info
    objective_name = None
    objective_magnification = None
    objective_na = None
    objective_immersion = None
    
    if instrument.objectives:
        obj = instrument.objectives[0]  # Use first objective
        objective_name = getattr(obj, 'manufacturer', None) or getattr(obj, 'model', None)
        objective_magnification = getattr(obj, 'nominal_magnification', None)
        objective_na = getattr(obj, 'lens_na', None)
        objective_immersion = getattr(obj, 'immersion', None)
        if objective_immersion and hasattr(objective_immersion, 'value'):
            objective_immersion = objective_immersion.value
    
    # Only return if we have some instrument data
    if any([microscope_name, microscope_model, microscope_serial, 
            objective_name, objective_magnification, objective_na]):
        return InstrumentInfo(
            microscope_name=microscope_name,
            microscope_model=microscope_model,
            microscope_serial=microscope_serial,
            objective_name=objective_name,
            objective_magnification=objective_magnification,
            objective_na=objective_na,
            objective_immersion=objective_immersion
        )
    
    return None


def extract_channels(ome: OME) -> List[ChannelInfo]:
    """
    Extract channel configurations from OME model.
    
    Args:
        ome: Parsed OME model
        
    Returns:
        List of ChannelInfo objects
    """
    channels = []
    
    if not ome or not ome.images:
        return channels
    
    image = ome.images[0]  # Use first image
    if not image.pixels or not image.pixels.channels:
        return channels
    
    for i, channel in enumerate(image.pixels.channels):
        # Extract basic channel info
        name = channel.name or f"Channel_{i}"
        fluorophore = getattr(channel, 'fluor', None)
        
        # Extract wavelengths
        excitation_wavelength = None
        emission_wavelength = None
        
        if hasattr(channel, 'excitation_wavelength') and channel.excitation_wavelength:
            try:
                excitation_wavelength = int(channel.excitation_wavelength.value)
            except AttributeError:
                excitation_wavelength = int(channel.excitation_wavelength)
        
        if hasattr(channel, 'emission_wavelength') and channel.emission_wavelength:
            try:
                emission_wavelength = int(channel.emission_wavelength.value)
            except AttributeError:
                emission_wavelength = int(channel.emission_wavelength)
        
        # Extract color (convert to hex if available)
        color = None
        if hasattr(channel, 'color') and channel.color:
            try:
                # Handle ome-types Color object
                if hasattr(channel.color, 'value'):
                    color_val = channel.color.value
                else:
                    color_val = int(channel.color)
                color = f"#{color_val:06x}"
            except (ValueError, TypeError, AttributeError):
                # If color conversion fails, leave as None
                color = None
        
        channels.append(ChannelInfo(
            name=name,
            fluorophore=fluorophore,
            excitation_wavelength_nm=excitation_wavelength,
            emission_wavelength_nm=emission_wavelength,
            color=color
        ))
    
    return channels


def extract_dimensions(ome: OME, tif: tifffile.TiffFile = None) -> ImageDimensions:
    """
    Extract image dimensions and pixel sizes.
    
    Args:
        ome: Parsed OME model
        tif: Optional TiffFile for fallback dimension extraction
        
    Returns:
        ImageDimensions object
    """
    if ome and ome.images and ome.images[0].pixels:
        pixels = ome.images[0].pixels
        
        size_x = pixels.size_x or 1
        size_y = pixels.size_y or 1
        size_z = pixels.size_z or 1
        size_c = pixels.size_c or 1
        size_t = pixels.size_t or 1
        
        # Extract pixel sizes
        pixel_size_x = None
        pixel_size_y = None
        pixel_size_z = None
        
        if pixels.physical_size_x:
            try:
                pixel_size_x = float(pixels.physical_size_x.value)
            except AttributeError:
                pixel_size_x = float(pixels.physical_size_x)
        if pixels.physical_size_y:
            try:
                pixel_size_y = float(pixels.physical_size_y.value)
            except AttributeError:
                pixel_size_y = float(pixels.physical_size_y)
        if pixels.physical_size_z:
            try:
                pixel_size_z = float(pixels.physical_size_z.value)
            except AttributeError:
                pixel_size_z = float(pixels.physical_size_z)
        
        # Extract pixel type
        pixel_type = "uint16"  # Default
        if pixels.type:
            pixel_type = str(pixels.type.value).lower()
    
    elif tif:
        # Fallback to TiffFile dimensions
        page = tif.pages[0]
        size_x = page.imagewidth
        size_y = page.imagelength
        size_z = 1
        size_c = 1
        size_t = 1
        
        # Try to determine from shape if multi-dimensional
        if len(tif.series) > 0:
            series = tif.series[0]
            shape = series.shape
            if len(shape) >= 4:  # CZYX or similar
                size_c = shape[0] if len(shape) >= 4 else 1
                size_z = shape[1] if len(shape) >= 4 else 1
        
        pixel_size_x = pixel_size_y = pixel_size_z = None
        pixel_type = str(page.dtype)
    
    else:
        # Minimal fallback
        size_x = size_y = size_z = size_c = size_t = 1
        pixel_size_x = pixel_size_y = pixel_size_z = None
        pixel_type = "uint16"
    
    return ImageDimensions(
        size_x=size_x,
        size_y=size_y,
        size_z=size_z,
        size_c=size_c,
        size_t=size_t,
        pixel_size_x_um=pixel_size_x,
        pixel_size_y_um=pixel_size_y,
        pixel_size_z_um=pixel_size_z,
        pixel_type=pixel_type
    )


def validate_ome_schema(xml_string: str) -> Tuple[bool, Optional[str]]:
    """
    Validate OME-XML against schema.
    
    Args:
        xml_string: OME-XML content
        
    Returns:
        Tuple of (is_valid, error_message)
    """
    try:
        # Try to parse with ome-types
        from_xml(xml_string)
        return True, None
    except Exception as e:
        return False, str(e)


def ome_to_json(ome: OME) -> Dict:
    """
    Convert OME model to JSON-serializable dict for DB storage.
    
    Args:
        ome: Parsed OME model
        
    Returns:
        JSON-serializable dictionary
    """
    try:
        # Convert to dict using pydantic's model_dump if available
        if hasattr(ome, 'model_dump'):
            return ome.model_dump(exclude_none=True)
        elif hasattr(ome, 'dict'):
            return ome.dict(exclude_none=True)
        else:
            # Fallback: extract key information manually
            result = {}
            
            if ome.uuid:
                result['uuid'] = str(ome.uuid)
            
            if ome.images:
                result['images'] = []
                for img in ome.images:
                    img_dict = {'id': img.id}
                    if img.name:
                        img_dict['name'] = img.name
                    if img.acquisition_date:
                        img_dict['acquisition_date'] = img.acquisition_date.isoformat()
                    result['images'].append(img_dict)
            
            if ome.instruments:
                result['instruments'] = []
                for inst in ome.instruments:
                    inst_dict = {'id': inst.id}
                    if inst.microscope:
                        inst_dict['microscope'] = {}
                        if hasattr(inst.microscope, 'manufacturer'):
                            inst_dict['microscope']['manufacturer'] = inst.microscope.manufacturer
                        if hasattr(inst.microscope, 'model'):
                            inst_dict['microscope']['model'] = inst.microscope.model
                    result['instruments'].append(inst_dict)
            
            return result
            
    except Exception:
        # Ultimate fallback
        return {'error': 'Failed to serialize OME model'}


def read_ome_tiff_data(path: Union[str, Path], 
                      c: int = 0, z: int = 0, t: int = 0) -> np.ndarray:
    """
    Read specific plane from OME-TIFF (for large files).
    
    Args:
        path: Path to OME-TIFF file
        c: Channel index
        z: Z-slice index
        t: Timepoint index
        
    Returns:
        2D numpy array for the specified plane
        
    Raises:
        ValueError: If indices are out of bounds
        FileNotFoundError: If file doesn't exist
    """
    path = Path(path)
    
    if not path.exists():
        raise FileNotFoundError(f"File not found: {path}")
    
    try:
        with tifffile.TiffFile(path) as tif:
            if len(tif.series) == 0:
                raise ValueError("No image series found in file")
            
            series = tif.series[0]
            
            # Handle different dimensionalities
            if series.shape == (series.shape[-2], series.shape[-1]):
                # 2D image (YX)
                if c != 0 or z != 0 or t != 0:
                    raise ValueError("Indices out of bounds for 2D image")
                return series.asarray()
            
            elif len(series.shape) == 3:
                # 3D image - could be CYX, ZYX, or TYX
                if c != 0 and z == 0 and t == 0:
                    # Assume CYX
                    if c >= series.shape[0]:
                        raise ValueError(f"Channel index {c} out of bounds (max {series.shape[0] - 1})")
                    return series.asarray()[c]
                elif z != 0 and c == 0 and t == 0:
                    # Assume ZYX
                    if z >= series.shape[0]:
                        raise ValueError(f"Z index {z} out of bounds (max {series.shape[0] - 1})")
                    return series.asarray()[z]
                else:
                    # Default to first plane
                    return series.asarray()[0]
            
            elif len(series.shape) >= 4:
                # Multi-dimensional (CZYX, TZYX, CTZYX, etc.)
                data = series.asarray()
                
                # Try to extract the requested plane
                # This is a simplified approach - real implementation would need
                # to understand the exact axis order from OME metadata
                try:
                    if len(series.shape) == 4:  # CZYX or TZYX
                        return data[c, z]  # Assume CZYX
                    elif len(series.shape) == 5:  # CTZYX
                        return data[c, t, z]
                    else:
                        # Fallback to first available plane
                        return data[0, 0] if len(data.shape) >= 4 else data[0]
                        
                except IndexError:
                    raise ValueError(f"Indices ({c}, {z}, {t}) out of bounds for shape {series.shape}")
            
            else:
                raise ValueError(f"Unsupported image dimensionality: {series.shape}")
                
    except Exception as e:
        if isinstance(e, (ValueError, FileNotFoundError)):
            raise
        else:
            raise ValueError(f"Failed to read OME-TIFF data: {e}")


def _create_fallback_metadata(path: Path, ome_xml: str, tif: tifffile.TiffFile) -> OMEMetadata:
    """Create minimal metadata when ome-types parsing fails."""
    # Extract basic dimensions from TiffFile
    page = tif.pages[0]
    dimensions = ImageDimensions(
        size_x=page.imagewidth,
        size_y=page.imagelength,
        size_z=1,
        size_c=1,
        size_t=1,
        pixel_type=str(page.dtype)
    )
    
    # Try to extract basic info from XML manually
    channels = []
    try:
        root = ET.fromstring(ome_xml)
        # Look for channel names in XML
        for channel_elem in root.findall('.//{*}Channel'):
            name = channel_elem.get('Name', f'Channel_{len(channels)}')
            channels.append(ChannelInfo(name=name))
    except Exception:
        # If XML parsing fails, create default channel
        channels = [ChannelInfo(name="Channel_0")]
    
    return OMEMetadata(
        filename=path.name,
        channels=channels,
        dimensions=dimensions,
        raw_xml=ome_xml,
        raw_json={'fallback': True, 'error': 'ome-types parsing failed'}
    )
