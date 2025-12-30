"""
SPR (Surface Plasmon Resonance) data parsers for Biacore instruments.

Supports parsing of Biacore CSV and TXT sensorgram files, with automatic
detection of injection phases and validation of array sizes.
"""

from __future__ import annotations

import csv
import logging
import re
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple, Union

import numpy as np

logger = logging.getLogger(__name__)

# Maximum data points per sensorgram for performance
MAX_ARRAY_SIZE = 50000


@dataclass
class Sensorgram:
    """Single sensorgram data from SPR experiment."""
    
    cycle: int
    concentration: float  # Analyte concentration (nM)
    time: np.ndarray  # Time points (seconds)
    response: np.ndarray  # Response values (RU)
    association_start: float  # Start time of association phase (s)
    dissociation_start: float  # Start time of dissociation phase (s)
    flow_cell: int = 1  # Flow cell number
    reference_subtracted: bool = False  # Whether reference was subtracted


@dataclass
class SPRData:
    """Complete SPR dataset from instrument file."""
    
    instrument: str  # Instrument model (e.g., "Biacore T200")
    chip_type: str  # Chip type (e.g., "CM5", "NTA")
    temperature: float  # Experiment temperature (°C)
    sensorgrams: List[Sensorgram]  # List of sensorgram cycles
    metadata: Dict[str, Any]  # Additional metadata from file
    flow_rate: Optional[float] = None  # Flow rate (μL/min)
    buffer: Optional[str] = None  # Running buffer composition
    ligand_name: Optional[str] = None  # Immobilized ligand name
    analyte_name: Optional[str] = None  # Injected analyte name


def validate_array_size(arr: np.ndarray, max_points: int = MAX_ARRAY_SIZE) -> None:
    """
    Validate that array doesn't exceed maximum size for database storage.
    
    Args:
        arr: Array to validate
        max_points: Maximum allowed array size
        
    Raises:
        ValueError: If array exceeds maximum size
    """
    if len(arr) > max_points:
        raise ValueError(
            f"Array size {len(arr)} exceeds maximum {max_points} points. "
            f"Consider downsampling the data."
        )


def detect_injection_phases(time: np.ndarray, response: np.ndarray) -> Tuple[float, float]:
    """
    Detect association and dissociation phase start times from sensorgram.
    
    Uses derivative analysis to find rapid changes in response that indicate
    injection start and end.
    
    Args:
        time: Time array (seconds)
        response: Response array (RU)
        
    Returns:
        Tuple of (association_start, dissociation_start) times in seconds
        
    Raises:
        ValueError: If phases cannot be detected
    """
    if len(time) != len(response) or len(time) < 10:
        raise ValueError("Invalid time/response arrays for phase detection")
    
    # Calculate first derivative to find rapid changes
    dt = np.diff(time)
    dr = np.diff(response)
    
    # Avoid division by zero
    dt[dt == 0] = np.finfo(float).eps
    derivative = dr / dt
    
    # Smooth derivative to reduce noise
    if len(derivative) >= 5:
        try:
            from scipy.ndimage import gaussian_filter1d
            derivative_smooth = gaussian_filter1d(derivative, sigma=2)
        except ImportError:
            # Fallback to simple moving average if scipy not available
            window = min(5, len(derivative))
            derivative_smooth = np.convolve(derivative, np.ones(window)/window, mode='same')
    else:
        derivative_smooth = derivative
    
    # Find association start (first significant positive derivative)
    association_start = 0.0
    pos_threshold = np.std(derivative_smooth) * 2
    pos_indices = np.where(derivative_smooth > pos_threshold)[0]
    if len(pos_indices) > 0:
        association_start = time[pos_indices[0]]
    
    # Find dissociation start (first significant negative derivative after association)
    dissociation_start = time[-1]  # Default to end if not found
    neg_threshold = -np.std(derivative_smooth) * 2
    
    # Look for negative derivative after association start
    assoc_idx = np.searchsorted(time, association_start)
    if assoc_idx < len(time) - 1:
        neg_indices = np.where(derivative_smooth[assoc_idx:] < neg_threshold)[0]
        if len(neg_indices) > 0:
            dissociation_start = time[assoc_idx + neg_indices[0]]
    
    logger.debug(f"Detected phases: association={association_start:.1f}s, dissociation={dissociation_start:.1f}s")
    
    return association_start, dissociation_start


def parse_biacore_csv(path: Union[str, Path]) -> SPRData:
    """
    Parse Biacore CSV export file.
    
    Expected format:
    - Tab-separated values
    - Header rows with metadata (instrument, chip type, etc.)
    - Data columns: Time, Response, Concentration, Cycle
    
    Args:
        path: Path to Biacore CSV file
        
    Returns:
        SPRData object with parsed sensorgrams
        
    Raises:
        FileNotFoundError: If file doesn't exist
        ValueError: If file format is invalid
    """
    path = Path(path)
    if not path.exists():
        raise FileNotFoundError(f"SPR file not found: {path}")
    
    logger.info(f"Parsing Biacore CSV: {path}")
    
    metadata = {}
    sensorgrams = []
    data_started = False
    
    try:
        with open(path, 'r', encoding='utf-8') as f:
            # First pass: extract metadata from header
            lines = f.readlines()
            
        header_lines = []
        data_lines = []
        
        for i, line in enumerate(lines):
            line = line.strip()
            if not line:
                continue
                
            # Check if this looks like a header line (contains column names)
            if not data_started and ('time' in line.lower() or 'response' in line.lower()):
                data_started = True
                data_lines.append(line)
            # Check if this looks like a data line (starts with number)
            elif re.match(r'^[\d\-\.]', line) and ('\t' in line or ',' in line):
                if not data_started:
                    # No header found, assume first data line and create default header
                    parts = line.split('\t') if '\t' in line else line.split(',')
                    if len(parts) >= 2:
                        data_lines.append("Time\tResponse\tConcentration\tCycle")  # Default header
                        data_started = True
                data_lines.append(line)
            elif not data_started:
                header_lines.append(line)
            else:
                data_lines.append(line)
        
        # Parse metadata from header
        instrument = "Biacore"
        chip_type = "CM5"
        temperature = 25.0
        flow_rate = None
        buffer = None
        ligand_name = None
        analyte_name = None
        
        for line in header_lines:
            if 'instrument' in line.lower():
                match = re.search(r'instrument[:\s]+([^,\n\r]+)', line, re.IGNORECASE)
                if match:
                    instrument = match.group(1).strip()
            elif 'chip' in line.lower():
                match = re.search(r'chip[:\s]+([^,\n\r]+)', line, re.IGNORECASE)
                if match:
                    chip_type = match.group(1).strip()
            elif 'temperature' in line.lower():
                match = re.search(r'temperature[:\s]+([\d\.]+)', line, re.IGNORECASE)
                if match:
                    temperature = float(match.group(1))
            elif 'flow' in line.lower():
                match = re.search(r'flow[^:]*:\s*([\d\.]+)', line, re.IGNORECASE)
                if match:
                    flow_rate = float(match.group(1))
            elif 'buffer' in line.lower():
                match = re.search(r'buffer[:\s]+([^,\n\r]+)', line, re.IGNORECASE)
                if match:
                    buffer = match.group(1).strip()
            elif 'ligand' in line.lower():
                match = re.search(r'ligand[:\s]+([^,\n\r]+)', line, re.IGNORECASE)
                if match:
                    ligand_name = match.group(1).strip()
            elif 'analyte' in line.lower():
                match = re.search(r'analyte[:\s]+([^,\n\r]+)', line, re.IGNORECASE)
                if match:
                    analyte_name = match.group(1).strip()
        
        # Parse data section
        if not data_lines:
            raise ValueError("No data section found in CSV file")
        
        # Parse CSV data
        csv_data = []
        for line in data_lines:
            if line.strip() and not line.startswith('#'):
                csv_data.append(line.split('\t'))
        
        if len(csv_data) < 2:
            raise ValueError("Insufficient data in CSV file")
        
        # Determine column indices (flexible header detection)
        header = csv_data[0]
        time_col = None
        response_col = None
        conc_col = None
        cycle_col = None
        
        logger.debug(f"CSV header: {header}")
        
        for i, col in enumerate(header):
            col_lower = col.lower().strip()
            logger.debug(f"Column {i}: '{col}' -> '{col_lower}'")
            if 'time' in col_lower or col_lower == 'time':
                time_col = i
                logger.debug(f"Found time column at index {i}")
            elif 'response' in col_lower or 'ru' in col_lower or col_lower == 'response':
                response_col = i
                logger.debug(f"Found response column at index {i}")
            elif 'concentration' in col_lower or 'conc' in col_lower or col_lower == 'concentration':
                conc_col = i
                logger.debug(f"Found concentration column at index {i}")
            elif 'cycle' in col_lower or col_lower == 'cycle':
                cycle_col = i
                logger.debug(f"Found cycle column at index {i}")
        
        if time_col is None or response_col is None:
            logger.error(f"Column detection failed. time_col={time_col}, response_col={response_col}, header={header}")
            raise ValueError("Required columns (Time, Response) not found in CSV")
        
        # Parse data rows
        data_rows = csv_data[1:]
        cycles = {}
        
        for row in data_rows:
            if len(row) <= max(time_col, response_col):
                continue
                
            try:
                time_val = float(row[time_col])
                response_val = float(row[response_col])
                
                # Get cycle and concentration
                cycle_val = 1
                if cycle_col is not None and len(row) > cycle_col:
                    try:
                        cycle_val = int(float(row[cycle_col]))
                    except (ValueError, IndexError):
                        pass
                
                conc_val = 0.0
                if conc_col is not None and len(row) > conc_col:
                    try:
                        conc_val = float(row[conc_col])
                    except (ValueError, IndexError):
                        pass
                
                if cycle_val not in cycles:
                    cycles[cycle_val] = {
                        'time': [],
                        'response': [],
                        'concentration': conc_val
                    }
                
                cycles[cycle_val]['time'].append(time_val)
                cycles[cycle_val]['response'].append(response_val)
                
            except (ValueError, IndexError) as e:
                logger.warning(f"Skipping invalid data row: {row[:3]}... ({e})")
                continue
        
        # Convert to sensorgrams
        for cycle_num, cycle_data in cycles.items():
            if len(cycle_data['time']) < 10:  # Skip cycles with too little data
                continue
                
            time_array = np.array(cycle_data['time'])
            response_array = np.array(cycle_data['response'])
            
            # Validate array size
            validate_array_size(time_array)
            validate_array_size(response_array)
            
            # Detect injection phases
            try:
                assoc_start, dissoc_start = detect_injection_phases(time_array, response_array)
            except Exception as e:
                logger.warning(f"Phase detection failed for cycle {cycle_num}: {e}")
                assoc_start = time_array[0]
                dissoc_start = time_array[-1]
            
            sensorgram = Sensorgram(
                cycle=cycle_num,
                concentration=cycle_data['concentration'],
                time=time_array,
                response=response_array,
                association_start=assoc_start,
                dissociation_start=dissoc_start,
                flow_cell=1,
                reference_subtracted=False
            )
            sensorgrams.append(sensorgram)
        
        if not sensorgrams:
            raise ValueError("No valid sensorgrams found in CSV file")
        
        # Sort sensorgrams by cycle number
        sensorgrams.sort(key=lambda x: x.cycle)
        
        spr_data = SPRData(
            instrument=instrument,
            chip_type=chip_type,
            temperature=temperature,
            sensorgrams=sensorgrams,
            metadata=metadata,
            flow_rate=flow_rate,
            buffer=buffer,
            ligand_name=ligand_name,
            analyte_name=analyte_name
        )
        
        logger.info(f"Successfully parsed {len(sensorgrams)} sensorgrams from {path}")
        return spr_data
        
    except Exception as e:
        logger.error(f"Failed to parse Biacore CSV {path}: {e}")
        raise ValueError(f"Invalid Biacore CSV format: {e}") from e


def parse_biacore_sensorgram_txt(path: Union[str, Path]) -> SPRData:
    """
    Parse Biacore sensorgram TXT export file.
    
    Expected format:
    - Space or tab-separated values
    - Header with metadata
    - Two columns: Time (s), Response (RU)
    
    Args:
        path: Path to Biacore TXT file
        
    Returns:
        SPRData object with single sensorgram
        
    Raises:
        FileNotFoundError: If file doesn't exist
        ValueError: If file format is invalid
    """
    path = Path(path)
    if not path.exists():
        raise FileNotFoundError(f"SPR file not found: {path}")
    
    logger.info(f"Parsing Biacore TXT: {path}")
    
    try:
        with open(path, 'r', encoding='utf-8') as f:
            lines = f.readlines()
        
        # Parse metadata from header (first few lines)
        instrument = "Biacore"
        chip_type = "CM5"
        temperature = 25.0
        concentration = 0.0
        metadata = {}
        
        data_start = 0
        for i, line in enumerate(lines):
            line = line.strip()
            if not line or line.startswith('#'):
                continue
            
            # Check if this is a data line (starts with number and has whitespace)
            if re.match(r'^[\d\-\.]+\s+[\d\-\.]', line):
                data_start = i
                break
            
            # Parse metadata
            if ':' in line:
                key, value = line.split(':', 1)
                key = key.strip().lower()
                value = value.strip()
                metadata[key] = value
                
                if 'instrument' in key:
                    instrument = value
                elif 'chip' in key:
                    chip_type = value
                elif 'temperature' in key:
                    try:
                        temperature = float(value.split()[0])
                    except (ValueError, IndexError):
                        pass
                elif 'concentration' in key:
                    try:
                        concentration = float(value.split()[0])
                    except (ValueError, IndexError):
                        pass
        
        # Parse data section
        time_values = []
        response_values = []
        
        for line in lines[data_start:]:
            line = line.strip()
            if not line or line.startswith('#'):
                continue
            
            parts = line.split()
            if len(parts) >= 2:
                try:
                    time_val = float(parts[0])
                    response_val = float(parts[1])
                    time_values.append(time_val)
                    response_values.append(response_val)
                except ValueError:
                    continue
        
        if len(time_values) < 10:
            raise ValueError("Insufficient data points in TXT file")
        
        time_array = np.array(time_values)
        response_array = np.array(response_values)
        
        # Validate array size
        validate_array_size(time_array)
        validate_array_size(response_array)
        
        # Detect injection phases
        try:
            assoc_start, dissoc_start = detect_injection_phases(time_array, response_array)
        except Exception as e:
            logger.warning(f"Phase detection failed: {e}")
            assoc_start = time_array[0]
            dissoc_start = time_array[-1]
        
        sensorgram = Sensorgram(
            cycle=1,
            concentration=concentration,
            time=time_array,
            response=response_array,
            association_start=assoc_start,
            dissociation_start=dissoc_start,
            flow_cell=1,
            reference_subtracted=False
        )
        
        spr_data = SPRData(
            instrument=instrument,
            chip_type=chip_type,
            temperature=temperature,
            sensorgrams=[sensorgram],
            metadata=metadata
        )
        
        logger.info(f"Successfully parsed sensorgram from {path}")
        return spr_data
        
    except Exception as e:
        logger.error(f"Failed to parse Biacore TXT {path}: {e}")
        raise ValueError(f"Invalid Biacore TXT format: {e}") from e


def extract_sensorgrams(data: SPRData) -> List[Sensorgram]:
    """
    Extract sensorgrams from SPRData object.
    
    Args:
        data: SPRData object
        
    Returns:
        List of Sensorgram objects
    """
    return data.sensorgrams


__all__ = [
    "Sensorgram",
    "SPRData",
    "parse_biacore_csv",
    "parse_biacore_sensorgram_txt",
    "extract_sensorgrams",
    "detect_injection_phases",
    "validate_array_size",
]
