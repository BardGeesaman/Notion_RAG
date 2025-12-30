"""
DSC (Differential Scanning Calorimetry) data parsers for MicroCal and TA Instruments.

Supports parsing of DSC CSV and TXT files with thermogram data and automatic
baseline correction algorithms.
"""

from __future__ import annotations

import csv
import logging
import re
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Dict, List, Optional, Union

import numpy as np
from scipy import signal

logger = logging.getLogger(__name__)

# Maximum data points per thermogram for performance  
MAX_ARRAY_SIZE = 50000


@dataclass
class Scan:
    """Single DSC thermogram scan data."""
    
    scan_number: int
    temperature: np.ndarray  # Temperature points (°C)
    heat_capacity: np.ndarray  # Heat capacity (kcal/mol/°C)
    baseline_subtracted: bool = False  # Whether baseline was subtracted
    scan_type: str = "sample"  # "sample", "reference", "buffer"
    scan_rate: float = 1.0  # Heating rate (°C/min)
    

@dataclass
class DSCData:
    """Complete DSC dataset from instrument file."""
    
    instrument: str  # Instrument model (e.g., "MicroCal PEAQ-DSC")
    scan_rate: float  # Heating rate (°C/min)
    protein_concentration: float  # Protein concentration (mg/mL or μM)
    scans: List[Scan]  # List of thermogram scans
    metadata: Dict[str, Any]  # Additional metadata from file
    buffer: Optional[str] = None  # Buffer composition
    ligand_concentration: float = 0.0  # Ligand concentration for binding studies
    cell_volume: Optional[float] = None  # Cell volume (μL)
    reference_buffer: Optional[str] = None  # Reference cell buffer


def baseline_correction(temp: np.ndarray, cp: np.ndarray, 
                       method: str = "linear") -> np.ndarray:
    """
    Apply baseline correction to DSC thermogram.
    
    Args:
        temp: Temperature array (°C)
        cp: Heat capacity array (kcal/mol/°C)
        method: Correction method ("linear", "polynomial", "spline")
        
    Returns:
        Baseline-corrected heat capacity array
        
    Raises:
        ValueError: If arrays have different lengths or invalid method
    """
    if len(temp) != len(cp):
        raise ValueError("Temperature and heat capacity arrays must have same length")
    
    if len(temp) < 10:
        raise ValueError("Insufficient data points for baseline correction")
    
    if method not in ["linear", "polynomial", "spline"]:
        raise ValueError(f"Invalid baseline correction method: {method}")
    
    logger.debug(f"Applying {method} baseline correction to {len(temp)} points")
    
    # Find baseline regions (typically first and last 10% of data)
    n_points = len(temp)
    baseline_fraction = 0.1
    n_baseline = max(int(n_points * baseline_fraction), 5)
    
    # Get baseline indices
    start_indices = np.arange(n_baseline)
    end_indices = np.arange(n_points - n_baseline, n_points)
    baseline_indices = np.concatenate([start_indices, end_indices])
    
    baseline_temp = temp[baseline_indices]
    baseline_cp = cp[baseline_indices]
    
    if method == "linear":
        # Linear baseline fit
        coeffs = np.polyfit(baseline_temp, baseline_cp, 1)
        baseline = np.polyval(coeffs, temp)
        
    elif method == "polynomial":
        # Polynomial baseline fit (degree 2)
        coeffs = np.polyfit(baseline_temp, baseline_cp, 2)
        baseline = np.polyval(coeffs, temp)
        
    elif method == "spline":
        # Spline baseline fit
        try:
            from scipy.interpolate import UnivariateSpline
            spline = UnivariateSpline(baseline_temp, baseline_cp, s=len(baseline_temp))
            baseline = spline(temp)
        except ImportError:
            logger.warning("scipy not available, falling back to linear baseline")
            coeffs = np.polyfit(baseline_temp, baseline_cp, 1)
            baseline = np.polyval(coeffs, temp)
    
    corrected_cp = cp - baseline
    
    logger.debug(f"Baseline correction complete. Mean correction: {np.mean(baseline):.3f}")
    
    return corrected_cp


def parse_microcal_csv(path: Union[str, Path]) -> DSCData:
    """
    Parse MicroCal DSC CSV export file.
    
    Expected format:
    - Comma-separated values
    - Header rows with metadata (instrument, scan rate, concentration)
    - Data columns: Temperature (°C), Cp (kcal/mol/°C)
    
    Args:
        path: Path to MicroCal CSV file
        
    Returns:
        DSCData object with parsed thermogram
        
    Raises:
        FileNotFoundError: If file doesn't exist
        ValueError: If file format is invalid
    """
    path = Path(path)
    if not path.exists():
        raise FileNotFoundError(f"DSC file not found: {path}")
    
    logger.info(f"Parsing MicroCal CSV: {path}")
    
    try:
        with open(path, 'r', encoding='utf-8') as f:
            lines = f.readlines()
        
        # Parse metadata from header
        instrument = "MicroCal PEAQ-DSC"
        scan_rate = 1.0
        protein_concentration = 0.0
        buffer = None
        ligand_concentration = 0.0
        cell_volume = None
        reference_buffer = None
        metadata = {}
        
        data_start = 0
        for i, line in enumerate(lines):
            line = line.strip()
            if not line:
                continue
            
            # Check if this looks like a header line (contains column names)
            if not data_start and ('temperature' in line.lower() or 'cp' in line.lower()):
                data_start = i
                break
            # Check if this looks like a data line (starts with number and has commas)
            elif re.match(r'^[\d\-\.]', line) and ',' in line:
                if data_start == 0:
                    # No header found, assume first data line and create default header
                    parts = line.split(',')
                    if len(parts) >= 2:
                        lines.insert(i, "Temperature,Cp")  # Default header
                        data_start = i
                        break
                else:
                    data_start = i
                    break
            
            # Parse metadata
            line_lower = line.lower()
            if ':' in line:
                key, value = line.split(':', 1)
                key = key.strip().lower()
                value = value.strip()
                metadata[key] = value
                
                if 'instrument' in key or 'model' in key:
                    instrument = value
                elif 'scan rate' in key or 'heating rate' in key:
                    try:
                        scan_rate = float(re.findall(r'[\d\.]+', value)[0])
                    except (IndexError, ValueError):
                        pass
                elif 'concentration' in key and 'protein' in key:
                    try:
                        protein_concentration = float(re.findall(r'[\d\.]+', value)[0])
                    except (IndexError, ValueError):
                        pass
                elif 'concentration' in key and ('ligand' in key or 'compound' in key):
                    try:
                        ligand_concentration = float(re.findall(r'[\d\.]+', value)[0])
                    except (IndexError, ValueError):
                        pass
                elif 'buffer' in key:
                    buffer = value
                elif 'cell volume' in key or 'volume' in key:
                    try:
                        cell_volume = float(re.findall(r'[\d\.]+', value)[0])
                    except (IndexError, ValueError):
                        pass
                elif 'reference' in key and 'buffer' in key:
                    reference_buffer = value
            
            # Also check for inline metadata
            elif 'instrument' in line_lower:
                match = re.search(r'instrument[:\s]+([^,\n\r]+)', line, re.IGNORECASE)
                if match:
                    instrument = match.group(1).strip()
            elif 'scan rate' in line_lower or 'heating rate' in line_lower:
                match = re.search(r'([\d\.]+)', line)
                if match:
                    scan_rate = float(match.group(1))
        
        # Parse CSV data
        csv_reader = csv.reader(lines[data_start:])
        rows = list(csv_reader)
        
        if len(rows) < 2:
            raise ValueError("Insufficient data in CSV file")
        
        # Determine column indices
        header = [col.strip().lower() for col in rows[0]]
        temp_col = None
        cp_col = None
        
        for i, col in enumerate(header):
            if 'temperature' in col or 'temp' in col or col == 'temperature':
                temp_col = i
            elif 'cp' in col or 'heat capacity' in col or 'capacity' in col or col == 'cp':
                cp_col = i
        
        if temp_col is None or cp_col is None:
            raise ValueError("Required columns (Temperature, Cp) not found in CSV")
        
        # Parse data rows
        temperatures = []
        heat_capacities = []
        
        for row in rows[1:]:
            if len(row) <= max(temp_col, cp_col):
                continue
            
            try:
                temp_val = float(row[temp_col])
                cp_val = float(row[cp_col])
                temperatures.append(temp_val)
                heat_capacities.append(cp_val)
            except (ValueError, IndexError) as e:
                logger.warning(f"Skipping invalid data row: {row[:2]}... ({e})")
                continue
        
        if len(temperatures) < 10:
            raise ValueError("Insufficient data points in CSV file")
        
        temp_array = np.array(temperatures)
        cp_array = np.array(heat_capacities)
        
        # Validate array size
        if len(temp_array) > MAX_ARRAY_SIZE:
            raise ValueError(
                f"Array size {len(temp_array)} exceeds maximum {MAX_ARRAY_SIZE} points. "
                f"Consider downsampling the data."
            )
        
        scan = Scan(
            scan_number=1,
            temperature=temp_array,
            heat_capacity=cp_array,
            baseline_subtracted=False,
            scan_type="sample",
            scan_rate=scan_rate
        )
        
        dsc_data = DSCData(
            instrument=instrument,
            scan_rate=scan_rate,
            protein_concentration=protein_concentration,
            scans=[scan],
            metadata=metadata,
            buffer=buffer,
            ligand_concentration=ligand_concentration,
            cell_volume=cell_volume,
            reference_buffer=reference_buffer
        )
        
        logger.info(f"Successfully parsed thermogram from {path}")
        return dsc_data
        
    except Exception as e:
        logger.error(f"Failed to parse MicroCal CSV {path}: {e}")
        raise ValueError(f"Invalid MicroCal CSV format: {e}") from e


def parse_ta_instruments_txt(path: Union[str, Path]) -> DSCData:
    """
    Parse TA Instruments DSC TXT export file.
    
    Expected format:
    - Tab or space-separated values
    - Header with metadata
    - Data columns: Temperature (°C), Heat Flow (mW), Cp (J/g/°C)
    
    Args:
        path: Path to TA Instruments TXT file
        
    Returns:
        DSCData object with parsed thermogram
        
    Raises:
        FileNotFoundError: If file doesn't exist
        ValueError: If file format is invalid
    """
    path = Path(path)
    if not path.exists():
        raise FileNotFoundError(f"DSC file not found: {path}")
    
    logger.info(f"Parsing TA Instruments TXT: {path}")
    
    try:
        with open(path, 'r', encoding='utf-8') as f:
            lines = f.readlines()
        
        # Parse metadata from header
        instrument = "TA Instruments DSC"
        scan_rate = 10.0  # TA typically uses higher scan rates
        protein_concentration = 0.0
        buffer = None
        ligand_concentration = 0.0
        cell_volume = None
        reference_buffer = None
        metadata = {}
        
        data_start = 0
        for i, line in enumerate(lines):
            line = line.strip()
            if not line or line.startswith('#'):
                continue
            
            # Check if this is a data line (starts with number)
            if re.match(r'^[\d\-\.]', line):
                data_start = i
                break
            
            # Parse metadata
            if ':' in line:
                key, value = line.split(':', 1)
                key = key.strip().lower()
                value = value.strip()
                metadata[key] = value
                
                if 'instrument' in key or 'model' in key:
                    instrument = value
                elif 'scan rate' in key or 'heating rate' in key or 'rate' in key:
                    try:
                        scan_rate = float(re.findall(r'[\d\.]+', value)[0])
                    except (IndexError, ValueError):
                        pass
                elif 'concentration' in key and 'protein' in key:
                    try:
                        protein_concentration = float(re.findall(r'[\d\.]+', value)[0])
                    except (IndexError, ValueError):
                        pass
                elif 'concentration' in key and ('ligand' in key or 'compound' in key):
                    try:
                        ligand_concentration = float(re.findall(r'[\d\.]+', value)[0])
                    except (IndexError, ValueError):
                        pass
                elif 'buffer' in key:
                    buffer = value
                elif 'cell volume' in key or 'volume' in key:
                    try:
                        cell_volume = float(re.findall(r'[\d\.]+', value)[0])
                    except (IndexError, ValueError):
                        pass
        
        # Parse data section
        temperatures = []
        heat_capacities = []
        
        # Determine column format from first data line
        first_data_line = lines[data_start].strip()
        delimiter = '\t' if '\t' in first_data_line else None
        
        for line in lines[data_start:]:
            line = line.strip()
            if not line or line.startswith('#'):
                continue
            
            if delimiter:
                parts = line.split(delimiter)
            else:
                parts = line.split()
            
            if len(parts) >= 2:
                try:
                    temp_val = float(parts[0])  # First column: Temperature
                    
                    # Try different column positions for heat capacity
                    cp_val = None
                    for col_idx in [1, 2]:  # Try columns 1 and 2
                        if len(parts) > col_idx:
                            try:
                                cp_candidate = float(parts[col_idx])
                                # TA Instruments often has heat flow in mW, convert to appropriate units
                                # Assume if values are very large (>1000), they're in different units
                                if abs(cp_candidate) > 1000:
                                    cp_candidate = cp_candidate / 1000.0  # Convert mW to W
                                cp_val = cp_candidate
                                break
                            except ValueError:
                                continue
                    
                    if cp_val is not None:
                        temperatures.append(temp_val)
                        heat_capacities.append(cp_val)
                        
                except ValueError:
                    continue
        
        if len(temperatures) < 10:
            raise ValueError("Insufficient data points in TXT file")
        
        temp_array = np.array(temperatures)
        cp_array = np.array(heat_capacities)
        
        # Validate array size
        if len(temp_array) > MAX_ARRAY_SIZE:
            raise ValueError(
                f"Array size {len(temp_array)} exceeds maximum {MAX_ARRAY_SIZE} points. "
                f"Consider downsampling the data."
            )
        
        scan = Scan(
            scan_number=1,
            temperature=temp_array,
            heat_capacity=cp_array,
            baseline_subtracted=False,
            scan_type="sample",
            scan_rate=scan_rate
        )
        
        dsc_data = DSCData(
            instrument=instrument,
            scan_rate=scan_rate,
            protein_concentration=protein_concentration,
            scans=[scan],
            metadata=metadata,
            buffer=buffer,
            ligand_concentration=ligand_concentration,
            cell_volume=cell_volume,
            reference_buffer=reference_buffer
        )
        
        logger.info(f"Successfully parsed thermogram from {path}")
        return dsc_data
        
    except Exception as e:
        logger.error(f"Failed to parse TA Instruments TXT {path}: {e}")
        raise ValueError(f"Invalid TA Instruments TXT format: {e}") from e


__all__ = [
    "Scan",
    "DSCData",
    "parse_microcal_csv",
    "parse_ta_instruments_txt",
    "baseline_correction",
]
