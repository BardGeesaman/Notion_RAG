"""
FCS file parsing and data I/O for flow cytometry analysis.

This module provides functions for reading FCS files, extracting metadata,
validating file format, and converting event data to/from Parquet format
for efficient storage and retrieval.

Supports FCS 2.0, 3.0, and 3.1 file formats.
"""

from __future__ import annotations

import os
import warnings
from dataclasses import dataclass
from datetime import datetime
from pathlib import Path
from typing import Dict, List, Optional, Tuple, Union

import numpy as np
import pandas as pd
import pyarrow as pa
import pyarrow.parquet as pq
from fcsparser import parse

from amprenta_rag.logging_utils import get_logger

logger = get_logger(__name__)


@dataclass
class FlowMetadata:
    """Structured metadata extracted from FCS files."""
    
    # File information
    filename: str
    file_size_bytes: int
    fcs_version: str
    
    # Acquisition information
    acquisition_date: Optional[datetime] = None
    cytometer_model: Optional[str] = None
    cytometer_serial: Optional[str] = None
    acquisition_software: Optional[str] = None
    
    # Sample information
    sample_id: Optional[str] = None
    sample_volume: Optional[float] = None
    
    # Data characteristics
    n_events: int = 0
    n_parameters: int = 0
    parameter_names: List[str] = None
    parameter_ranges: Dict[str, Tuple[float, float]] = None
    
    # Compensation and transformation
    spillover_matrix: Optional[np.ndarray] = None
    compensation_applied: bool = False
    
    # Additional metadata
    keywords: Dict[str, str] = None
    
    def __post_init__(self):
        """Initialize default values for mutable fields."""
        if self.parameter_names is None:
            self.parameter_names = []
        if self.parameter_ranges is None:
            self.parameter_ranges = {}
        if self.keywords is None:
            self.keywords = {}


def load_fcs(path: Union[str, Path]) -> Tuple[Dict, np.ndarray]:
    """
    Load FCS file and return metadata and event data.
    
    Args:
        path: Path to FCS file
        
    Returns:
        Tuple of (metadata_dict, events_array)
        
    Raises:
        FileNotFoundError: If file doesn't exist
        ValueError: If file is not a valid FCS format
        RuntimeError: If parsing fails
    """
    path = Path(path)
    
    if not path.exists():
        raise FileNotFoundError(f"FCS file not found: {path}")
    
    if not path.suffix.lower() in ['.fcs', '.lmd']:
        logger.warning(f"File extension {path.suffix} is not typical for FCS files")
    
    try:
        logger.info(f"Loading FCS file: {path}")
        
        # Parse FCS file using fcsparser
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")  # Suppress fcsparser warnings
            meta, data = parse(str(path), reformat_meta=True)
        
        # Convert data to numpy array
        if isinstance(data, pd.DataFrame):
            events = data.values.astype(np.float32)
        else:
            events = np.array(data, dtype=np.float32)
        
        logger.info(f"Loaded {events.shape[0]} events with {events.shape[1]} parameters")
        
        return meta, events
        
    except Exception as e:
        raise RuntimeError(f"Failed to parse FCS file {path}: {e}") from e


def extract_metadata(meta: Dict, filename: str = "", file_size: int = 0) -> FlowMetadata:
    """
    Extract structured metadata from FCS metadata dictionary.
    
    Args:
        meta: Metadata dictionary from fcsparser
        filename: Original filename
        file_size: File size in bytes
        
    Returns:
        FlowMetadata object with structured information
    """
    # Extract basic file info
    fcs_version = meta.get('$FIL', meta.get('FIL', 'Unknown'))
    n_events = int(meta.get('$TOT', meta.get('TOT', 0)))
    n_parameters = int(meta.get('$PAR', meta.get('PAR', 0)))
    
    # Extract parameter information
    parameter_names = []
    parameter_ranges = {}
    
    for i in range(1, n_parameters + 1):
        # Parameter name (short name preferred)
        param_name = (
            meta.get(f'$P{i}S', meta.get(f'P{i}S', '')) or
            meta.get(f'$P{i}N', meta.get(f'P{i}N', f'P{i}'))
        )
        parameter_names.append(param_name)
        
        # Parameter range
        param_range = float(meta.get(f'$P{i}R', meta.get(f'P{i}R', 0)))
        if param_range > 0:
            parameter_ranges[param_name] = (0.0, param_range)
    
    # Extract acquisition date
    acquisition_date = None
    date_str = meta.get('$DATE', meta.get('DATE', ''))
    time_str = meta.get('$BTIM', meta.get('BTIM', ''))
    
    if date_str:
        try:
            if time_str:
                datetime_str = f"{date_str} {time_str}"
                acquisition_date = datetime.strptime(datetime_str, '%d-%b-%Y %H:%M:%S')
            else:
                acquisition_date = datetime.strptime(date_str, '%d-%b-%Y')
        except ValueError:
            logger.warning(f"Could not parse acquisition date: {date_str} {time_str}")
    
    # Extract instrument information
    cytometer_model = meta.get('$CYT', meta.get('CYT', None))
    cytometer_serial = meta.get('$CYTSN', meta.get('CYTSN', None))
    acquisition_software = meta.get('$SRC', meta.get('SRC', None))
    
    # Extract sample information
    sample_id = meta.get('$SMNO', meta.get('SMNO', None))
    sample_volume_str = meta.get('$VOL', meta.get('VOL', ''))
    sample_volume = None
    if sample_volume_str:
        try:
            sample_volume = float(sample_volume_str)
        except ValueError:
            pass
    
    # Extract spillover/compensation matrix
    spillover_matrix = None
    compensation_applied = False
    
    # Look for spillover matrix in various formats
    spill_key = None
    for key in ['$SPILLOVER', 'SPILLOVER', '$SPILL', 'SPILL']:
        if key in meta:
            spill_key = key
            break
    
    if spill_key and meta[spill_key]:
        try:
            spill_data = meta[spill_key]
            if isinstance(spill_data, str):
                # Parse comma-separated spillover matrix
                parts = spill_data.split(',')
                n_markers = int(parts[0])
                
                # Extract numeric values only (skip parameter names)
                numeric_values = []
                for part in parts[1:]:  # Skip the first element (n_markers)
                    try:
                        numeric_values.append(float(part))
                    except ValueError:
                        # Skip non-numeric parts (parameter names)
                        continue
                
                # Need exactly n_markers^2 values for the matrix
                if len(numeric_values) >= n_markers * n_markers:
                    matrix_values = numeric_values[:n_markers * n_markers]
                    spillover_matrix = np.array(matrix_values).reshape(n_markers, n_markers)
            elif isinstance(spill_data, (list, np.ndarray)):
                spillover_matrix = np.array(spill_data)
        except (ValueError, IndexError) as e:
            logger.warning(f"Could not parse spillover matrix: {e}")
    
    # Check if compensation was already applied
    comp_applied = meta.get('$COMP', meta.get('COMP', ''))
    if comp_applied and comp_applied.upper() in ['TRUE', '1', 'YES']:
        compensation_applied = True
    
    return FlowMetadata(
        filename=filename,
        file_size_bytes=file_size,
        fcs_version=fcs_version,
        acquisition_date=acquisition_date,
        cytometer_model=cytometer_model,
        cytometer_serial=cytometer_serial,
        acquisition_software=acquisition_software,
        sample_id=sample_id,
        sample_volume=sample_volume,
        n_events=n_events,
        n_parameters=n_parameters,
        parameter_names=parameter_names,
        parameter_ranges=parameter_ranges,
        spillover_matrix=spillover_matrix,
        compensation_applied=compensation_applied,
        keywords=dict(meta),  # Store all keywords
    )


def validate_fcs(path: Union[str, Path]) -> List[str]:
    """
    Validate FCS file format and return list of issues found.
    
    Args:
        path: Path to FCS file
        
    Returns:
        List of validation error/warning messages (empty if valid)
    """
    issues = []
    path = Path(path)
    
    # Check file existence
    if not path.exists():
        issues.append(f"File not found: {path}")
        return issues
    
    # Check file size
    if path.stat().st_size == 0:
        issues.append("File is empty")
        return issues
    
    try:
        # Try to load the file
        meta, events = load_fcs(path)
        
        # Validate metadata
        n_events = int(meta.get('$TOT', meta.get('TOT', 0)))
        n_parameters = int(meta.get('$PAR', meta.get('PAR', 0)))
        
        if n_events <= 0:
            issues.append(f"Invalid event count: {n_events}")
        
        if n_parameters <= 0:
            issues.append(f"Invalid parameter count: {n_parameters}")
        
        # Validate event data shape
        if events.shape[0] != n_events:
            issues.append(f"Event count mismatch: header={n_events}, data={events.shape[0]}")
        
        if events.shape[1] != n_parameters:
            issues.append(f"Parameter count mismatch: header={n_parameters}, data={events.shape[1]}")
        
        # Check for missing parameter names
        missing_names = 0
        for i in range(1, n_parameters + 1):
            param_name = meta.get(f'$P{i}N', meta.get(f'P{i}N', ''))
            if not param_name:
                missing_names += 1
        
        if missing_names > 0:
            issues.append(f"{missing_names} parameters missing names")
        
        # Check for infinite or NaN values
        if np.any(np.isinf(events)):
            issues.append("Data contains infinite values")
        
        if np.any(np.isnan(events)):
            issues.append("Data contains NaN values")
        
        # Validate FCS version
        fcs_version = meta.get('$FIL', meta.get('FIL', ''))
        if not fcs_version.startswith('FCS'):
            issues.append(f"Unrecognized FCS version: {fcs_version}")
        
        logger.info(f"FCS validation complete: {len(issues)} issues found")
        
    except Exception as e:
        issues.append(f"Failed to parse file: {e}")
    
    return issues


def save_events_parquet(
    events: np.ndarray, 
    param_names: List[str], 
    path: Union[str, Path],
    metadata: Optional[Dict] = None
) -> None:
    """
    Save event data to Parquet format for efficient storage.
    
    Args:
        events: Event data array (n_events x n_parameters)
        param_names: List of parameter names
        path: Output Parquet file path
        metadata: Optional metadata to store with file
        
    Raises:
        ValueError: If events and param_names dimensions don't match
        OSError: If file cannot be written
    """
    if events.shape[1] != len(param_names):
        raise ValueError(
            f"Parameter count mismatch: events={events.shape[1]}, names={len(param_names)}"
        )
    
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    
    try:
        # Create DataFrame with proper column names
        df = pd.DataFrame(events, columns=param_names)
        
        # Convert to PyArrow table for efficient storage
        table = pa.Table.from_pandas(df, preserve_index=False)
        
        # Add metadata if provided
        if metadata:
            # Convert metadata to string format for Parquet
            meta_dict = {str(k): str(v) for k, v in metadata.items()}
            existing_meta = table.schema.metadata or {}
            existing_meta.update(meta_dict)
            table = table.replace_schema_metadata(existing_meta)
        
        # Write to Parquet with compression
        pq.write_table(
            table, 
            str(path),
            compression='snappy',
            use_dictionary=True,
            row_group_size=50000
        )
        
        logger.info(f"Saved {events.shape[0]} events to {path}")
        
    except Exception as e:
        raise OSError(f"Failed to save Parquet file {path}: {e}") from e


def load_events_parquet(path: Union[str, Path]) -> pd.DataFrame:
    """
    Load event data from Parquet format.
    
    Args:
        path: Path to Parquet file
        
    Returns:
        DataFrame with event data
        
    Raises:
        FileNotFoundError: If file doesn't exist
        OSError: If file cannot be read
    """
    path = Path(path)
    
    if not path.exists():
        raise FileNotFoundError(f"Parquet file not found: {path}")
    
    try:
        # Read Parquet file
        table = pq.read_table(str(path))
        df = table.to_pandas()
        
        logger.info(f"Loaded {len(df)} events from {path}")
        return df
        
    except Exception as e:
        raise OSError(f"Failed to load Parquet file {path}: {e}") from e
