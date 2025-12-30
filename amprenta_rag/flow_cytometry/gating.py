"""
Gating algorithms for flow cytometry analysis.

This module provides functions for applying different types of gates to flow cytometry data:
- Polygon gates using matplotlib path containment
- Rectangle gates with boundary checking
- Quadrant gates returning Q1-Q4 populations
- Boolean gates for combining gate results
- Population statistics computation
"""

from __future__ import annotations

import logging
from dataclasses import dataclass
from typing import Dict, List, Optional
from uuid import UUID

import numpy as np
from matplotlib.path import Path

logger = logging.getLogger(__name__)


@dataclass
class PopulationStats:
    """Statistics for a gated cell population."""
    
    n_events: int
    pct_of_parent: Optional[float] = None
    pct_of_total: Optional[float] = None
    median_values: Optional[Dict[str, float]] = None
    mean_values: Optional[Dict[str, float]] = None
    cv_values: Optional[Dict[str, float]] = None  # Coefficient of variation


def apply_polygon_gate(
    events: np.ndarray,
    x_idx: int,
    y_idx: int,
    vertices: List[List[float]]
) -> np.ndarray:
    """
    Apply polygon gate to flow cytometry events.
    
    Args:
        events: Event data array (n_events, n_parameters)
        x_idx: Column index for X parameter
        y_idx: Column index for Y parameter  
        vertices: List of [x, y] coordinate pairs defining polygon
        
    Returns:
        Boolean mask indicating events inside polygon
        
    Raises:
        ValueError: If vertices has fewer than 3 points or invalid indices
    """
    if len(vertices) < 3:
        raise ValueError("Polygon gate requires at least 3 vertices")
    
    if x_idx >= events.shape[1] or y_idx >= events.shape[1]:
        raise ValueError(f"Parameter indices out of bounds: x_idx={x_idx}, y_idx={y_idx}, n_params={events.shape[1]}")
    
    if x_idx < 0 or y_idx < 0:
        raise ValueError("Parameter indices must be non-negative")
    
    # Extract x, y coordinates
    x_values = events[:, x_idx]
    y_values = events[:, y_idx]
    points = np.column_stack([x_values, y_values])
    
    # Create matplotlib Path and test containment
    polygon_path = Path(vertices)
    mask = polygon_path.contains_points(points)
    
    logger.debug(f"Polygon gate: {np.sum(mask)}/{len(mask)} events inside ({len(vertices)} vertices)")
    
    return mask


def apply_rectangle_gate(
    events: np.ndarray,
    x_idx: int,
    y_idx: int,
    bounds: Dict[str, float]
) -> np.ndarray:
    """
    Apply rectangular gate to flow cytometry events.
    
    Args:
        events: Event data array (n_events, n_parameters)
        x_idx: Column index for X parameter
        y_idx: Column index for Y parameter
        bounds: Dictionary with keys 'x_min', 'x_max', 'y_min', 'y_max'
        
    Returns:
        Boolean mask indicating events inside rectangle
        
    Raises:
        ValueError: If bounds are invalid or indices out of range
    """
    required_keys = {'x_min', 'x_max', 'y_min', 'y_max'}
    if not required_keys.issubset(bounds.keys()):
        raise ValueError(f"Rectangle bounds must contain keys: {required_keys}")
    
    if x_idx >= events.shape[1] or y_idx >= events.shape[1]:
        raise ValueError(f"Parameter indices out of bounds: x_idx={x_idx}, y_idx={y_idx}, n_params={events.shape[1]}")
    
    if x_idx < 0 or y_idx < 0:
        raise ValueError("Parameter indices must be non-negative")
    
    x_min, x_max = bounds['x_min'], bounds['x_max']
    y_min, y_max = bounds['y_min'], bounds['y_max']
    
    if x_min >= x_max or y_min >= y_max:
        raise ValueError("Invalid bounds: min values must be less than max values")
    
    # Extract coordinates
    x_values = events[:, x_idx]
    y_values = events[:, y_idx]
    
    # Apply bounds (inclusive boundaries)
    x_mask = (x_values >= x_min) & (x_values <= x_max)
    y_mask = (y_values >= y_min) & (y_values <= y_max)
    mask = x_mask & y_mask
    
    logger.debug(f"Rectangle gate: {np.sum(mask)}/{len(mask)} events inside bounds")
    
    return mask


def apply_quadrant_gate(
    events: np.ndarray,
    x_idx: int,
    y_idx: int,
    x_thresh: float,
    y_thresh: float
) -> Dict[str, np.ndarray]:
    """
    Apply quadrant gate to flow cytometry events.
    
    Args:
        events: Event data array (n_events, n_parameters)
        x_idx: Column index for X parameter
        y_idx: Column index for Y parameter
        x_thresh: X-axis threshold value
        y_thresh: Y-axis threshold value
        
    Returns:
        Dictionary with keys "Q1", "Q2", "Q3", "Q4" containing boolean masks
        Q1: x < x_thresh, y >= y_thresh (upper left)
        Q2: x >= x_thresh, y >= y_thresh (upper right)  
        Q3: x < x_thresh, y < y_thresh (lower left)
        Q4: x >= x_thresh, y < y_thresh (lower right)
        
    Raises:
        ValueError: If indices are out of bounds
    """
    if x_idx >= events.shape[1] or y_idx >= events.shape[1]:
        raise ValueError(f"Parameter indices out of bounds: x_idx={x_idx}, y_idx={y_idx}, n_params={events.shape[1]}")
    
    if x_idx < 0 or y_idx < 0:
        raise ValueError("Parameter indices must be non-negative")
    
    # Extract coordinates
    x_values = events[:, x_idx]
    y_values = events[:, y_idx]
    
    # Define quadrants
    q1_mask = (x_values < x_thresh) & (y_values >= y_thresh)  # Upper left
    q2_mask = (x_values >= x_thresh) & (y_values >= y_thresh)  # Upper right
    q3_mask = (x_values < x_thresh) & (y_values < y_thresh)  # Lower left
    q4_mask = (x_values >= x_thresh) & (y_values < y_thresh)  # Lower right
    
    quadrants = {
        "Q1": q1_mask,
        "Q2": q2_mask, 
        "Q3": q3_mask,
        "Q4": q4_mask
    }
    
    for q, mask in quadrants.items():
        logger.debug(f"Quadrant {q}: {np.sum(mask)}/{len(mask)} events")
    
    return quadrants


def apply_boolean_gate(
    masks: Dict[UUID, np.ndarray],
    operator: str,
    operand_ids: List[UUID]
) -> np.ndarray:
    """
    Apply boolean logic to combine gate masks.
    
    Args:
        masks: Dictionary mapping gate IDs to boolean masks
        operator: Boolean operator ("AND", "OR", "NOT")
        operand_ids: List of gate IDs to combine
        
    Returns:
        Combined boolean mask
        
    Raises:
        ValueError: If operator is invalid or operand gates not found
    """
    valid_operators = {"AND", "OR", "NOT"}
    if operator not in valid_operators:
        raise ValueError(f"Invalid operator '{operator}'. Must be one of: {valid_operators}")
    
    if not operand_ids:
        raise ValueError("At least one operand gate ID required")
    
    # Validate all operand gates exist
    missing_gates = [gate_id for gate_id in operand_ids if gate_id not in masks]
    if missing_gates:
        raise ValueError(f"Gate IDs not found: {missing_gates}")
    
    operand_masks = [masks[gate_id] for gate_id in operand_ids]
    
    # Ensure all masks have same length
    mask_lengths = [len(mask) for mask in operand_masks]
    if len(set(mask_lengths)) > 1:
        raise ValueError(f"Mask length mismatch: {mask_lengths}")
    
    if operator == "AND":
        # Intersection of all masks
        result = operand_masks[0].copy()
        for mask in operand_masks[1:]:
            result &= mask
    
    elif operator == "OR":
        # Union of all masks
        result = operand_masks[0].copy()
        for mask in operand_masks[1:]:
            result |= mask
    
    elif operator == "NOT":
        # Negation (only use first operand)
        if len(operand_ids) > 1:
            logger.warning(f"NOT operator only uses first operand, ignoring {len(operand_ids)-1} others")
        result = ~operand_masks[0]
    
    logger.debug(f"Boolean gate {operator}: {np.sum(result)}/{len(result)} events")
    
    return result


def compute_population_stats(
    events: np.ndarray,
    mask: np.ndarray,
    param_names: List[str],
    parent_mask: Optional[np.ndarray] = None,
    total_events: Optional[int] = None
) -> PopulationStats:
    """
    Compute statistics for a gated cell population.
    
    Args:
        events: Event data array (n_events, n_parameters)
        mask: Boolean mask indicating population events
        param_names: List of parameter names
        parent_mask: Boolean mask for parent population (for percentage calculation)
        total_events: Total event count (for percentage calculation)
        
    Returns:
        PopulationStats object with computed statistics
        
    Raises:
        ValueError: If mask length doesn't match events or param names mismatch
    """
    if len(mask) != events.shape[0]:
        raise ValueError(f"Mask length {len(mask)} doesn't match event count {events.shape[0]}")
    
    if len(param_names) != events.shape[1]:
        raise ValueError(f"Parameter name count {len(param_names)} doesn't match data columns {events.shape[1]}")
    
    n_events = int(np.sum(mask))
    
    # Calculate percentages
    pct_of_parent = None
    pct_of_total = None
    
    if parent_mask is not None:
        parent_count = int(np.sum(parent_mask))
        pct_of_parent = (n_events / parent_count * 100.0) if parent_count > 0 else 0.0
    
    if total_events is not None:
        pct_of_total = (n_events / total_events * 100.0) if total_events > 0 else 0.0
    
    # Compute parameter statistics for gated events
    median_values = {}
    mean_values = {}
    cv_values = {}
    
    if n_events > 0:
        gated_events = events[mask]
        
        for i, param_name in enumerate(param_names):
            param_data = gated_events[:, i]
            
            median_values[param_name] = float(np.median(param_data))
            mean_values[param_name] = float(np.mean(param_data))
            
            # Coefficient of variation (CV = std / mean)
            std_val = float(np.std(param_data))
            mean_val = mean_values[param_name]
            cv_values[param_name] = (std_val / abs(mean_val) * 100.0) if mean_val != 0 else 0.0
    
    else:
        # Empty population - set all stats to zero/None
        for param_name in param_names:
            median_values[param_name] = 0.0
            mean_values[param_name] = 0.0
            cv_values[param_name] = 0.0
    
    parent_str = f"{pct_of_parent:.1f}%" if pct_of_parent is not None else "N/A"
    total_str = f"{pct_of_total:.1f}%" if pct_of_total is not None else "N/A"
    logger.debug(f"Population stats: {n_events} events, {parent_str} of parent, {total_str} of total")
    
    return PopulationStats(
        n_events=n_events,
        pct_of_parent=pct_of_parent,
        pct_of_total=pct_of_total,
        median_values=median_values,
        mean_values=mean_values,
        cv_values=cv_values
    )


__all__ = [
    "PopulationStats",
    "apply_polygon_gate",
    "apply_rectangle_gate", 
    "apply_quadrant_gate",
    "apply_boolean_gate",
    "compute_population_stats",
]
