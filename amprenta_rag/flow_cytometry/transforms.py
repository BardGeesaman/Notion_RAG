"""
Mathematical transformations for flow cytometry data analysis.

This module provides standard transformations used in flow cytometry:
- Logicle transformation for handling negative values and wide dynamic range
- Arcsinh transformation for CyTOF data
- Compensation matrix application for spectral overlap correction
- Event subsampling for performance and visualization

All transformations follow established flow cytometry standards and best practices.
"""

from __future__ import annotations

from typing import Optional, Tuple

import numpy as np

from amprenta_rag.logging_utils import get_logger

logger = get_logger(__name__)


def logicle_transform(
    data: np.ndarray, 
    T: float, 
    M: float, 
    W: float, 
    A: float = 0.0
) -> np.ndarray:
    """
    Apply logicle (biexponential) transformation to flow cytometry data.
    
    The logicle transformation provides a smooth transition between linear 
    and logarithmic scales, handling both positive and negative values.
    
    Args:
        data: Input data array
        T: Top of scale value (maximum data value)
        M: Number of decades in log portion
        W: Width of linear portion in decades
        A: Additional decades of negative data
        
    Returns:
        Transformed data array
        
    Raises:
        ValueError: If parameters violate logicle bounds
        
    References:
        Parks, Roederer, Moore (2006) Cytometry A 69:541-551
    """
    # Validate parameters
    if T <= 0:
        raise ValueError(f"T must be positive, got {T}")
    if M <= 0:
        raise ValueError(f"M must be positive, got {M}")
    if W <= 0 or W >= M/2:
        raise ValueError(f"W must be in range (0, M/2), got W={W}, M={M}")
    if A < 0:
        raise ValueError(f"A must be non-negative, got {A}")
    
    # Calculate derived parameters
    w = W / (M + A)
    x2 = A / (M + A)
    x1 = x2 + w
    x0 = x2 + 2*w
    b = (M + A) * np.log(10)
    
    # Calculate constants
    d = _solve_logicle_d(b, w)
    c_a = np.exp(x0 * (b + d))
    m_f_a = np.exp(b * x1) - c_a / np.exp(d * x1)
    
    # Normalize input data
    x = data / T
    
    # Apply transformation
    result = np.zeros_like(x)
    
    # Linear region: x <= x1
    linear_mask = x <= x1
    if np.any(linear_mask):
        result[linear_mask] = (
            (np.exp(b * x[linear_mask]) - c_a / np.exp(d * x[linear_mask])) / m_f_a
        ) * (M + A)
    
    # Log region: x > x1  
    log_mask = x > x1
    if np.any(log_mask):
        result[log_mask] = (
            (np.exp(b * x1) - c_a / np.exp(d * x1)) / m_f_a +
            b * (x[log_mask] - x1) / m_f_a
        ) * (M + A)
    
    return result


def _solve_logicle_d(b: float, w: float) -> float:
    """
    Solve for parameter d in logicle transformation using iterative method.
    
    Args:
        b: Derived parameter b
        w: Derived parameter w
        
    Returns:
        Parameter d value
    """
    # Initial guess - ensure positive value
    d = max(0.1, 2 * (np.log(max(w, 1e-10)) - np.log(2)) / b)
    
    # Newton-Raphson iteration
    for _ in range(10):  # Usually converges quickly
        if d <= 0:
            d = 0.1  # Reset if d becomes non-positive
            
        f = 2 * np.log(d) - b * d + 2 * np.log(max(w, 1e-10)) - 2 * np.log(2)
        df = 2/d - b
        
        if abs(df) < 1e-10:
            break
            
        d_new = d - f/df
        
        # Ensure d stays positive
        if d_new <= 0:
            d_new = d / 2
            
        if abs(d_new - d) < 1e-10:
            break
            
        d = d_new
    
    return max(d, 0.1)  # Ensure positive return value


def arcsinh_transform(data: np.ndarray, cofactor: float = 150.0) -> np.ndarray:
    """
    Apply arcsinh (inverse hyperbolic sine) transformation.
    
    Commonly used for CyTOF (mass cytometry) data to handle the wide
    dynamic range and reduce the influence of outliers.
    
    Args:
        data: Input data array
        cofactor: Scaling factor (default 150 for CyTOF)
        
    Returns:
        Transformed data array
        
    Raises:
        ValueError: If cofactor is not positive
    """
    if cofactor <= 0:
        raise ValueError(f"Cofactor must be positive, got {cofactor}")
    
    return np.arcsinh(data / cofactor)


def apply_compensation(data: np.ndarray, spillover: np.ndarray) -> np.ndarray:
    """
    Apply spectral compensation using spillover matrix.
    
    Corrects for spectral overlap between fluorescent channels by
    applying the inverse of the spillover matrix.
    
    Args:
        data: Event data (n_events x n_channels)
        spillover: Spillover matrix (n_channels x n_channels)
        
    Returns:
        Compensated data array
        
    Raises:
        ValueError: If matrix dimensions don't match
        np.linalg.LinAlgError: If spillover matrix is singular
    """
    if data.shape[1] != spillover.shape[0]:
        raise ValueError(
            f"Data channels ({data.shape[1]}) must match spillover matrix size ({spillover.shape[0]})"
        )
    
    if spillover.shape[0] != spillover.shape[1]:
        raise ValueError(
            f"Spillover matrix must be square, got shape {spillover.shape}"
        )
    
    try:
        # Compute inverse spillover matrix
        inv_spillover = np.linalg.inv(spillover)
        
        # Apply compensation: compensated = data * inv_spillover^T
        compensated = data @ inv_spillover.T
        
        logger.info(f"Applied compensation to {data.shape[0]} events")
        return compensated
        
    except np.linalg.LinAlgError as e:
        raise np.linalg.LinAlgError(f"Cannot invert spillover matrix: {e}") from e


def auto_logicle_params(data: np.ndarray, percentile: float = 0.05) -> Tuple[float, float, float, float]:
    """
    Automatically estimate logicle transformation parameters from data.
    
    Args:
        data: Input data array (single channel or full dataset)
        percentile: Percentile for negative data estimation (default 0.05)
        
    Returns:
        Tuple of (T, M, W, A) parameters
        
    Raises:
        ValueError: If data is empty or percentile invalid
    """
    if data.size == 0:
        raise ValueError("Data array cannot be empty")
    
    if not 0 < percentile < 1:
        raise ValueError(f"Percentile must be in (0, 1), got {percentile}")
    
    # Flatten data if multidimensional
    flat_data = data.flatten()
    
    # Estimate T (top of scale) from positive data
    positive_data = flat_data[flat_data > 0]
    if len(positive_data) > 0:
        T = float(np.percentile(positive_data, 99.9))
    else:
        T = 1000.0  # Default fallback
    
    # Ensure T is reasonable
    T = max(T, 100.0)
    
    # Standard values for flow cytometry
    M = 4.5  # ~4.5 decades is typical
    
    # Estimate W from data spread around zero
    near_zero = flat_data[np.abs(flat_data) < T * 0.1]
    if len(near_zero) > 10:
        spread = np.std(near_zero)
        W = max(0.1, min(1.0, spread / T * M))
    else:
        W = 0.5  # Default
    
    # Estimate A from negative data
    negative_data = flat_data[flat_data < 0]
    if len(negative_data) > 0:
        neg_extent = abs(float(np.percentile(negative_data, percentile * 100)))
        A = max(0.0, min(2.0, np.log10(neg_extent / T) + M))
    else:
        A = 0.0
    
    # Validate bounds
    W = min(W, M/2 - 0.01)  # Ensure W < M/2
    A = max(A, 0.0)         # Ensure A >= 0
    
    logger.info(f"Auto-estimated logicle params: T={T:.1f}, M={M:.1f}, W={W:.2f}, A={A:.2f}")
    
    return T, M, W, A


def subsample_events(
    events: np.ndarray, 
    max_events: int = 50000, 
    random_state: Optional[int] = None
) -> np.ndarray:
    """
    Randomly subsample events for performance and visualization.
    
    Args:
        events: Event data array (n_events x n_parameters)
        max_events: Maximum number of events to keep
        random_state: Random seed for reproducibility
        
    Returns:
        Subsampled event array
        
    Raises:
        ValueError: If max_events is not positive
    """
    if max_events <= 0:
        raise ValueError(f"max_events must be positive, got {max_events}")
    
    n_events = events.shape[0]
    
    if n_events <= max_events:
        logger.info(f"No subsampling needed: {n_events} <= {max_events}")
        return events
    
    # Set random seed if provided
    if random_state is not None:
        np.random.seed(random_state)
    
    # Random sampling without replacement
    indices = np.random.choice(n_events, size=max_events, replace=False)
    indices.sort()  # Keep temporal order
    
    subsampled = events[indices]
    
    logger.info(f"Subsampled {n_events} events to {max_events}")
    
    return subsampled
