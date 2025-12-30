"""
MST (Microscale Thermophoresis) affinity analysis algorithms.

Provides dose-response curve fitting for MST data using Hill equation and
quality assessment metrics for binding affinity determination.
"""

from __future__ import annotations

import logging
import warnings
from dataclasses import dataclass
from typing import List, Tuple

import numpy as np
from lmfit import Model, Parameters

logger = logging.getLogger(__name__)


@dataclass
class AffinityFit:
    """Results from MST dose-response fitting."""
    
    kd: float  # Equilibrium dissociation constant (M)
    kd_error: float  # Standard error of KD
    hill_coefficient: float  # Hill coefficient (cooperativity)
    amplitude: float  # Response amplitude (‰)
    baseline: float  # Response baseline (‰)
    r_squared: float  # Coefficient of determination
    fitted_curve: np.ndarray  # Fitted dose-response curve


@dataclass
class QualityMetrics:
    """Quality assessment metrics for MST data."""
    
    aggregation_detected: bool  # Whether protein aggregation is detected
    photobleaching_percent: float  # Photobleaching percentage
    signal_to_noise: float  # Signal-to-noise ratio
    response_amplitude: float  # Total response amplitude (‰)
    data_quality_score: float  # Overall quality score (0-100)


def _hill_equation(concentration: np.ndarray, kd: float, amplitude: float, 
                   baseline: float, hill_coeff: float = 1.0) -> np.ndarray:
    """
    Hill equation for dose-response fitting.
    
    Fnorm = baseline + amplitude / (1 + (KD/C)^n)
    
    Args:
        concentration: Array of ligand concentrations (M)
        kd: Equilibrium dissociation constant (M)
        amplitude: Response amplitude (‰)
        baseline: Response baseline (‰)
        hill_coeff: Hill coefficient (cooperativity factor)
        
    Returns:
        Array of normalized fluorescence values (‰)
    """
    # Avoid division by zero
    concentration = np.maximum(concentration, 1e-15)
    
    # Hill equation
    binding_fraction = 1.0 / (1.0 + (kd / concentration) ** hill_coeff)
    return baseline + amplitude * binding_fraction


def fit_dose_response(concentrations: np.ndarray, fnorm: np.ndarray, 
                     errors: np.ndarray = None) -> AffinityFit:
    """
    Fit MST dose-response data to Hill equation.
    
    Fits the equation: Fnorm = baseline + amplitude / (1 + (KD/C)^n)
    where n is the Hill coefficient.
    
    Args:
        concentrations: Array of ligand concentrations (M or nM)
        fnorm: Array of normalized fluorescence values (‰)
        errors: Optional array of measurement errors (‰)
        
    Returns:
        AffinityFit object with fitted parameters and quality metrics
        
    Raises:
        ValueError: If arrays have different lengths or insufficient data
    """
    if len(concentrations) != len(fnorm):
        raise ValueError("Concentration and Fnorm arrays must have same length")
    
    if len(concentrations) < 4:
        raise ValueError("At least 4 data points required for dose-response fitting")
    
    # Convert concentrations from nM to M if needed (assume nM if > 1e-6)
    if np.max(concentrations) > 1e-6:
        concentrations = concentrations * 1e-9
        logger.debug("Converted concentrations from nM to M")
    
    # Remove any invalid data points
    valid_mask = np.isfinite(concentrations) & np.isfinite(fnorm) & (concentrations >= 0)
    concentrations = concentrations[valid_mask]
    fnorm = fnorm[valid_mask]
    
    if errors is not None:
        errors = errors[valid_mask]
    
    if len(concentrations) < 4:
        raise ValueError("Insufficient valid data points for fitting")
    
    # Sort by concentration
    sort_idx = np.argsort(concentrations)
    concentrations = concentrations[sort_idx]
    fnorm = fnorm[sort_idx]
    if errors is not None:
        errors = errors[sort_idx]
    
    logger.info(f"Fitting dose-response curve: {len(concentrations)} points")
    
    try:
        # Initial parameter estimates
        baseline_est = fnorm[0] if concentrations[0] == 0 else np.min(fnorm)
        amplitude_est = np.max(fnorm) - baseline_est
        
        # Estimate KD as concentration at half-maximal response
        half_response = baseline_est + amplitude_est / 2
        kd_est = np.interp(half_response, fnorm, concentrations)
        
        # Handle case where interpolation fails
        if not np.isfinite(kd_est) or kd_est <= 0:
            kd_est = np.median(concentrations[concentrations > 0])
        
        logger.debug(f"Initial estimates: KD={kd_est*1e9:.1f} nM, amplitude={amplitude_est:.1f}‰")
        
        # Create lmfit model
        model = Model(_hill_equation)
        params = Parameters()
        
        # Set parameter bounds and initial values
        params.add('kd', value=kd_est, min=concentrations.min()/1000, max=concentrations.max()*1000)
        params.add('amplitude', value=amplitude_est, min=0.1, max=abs(amplitude_est)*10)
        params.add('baseline', value=baseline_est, min=baseline_est-50, max=baseline_est+50)
        params.add('hill_coeff', value=1.0, min=0.1, max=5.0)
        
        # Perform weighted fit if errors provided
        weights = None
        if errors is not None and np.any(errors > 0):
            weights = 1.0 / np.maximum(errors, 0.1)  # Avoid division by zero
        
        result = model.fit(fnorm, params, concentration=concentrations, weights=weights)
        
        if not result.success:
            logger.warning("Dose-response fit failed to converge, using initial estimates")
            kd_fit = kd_est
            kd_error = kd_est * 0.1  # 10% error estimate
            amplitude_fit = amplitude_est
            baseline_fit = baseline_est
            hill_coeff_fit = 1.0
            fitted_curve = _hill_equation(concentrations, kd_fit, amplitude_fit, baseline_fit, hill_coeff_fit)
            r_squared = 0.5  # Poor fit indicator
        else:
            kd_fit = result.params['kd'].value
            kd_error = result.params['kd'].stderr if result.params['kd'].stderr is not None else kd_fit * 0.1
            amplitude_fit = result.params['amplitude'].value
            baseline_fit = result.params['baseline'].value
            hill_coeff_fit = result.params['hill_coeff'].value
            fitted_curve = result.best_fit
            
            # Calculate R-squared
            ss_res = np.sum((fnorm - fitted_curve) ** 2)
            ss_tot = np.sum((fnorm - np.mean(fnorm)) ** 2)
            r_squared = 1 - (ss_res / ss_tot) if ss_tot > 0 else 0
        
        logger.info(f"Fit complete: KD={kd_fit*1e9:.2f}±{kd_error*1e9:.2f} nM, "
                   f"Hill={hill_coeff_fit:.2f}, R²={r_squared:.3f}")
        
        return AffinityFit(
            kd=kd_fit,
            kd_error=kd_error,
            hill_coefficient=hill_coeff_fit,
            amplitude=amplitude_fit,
            baseline=baseline_fit,
            r_squared=r_squared,
            fitted_curve=fitted_curve
        )
        
    except Exception as e:
        logger.error(f"Dose-response fitting failed: {e}")
        raise ValueError(f"Dose-response fitting failed: {e}") from e


def calculate_kd_hill(concs: np.ndarray, response: np.ndarray) -> Tuple[float, float, float]:
    """
    Calculate KD and Hill coefficient from dose-response data.
    
    This is a convenience function that wraps fit_dose_response and returns
    the key parameters as a tuple.
    
    Args:
        concs: Array of concentrations (M or nM)
        response: Array of response values (‰)
        
    Returns:
        Tuple of (KD in M, Hill coefficient, R²)
        
    Raises:
        ValueError: If fitting fails
    """
    fit_result = fit_dose_response(concs, response)
    return fit_result.kd, fit_result.hill_coefficient, fit_result.r_squared


def quality_check(fit: AffinityFit, concentrations: np.ndarray = None, 
                  fnorm: np.ndarray = None) -> QualityMetrics:
    """
    Assess the quality of MST dose-response data and fitting.
    
    Evaluates various quality metrics including aggregation detection,
    signal-to-noise ratio, and overall data quality.
    
    Args:
        fit: AffinityFit object from dose-response fitting
        concentrations: Optional array of concentrations for additional analysis
        fnorm: Optional array of Fnorm values for additional analysis
        
    Returns:
        QualityMetrics object with quality assessment
    """
    logger.debug("Performing MST data quality assessment")
    
    # Initialize quality metrics
    aggregation_detected = False
    photobleaching_percent = 0.0
    signal_to_noise = 0.0
    response_amplitude = abs(fit.amplitude)
    
    # Assess fit quality
    fit_quality_score = min(fit.r_squared * 100, 100)  # Convert R² to 0-100 scale
    
    # Check for aggregation based on response pattern
    if fnorm is not None and concentrations is not None:
        aggregation_detected = detect_aggregation(fnorm)
        
        # Calculate signal-to-noise ratio
        if len(fnorm) > 3:
            # Use first few points (low concentration) to estimate noise
            noise_level = np.std(fnorm[:3]) if len(fnorm) >= 3 else np.std(fnorm)
            signal_to_noise = response_amplitude / max(noise_level, 0.1)
        
        # Estimate photobleaching (simplified)
        # In practice, this would analyze time-course data within each measurement
        photobleaching_percent = min(abs(fit.baseline) / 10, 20)  # Rough estimate
    
    # Calculate overall quality score
    quality_factors = [
        fit_quality_score * 0.4,  # Fit quality (40%)
        min(signal_to_noise * 5, 30),  # S/N ratio (30%, capped at 30)
        min(response_amplitude, 20),  # Response amplitude (20%, capped at 20)
        10 if not aggregation_detected else 0,  # Aggregation penalty (10%)
    ]
    
    data_quality_score = sum(quality_factors)
    data_quality_score = max(0, min(data_quality_score, 100))  # Clamp to 0-100
    
    logger.info(f"Quality assessment: Score={data_quality_score:.1f}, "
               f"S/N={signal_to_noise:.1f}, Aggregation={aggregation_detected}")
    
    return QualityMetrics(
        aggregation_detected=aggregation_detected,
        photobleaching_percent=photobleaching_percent,
        signal_to_noise=signal_to_noise,
        response_amplitude=response_amplitude,
        data_quality_score=data_quality_score
    )


def detect_aggregation(fnorm_values: np.ndarray, threshold: float = 5.0) -> bool:
    """
    Detect protein aggregation in MST data.
    
    Aggregation typically manifests as:
    1. Large positive Fnorm changes at high concentrations
    2. Non-monotonic binding curves
    3. Unusually large response amplitudes
    
    Args:
        fnorm_values: Array of normalized fluorescence values (‰)
        threshold: Threshold for aggregation detection (‰)
        
    Returns:
        True if aggregation is detected
    """
    if len(fnorm_values) < 4:
        return False
    
    # Check for large positive changes (aggregation signature)
    max_response = np.max(fnorm_values)
    min_response = np.min(fnorm_values)
    
    # Aggregation often causes large positive Fnorm changes
    if max_response > threshold and (max_response - min_response) > threshold * 2:
        logger.warning(f"Potential aggregation detected: max Fnorm = {max_response:.1f}‰")
        return True
    
    # Check for non-monotonic behavior (rough heuristic)
    # Calculate differences between consecutive points
    diffs = np.diff(fnorm_values)
    
    # If there are large oscillations, it might indicate aggregation
    if len(diffs) > 3:
        large_changes = np.abs(diffs) > threshold
        if np.sum(large_changes) > len(diffs) * 0.3:  # >30% of points show large changes
            logger.warning("Potential aggregation detected: non-monotonic binding curve")
            return True
    
    # Check for unusually large total response amplitude
    total_amplitude = max_response - min_response
    if total_amplitude > 50:  # Arbitrary threshold for very large responses
        logger.warning(f"Potential aggregation detected: large amplitude = {total_amplitude:.1f}‰")
        return True
    
    return False


# Update module exports
__all__ = [
    "AffinityFit",
    "QualityMetrics",
    "fit_dose_response", 
    "calculate_kd_hill",
    "quality_check",
    "detect_aggregation",
]
