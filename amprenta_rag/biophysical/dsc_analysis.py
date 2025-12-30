"""
DSC (Differential Scanning Calorimetry) thermal analysis algorithms.

Provides thermal unfolding analysis, peak detection, and transition deconvolution
for DSC thermogram data.
"""

from __future__ import annotations

import logging
import warnings
from dataclasses import dataclass
from typing import List, Tuple

import numpy as np
from lmfit import Model, Parameters
from scipy import signal

logger = logging.getLogger(__name__)


@dataclass
class ThermalFit:
    """Results from DSC thermal transition fitting."""
    
    tm: float  # Melting temperature (°C)
    tm_error: float  # Standard error of Tm
    delta_h: float  # Enthalpy change (kcal/mol)
    delta_cp: float  # Heat capacity change (kcal/mol/°C)
    onset_temp: float  # Transition onset temperature (°C)
    cooperativity: float  # van't Hoff ratio (cooperativity factor)
    fitted_curve: np.ndarray  # Fitted thermogram curve


@dataclass
class Peak:
    """Individual peak detected in DSC thermogram."""
    
    temperature: float  # Peak temperature (°C)
    height: float  # Peak height (kcal/mol/°C)
    width: float  # Peak width at half maximum (°C)
    area: float  # Peak area (kcal/mol)


def _two_state_model(temperature: np.ndarray, tm: float, delta_h: float, 
                     delta_cp: float, cp_baseline: float) -> np.ndarray:
    """
    Two-state thermal unfolding model.
    
    Cp = Cp_baseline + delta_Cp * fU + delta_H * dfU/dT
    
    where fU is the fraction unfolded:
    fU = 1 / (1 + exp(delta_H/R * (1/T - 1/Tm)))
    
    Args:
        temperature: Array of temperatures (°C)
        tm: Melting temperature (°C)
        delta_h: Enthalpy change (kcal/mol)
        delta_cp: Heat capacity change (kcal/mol/°C)
        cp_baseline: Baseline heat capacity (kcal/mol/°C)
        
    Returns:
        Array of heat capacity values (kcal/mol/°C)
    """
    # Convert temperature to Kelvin
    T_K = temperature + 273.15
    Tm_K = tm + 273.15
    
    # Gas constant in kcal/mol/K
    R = 1.987e-3  # kcal/mol/K
    
    # Calculate fraction unfolded
    # Avoid overflow by limiting the exponent
    exponent = (delta_h / R) * (1/T_K - 1/Tm_K)
    exponent = np.clip(exponent, -500, 500)  # Prevent overflow
    
    fU = 1.0 / (1.0 + np.exp(-exponent))
    
    # Calculate derivative dfU/dT for the transition peak
    dfU_dT = (delta_h / (R * T_K**2)) * fU * (1 - fU)
    
    # Two-state heat capacity
    cp = cp_baseline + delta_cp * fU + delta_h * dfU_dT
    
    return cp


def _gaussian_peak(temperature: np.ndarray, center: float, height: float, 
                   width: float, baseline: float = 0.0) -> np.ndarray:
    """
    Gaussian peak model for DSC peak fitting.
    
    Args:
        temperature: Array of temperatures (°C)
        center: Peak center temperature (°C)
        height: Peak height (kcal/mol/°C)
        width: Peak width parameter (°C)
        baseline: Baseline offset (kcal/mol/°C)
        
    Returns:
        Array of heat capacity values (kcal/mol/°C)
    """
    return baseline + height * np.exp(-0.5 * ((temperature - center) / width)**2)


def fit_two_state_unfolding(temp: np.ndarray, cp: np.ndarray) -> ThermalFit:
    """
    Fit DSC thermogram to two-state unfolding model.
    
    Fits the thermodynamic model for protein unfolding:
    N ⇌ U (native ⇌ unfolded)
    
    Args:
        temp: Array of temperatures (°C)
        cp: Array of heat capacity values (kcal/mol/°C)
        
    Returns:
        ThermalFit object with fitted thermodynamic parameters
        
    Raises:
        ValueError: If arrays have different lengths or insufficient data
    """
    if len(temp) != len(cp):
        raise ValueError("Temperature and heat capacity arrays must have same length")
    
    if len(temp) < 20:
        raise ValueError("At least 20 data points required for thermal fitting")
    
    # Remove any invalid data points
    valid_mask = np.isfinite(temp) & np.isfinite(cp)
    temp = temp[valid_mask]
    cp = cp[valid_mask]
    
    if len(temp) < 20:
        raise ValueError("Insufficient valid data points for fitting")
    
    # Sort by temperature
    sort_idx = np.argsort(temp)
    temp = temp[sort_idx]
    cp = cp[sort_idx]
    
    logger.info(f"Fitting two-state unfolding model: {len(temp)} points, "
               f"T range {temp[0]:.1f}-{temp[-1]:.1f}°C")
    
    try:
        # Initial parameter estimates
        cp_baseline = np.mean(cp[:5])  # Baseline from first few points
        cp_max = np.max(cp)
        tm_est = temp[np.argmax(cp)]  # Temperature at maximum Cp
        
        # Estimate delta_H from peak area (rough approximation)
        delta_h_est = 50.0  # Typical value for protein unfolding (kcal/mol)
        
        # Estimate delta_Cp from difference between folded and unfolded baselines
        if len(cp) > 10:
            cp_final = np.mean(cp[-5:])  # Final baseline
            delta_cp_est = cp_final - cp_baseline
        else:
            delta_cp_est = 0.0
        
        logger.debug(f"Initial estimates: Tm={tm_est:.1f}°C, ΔH={delta_h_est:.1f} kcal/mol")
        
        # Create lmfit model
        model = Model(_two_state_model)
        params = Parameters()
        
        # Set parameter bounds and initial values
        params.add('tm', value=tm_est, min=temp.min(), max=temp.max())
        params.add('delta_h', value=delta_h_est, min=5.0, max=200.0)
        params.add('delta_cp', value=delta_cp_est, min=-5.0, max=5.0)
        params.add('cp_baseline', value=cp_baseline, min=cp_baseline-1.0, max=cp_baseline+1.0)
        
        # Perform fit
        result = model.fit(cp, params, temperature=temp)
        
        if not result.success:
            logger.warning("Two-state fit failed to converge, using peak detection fallback")
            # Fallback to simple peak detection
            peaks = detect_peaks(temp, cp)
            if peaks:
                main_peak = max(peaks, key=lambda p: p.height)
                tm_fit = main_peak.temperature
                tm_error = 1.0
                delta_h_fit = main_peak.area
                delta_cp_fit = delta_cp_est
                fitted_curve = _gaussian_peak(temp, tm_fit, main_peak.height, 
                                            main_peak.width, cp_baseline)
            else:
                tm_fit = tm_est
                tm_error = 5.0
                delta_h_fit = delta_h_est
                delta_cp_fit = delta_cp_est
                fitted_curve = np.full_like(cp, cp_baseline)
        else:
            tm_fit = result.params['tm'].value
            tm_error = result.params['tm'].stderr if result.params['tm'].stderr is not None else 1.0
            delta_h_fit = result.params['delta_h'].value
            delta_cp_fit = result.params['delta_cp'].value
            fitted_curve = result.best_fit
        
        # Calculate onset temperature (temperature where transition begins)
        # Defined as temperature where Cp exceeds baseline by 10%
        cp_threshold = cp_baseline + 0.1 * (np.max(fitted_curve) - cp_baseline)
        onset_indices = np.where(fitted_curve > cp_threshold)[0]
        onset_temp = temp[onset_indices[0]] if len(onset_indices) > 0 else tm_fit - 10
        
        # Calculate cooperativity (van't Hoff ratio)
        # Simplified calculation based on peak width
        peak_height = np.max(fitted_curve) - cp_baseline
        if peak_height > 0:
            # Find half-maximum points
            half_max = cp_baseline + peak_height / 2
            half_max_indices = np.where(np.abs(fitted_curve - half_max) < peak_height * 0.1)[0]
            if len(half_max_indices) >= 2:
                width_at_half_max = temp[half_max_indices[-1]] - temp[half_max_indices[0]]
                # Cooperativity inversely related to peak width
                cooperativity = 50.0 / max(width_at_half_max, 1.0)  # Empirical scaling
            else:
                cooperativity = 1.0
        else:
            cooperativity = 1.0
        
        cooperativity = max(0.1, min(cooperativity, 10.0))  # Clamp to reasonable range
        
        logger.info(f"Fit complete: Tm={tm_fit:.2f}±{tm_error:.2f}°C, "
                   f"ΔH={delta_h_fit:.1f} kcal/mol, cooperativity={cooperativity:.2f}")
        
        return ThermalFit(
            tm=tm_fit,
            tm_error=tm_error,
            delta_h=delta_h_fit,
            delta_cp=delta_cp_fit,
            onset_temp=onset_temp,
            cooperativity=cooperativity,
            fitted_curve=fitted_curve
        )
        
    except Exception as e:
        logger.error(f"Thermal fitting failed: {e}")
        raise ValueError(f"Thermal fitting failed: {e}") from e


def detect_peaks(temp: np.ndarray, cp: np.ndarray, min_height: float = None, 
                 min_distance: float = 5.0) -> List[Peak]:
    """
    Detect peaks in DSC thermogram using scipy peak detection.
    
    Args:
        temp: Array of temperatures (°C)
        cp: Array of heat capacity values (kcal/mol/°C)
        min_height: Minimum peak height (auto-detected if None)
        min_distance: Minimum distance between peaks (°C)
        
    Returns:
        List of Peak objects for detected peaks
        
    Raises:
        ValueError: If arrays have different lengths
    """
    if len(temp) != len(cp):
        raise ValueError("Temperature and heat capacity arrays must have same length")
    
    if len(temp) < 10:
        return []
    
    # Sort by temperature
    sort_idx = np.argsort(temp)
    temp = temp[sort_idx]
    cp = cp[sort_idx]
    
    logger.debug(f"Detecting peaks in thermogram: {len(temp)} points")
    
    try:
        # Auto-detect minimum height if not provided
        if min_height is None:
            baseline = np.percentile(cp, 10)  # Use 10th percentile as baseline
            noise_level = np.std(cp[:10]) if len(cp) > 10 else np.std(cp) * 0.1
            min_height = baseline + 3 * noise_level  # 3-sigma above baseline
        
        # Convert minimum distance from temperature to index units
        temp_spacing = np.mean(np.diff(temp))
        min_distance_idx = max(1, int(min_distance / temp_spacing))
        
        # Find peaks using scipy
        peak_indices, peak_properties = signal.find_peaks(
            cp, 
            height=min_height,
            distance=min_distance_idx,
            width=1  # Minimum width of 1 data point
        )
        
        peaks = []
        
        for i, peak_idx in enumerate(peak_indices):
            peak_temp = temp[peak_idx]
            peak_height = cp[peak_idx]
            
            # Calculate peak width at half maximum
            if 'widths' in peak_properties:
                width_samples = peak_properties['widths'][i]
                width_temp = width_samples * temp_spacing
            else:
                # Fallback width calculation
                half_height = peak_height / 2
                left_idx = peak_idx
                right_idx = peak_idx
                
                # Find left half-maximum
                while left_idx > 0 and cp[left_idx] > half_height:
                    left_idx -= 1
                
                # Find right half-maximum
                while right_idx < len(cp) - 1 and cp[right_idx] > half_height:
                    right_idx += 1
                
                width_temp = temp[right_idx] - temp[left_idx]
            
            # Calculate peak area (approximate integration)
            # Use a window around the peak
            window_size = max(5, int(width_temp / temp_spacing))
            left_bound = max(0, peak_idx - window_size)
            right_bound = min(len(cp), peak_idx + window_size + 1)
            
            # Simple trapezoidal integration
            temp_window = temp[left_bound:right_bound]
            cp_window = cp[left_bound:right_bound] - np.min(cp)  # Baseline subtraction
            peak_area = np.trapz(cp_window, temp_window)
            
            peak = Peak(
                temperature=peak_temp,
                height=peak_height,
                width=width_temp,
                area=peak_area
            )
            peaks.append(peak)
        
        logger.info(f"Detected {len(peaks)} peaks")
        for i, peak in enumerate(peaks):
            logger.debug(f"Peak {i+1}: T={peak.temperature:.1f}°C, "
                        f"height={peak.height:.3f}, width={peak.width:.1f}°C")
        
        return peaks
        
    except Exception as e:
        logger.error(f"Peak detection failed: {e}")
        return []


def calculate_reversibility(scan1, scan2, temperature_tolerance: float = 2.0) -> float:
    """
    Calculate thermal reversibility by comparing two DSC scans.
    
    Reversibility is assessed by comparing peak positions and heights
    between heating and cooling scans, or between first and second heating scans.
    
    Args:
        scan1: First Scan object (e.g., first heating)
        scan2: Second Scan object (e.g., cooling or second heating)
        temperature_tolerance: Temperature tolerance for peak matching (°C)
        
    Returns:
        Reversibility percentage (0-100%)
        
    Raises:
        ValueError: If scans have insufficient data
    """
    if len(scan1.temperature) < 10 or len(scan2.temperature) < 10:
        raise ValueError("Insufficient data in scans for reversibility calculation")
    
    logger.debug("Calculating thermal reversibility between scans")
    
    try:
        # Detect peaks in both scans
        peaks1 = detect_peaks(scan1.temperature, scan1.heat_capacity)
        peaks2 = detect_peaks(scan2.temperature, scan2.heat_capacity)
        
        if not peaks1 or not peaks2:
            logger.warning("No peaks detected in one or both scans")
            return 0.0
        
        # Match peaks between scans based on temperature
        matched_peaks = []
        
        for peak1 in peaks1:
            best_match = None
            min_temp_diff = float('inf')
            
            for peak2 in peaks2:
                temp_diff = abs(peak1.temperature - peak2.temperature)
                if temp_diff < temperature_tolerance and temp_diff < min_temp_diff:
                    min_temp_diff = temp_diff
                    best_match = peak2
            
            if best_match is not None:
                matched_peaks.append((peak1, best_match))
        
        if not matched_peaks:
            logger.warning("No matching peaks found between scans")
            return 0.0
        
        # Calculate reversibility based on peak height conservation
        reversibilities = []
        
        for peak1, peak2 in matched_peaks:
            # Height reversibility (how well peak height is preserved)
            max_height = max(peak1.height, peak2.height)
            min_height = min(peak1.height, peak2.height)
            height_reversibility = (min_height / max_height) * 100 if max_height > 0 else 0
            
            # Area reversibility (how well peak area is preserved)
            max_area = max(peak1.area, peak2.area)
            min_area = min(peak1.area, peak2.area)
            area_reversibility = (min_area / max_area) * 100 if max_area > 0 else 0
            
            # Combined reversibility (average of height and area)
            peak_reversibility = (height_reversibility + area_reversibility) / 2
            reversibilities.append(peak_reversibility)
        
        # Overall reversibility is the average of all matched peaks
        overall_reversibility = np.mean(reversibilities)
        
        logger.info(f"Reversibility analysis: {len(matched_peaks)} matched peaks, "
                   f"{overall_reversibility:.1f}% reversible")
        
        return overall_reversibility
        
    except Exception as e:
        logger.error(f"Reversibility calculation failed: {e}")
        return 0.0


def deconvolute_transitions(temp: np.ndarray, cp: np.ndarray, n_peaks: int) -> List[ThermalFit]:
    """
    Deconvolute multiple thermal transitions in DSC data.
    
    Fits multiple overlapping thermal transitions using a combination of
    Gaussian peaks or two-state models.
    
    Args:
        temp: Array of temperatures (°C)
        cp: Array of heat capacity values (kcal/mol/°C)
        n_peaks: Number of transitions to fit
        
    Returns:
        List of ThermalFit objects for each deconvoluted transition
        
    Raises:
        ValueError: If arrays have different lengths or invalid n_peaks
    """
    if len(temp) != len(cp):
        raise ValueError("Temperature and heat capacity arrays must have same length")
    
    if n_peaks < 1 or n_peaks > 5:
        raise ValueError("Number of peaks must be between 1 and 5")
    
    if len(temp) < n_peaks * 10:
        raise ValueError(f"Insufficient data points for {n_peaks} peak deconvolution")
    
    logger.info(f"Deconvoluting {n_peaks} thermal transitions")
    
    try:
        # First, detect peaks to get initial estimates
        detected_peaks = detect_peaks(temp, cp)
        
        if len(detected_peaks) < n_peaks:
            logger.warning(f"Only {len(detected_peaks)} peaks detected, requested {n_peaks}")
            # Pad with additional peaks if needed
            while len(detected_peaks) < n_peaks:
                # Add peaks at evenly spaced temperatures
                temp_range = temp[-1] - temp[0]
                new_temp = temp[0] + (len(detected_peaks) + 1) * temp_range / (n_peaks + 1)
                new_peak = Peak(
                    temperature=new_temp,
                    height=np.mean(cp),
                    width=temp_range / (n_peaks * 2),
                    area=np.mean(cp) * temp_range / n_peaks
                )
                detected_peaks.append(new_peak)
        
        # Take the n_peaks highest peaks
        detected_peaks.sort(key=lambda p: p.height, reverse=True)
        selected_peaks = detected_peaks[:n_peaks]
        selected_peaks.sort(key=lambda p: p.temperature)  # Sort by temperature
        
        # For simplicity, fit each peak region individually
        # In a full implementation, this would fit all peaks simultaneously
        thermal_fits = []
        
        for i, peak in enumerate(selected_peaks):
            # Define region around this peak
            peak_temp = peak.temperature
            peak_width = peak.width
            
            # Create a window around the peak (±2 widths)
            temp_window = 2 * peak_width
            temp_mask = (temp >= peak_temp - temp_window) & (temp <= peak_temp + temp_window)
            
            if np.sum(temp_mask) < 10:
                # Expand window if too few points
                temp_mask = (temp >= peak_temp - 10) & (temp <= peak_temp + 10)
            
            if np.sum(temp_mask) < 5:
                logger.warning(f"Insufficient data for peak {i+1} at {peak_temp:.1f}°C")
                continue
            
            temp_region = temp[temp_mask]
            cp_region = cp[temp_mask]
            
            try:
                # Fit this region as a single transition
                fit = fit_two_state_unfolding(temp_region, cp_region)
                thermal_fits.append(fit)
            except Exception as e:
                logger.warning(f"Failed to fit peak {i+1}: {e}")
                # Create a fallback fit based on the detected peak
                fallback_fit = ThermalFit(
                    tm=peak.temperature,
                    tm_error=2.0,
                    delta_h=peak.area,
                    delta_cp=0.0,
                    onset_temp=peak.temperature - peak.width,
                    cooperativity=1.0,
                    fitted_curve=_gaussian_peak(temp_region, peak.temperature, 
                                              peak.height, peak.width, np.min(cp_region))
                )
                thermal_fits.append(fallback_fit)
        
        logger.info(f"Deconvolution complete: {len(thermal_fits)} transitions fitted")
        
        return thermal_fits
        
    except Exception as e:
        logger.error(f"Transition deconvolution failed: {e}")
        raise ValueError(f"Transition deconvolution failed: {e}") from e


# Update module exports
__all__ = [
    "ThermalFit",
    "Peak",
    "fit_two_state_unfolding",
    "detect_peaks",
    "calculate_reversibility", 
    "deconvolute_transitions",
]
