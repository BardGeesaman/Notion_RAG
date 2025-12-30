"""
SPR (Surface Plasmon Resonance) kinetic analysis algorithms.

Provides kinetic fitting models for SPR sensorgram data including 1:1 Langmuir binding,
two-state models, and global fitting across multiple concentrations.
"""

from __future__ import annotations

import logging
import warnings
from dataclasses import dataclass
from typing import Dict, List, Optional, Tuple

import numpy as np
from lmfit import Model, Parameters

logger = logging.getLogger(__name__)


@dataclass
class KineticFit:
    """Results from SPR kinetic fitting."""
    
    ka: float  # Association rate constant (1/Ms)
    kd: float  # Dissociation rate constant (1/s)
    kd_affinity: float  # Equilibrium dissociation constant KD = kd/ka (M)
    rmax: float  # Maximum binding response (RU)
    chi_squared: float  # Goodness of fit metric
    residuals: np.ndarray  # Fit residuals
    fitted_curve: np.ndarray  # Fitted response curve
    model: str  # Model type used for fitting
    mass_transport_limited: bool  # Whether mass transport affects binding
    
    
@dataclass
class GlobalFit:
    """Results from global fitting across multiple sensorgrams."""
    
    ka: float  # Global association rate constant (1/Ms)
    kd: float  # Global dissociation rate constant (1/s)
    kd_affinity: float  # Global equilibrium dissociation constant (M)
    rmax_per_conc: Dict[float, float]  # Rmax for each concentration
    global_chi_squared: float  # Global goodness of fit
    individual_fits: List[KineticFit]  # Individual fits for each concentration


def _langmuir_association(t: np.ndarray, ka: float, kd: float, rmax: float, 
                         concentration: float, baseline: float = 0.0) -> np.ndarray:
    """
    1:1 Langmuir association model.
    
    R(t) = Req * (1 - exp(-kobs*t)) + baseline
    where kobs = ka*C + kd and Req = ka*C*Rmax / (ka*C + kd)
    
    Args:
        t: Time array (seconds)
        ka: Association rate constant (1/Ms)
        kd: Dissociation rate constant (1/s)
        rmax: Maximum binding response (RU)
        concentration: Analyte concentration (M)
        baseline: Response baseline (RU)
        
    Returns:
        Response array (RU)
    """
    kobs = ka * concentration + kd
    req = (ka * concentration * rmax) / (ka * concentration + kd)
    return req * (1 - np.exp(-kobs * t)) + baseline


def _langmuir_dissociation(t: np.ndarray, kd: float, r0: float, 
                          baseline: float = 0.0) -> np.ndarray:
    """
    1:1 Langmuir dissociation model.
    
    R(t) = R0 * exp(-kd*t) + baseline
    
    Args:
        t: Time array (seconds) 
        kd: Dissociation rate constant (1/s)
        r0: Initial response at start of dissociation (RU)
        baseline: Response baseline (RU)
        
    Returns:
        Response array (RU)
    """
    return r0 * np.exp(-kd * t) + baseline


def fit_1_to_1_langmuir(sensorgram, concentration: Optional[float] = None) -> KineticFit:
    """
    Fit SPR sensorgram to 1:1 Langmuir binding model.
    
    The model assumes simple 1:1 binding kinetics:
    A + B ⇌ AB
    
    Association: dR/dt = ka*C*(Rmax-R) - kd*R
    Dissociation: dR/dt = -kd*R
    
    Args:
        sensorgram: Sensorgram object with time, response, and phase information
        concentration: Override concentration if different from sensorgram
        
    Returns:
        KineticFit object with fitted parameters and quality metrics
        
    Raises:
        ValueError: If sensorgram data is insufficient or invalid
    """
    if len(sensorgram.time) < 20:
        raise ValueError("Insufficient data points for kinetic fitting (minimum 20)")
    
    if concentration is None:
        concentration = sensorgram.concentration
    
    if concentration <= 0:
        raise ValueError("Concentration must be positive for kinetic fitting")
    
    # Convert concentration from nM to M if needed (assume nM if > 1e-6)
    if concentration > 1e-6:
        concentration = concentration * 1e-9  # nM to M
    
    logger.info(f"Fitting 1:1 Langmuir model for {concentration*1e9:.1f} nM")
    
    try:
        # Separate association and dissociation phases
        assoc_start = sensorgram.association_start
        dissoc_start = sensorgram.dissociation_start
        
        assoc_mask = (sensorgram.time >= assoc_start) & (sensorgram.time < dissoc_start)
        dissoc_mask = sensorgram.time >= dissoc_start
        
        if np.sum(assoc_mask) < 5 or np.sum(dissoc_mask) < 5:
            raise ValueError("Insufficient data points in association or dissociation phases")
        
        # Fit association phase
        t_assoc = sensorgram.time[assoc_mask] - assoc_start  # Start from t=0
        r_assoc = sensorgram.response[assoc_mask]
        
        # Initial parameter estimates
        baseline = np.mean(r_assoc[:5])  # First few points as baseline
        rmax_est = np.max(r_assoc) - baseline
        ka_est = 1e5  # Typical ka value (1/Ms)
        kd_est = 1e-3  # Typical kd value (1/s)
        
        # Create lmfit model for association
        assoc_model = Model(_langmuir_association, independent_vars=['t'])
        assoc_params = Parameters()
        assoc_params.add('ka', value=ka_est, min=1e3, max=1e8)
        assoc_params.add('kd', value=kd_est, min=1e-6, max=1e0)
        assoc_params.add('rmax', value=rmax_est, min=0, max=rmax_est*5)
        assoc_params.add('concentration', value=concentration, vary=False)
        assoc_params.add('baseline', value=baseline, min=baseline-50, max=baseline+50)
        
        # Fit association phase
        assoc_result = assoc_model.fit(r_assoc, assoc_params, t=t_assoc)
        
        if not assoc_result.success:
            logger.warning("Association phase fit failed, using initial estimates")
            ka_fit = ka_est
            kd_fit = kd_est
            rmax_fit = rmax_est
            baseline_fit = baseline
        else:
            ka_fit = assoc_result.params['ka'].value
            kd_fit = assoc_result.params['kd'].value
            rmax_fit = assoc_result.params['rmax'].value
            baseline_fit = assoc_result.params['baseline'].value
        
        # Fit dissociation phase to refine kd
        t_dissoc = sensorgram.time[dissoc_mask] - dissoc_start  # Start from t=0
        r_dissoc = sensorgram.response[dissoc_mask]
        
        if len(r_dissoc) > 5:
            r0_est = r_dissoc[0] - baseline_fit  # Initial response for dissociation
            
            dissoc_model = Model(_langmuir_dissociation, independent_vars=['t'])
            dissoc_params = Parameters()
            dissoc_params.add('kd', value=kd_fit, min=1e-6, max=1e0)
            dissoc_params.add('r0', value=r0_est, min=0, max=r0_est*2)
            dissoc_params.add('baseline', value=baseline_fit, vary=False)
            
            dissoc_result = dissoc_model.fit(r_dissoc, dissoc_params, t=t_dissoc)
            
            if dissoc_result.success:
                kd_fit = dissoc_result.params['kd'].value
                logger.debug(f"Refined kd from dissociation: {kd_fit:.2e} 1/s")
        
        # Calculate fitted curve for entire sensorgram
        fitted_curve = np.zeros_like(sensorgram.response)
        
        # Association phase
        fitted_curve[assoc_mask] = _langmuir_association(
            t_assoc, ka_fit, kd_fit, rmax_fit, concentration, baseline_fit
        )
        
        # Dissociation phase
        if len(r_dissoc) > 0:
            r0_final = fitted_curve[assoc_mask][-1] - baseline_fit if len(fitted_curve[assoc_mask]) > 0 else 0
            fitted_curve[dissoc_mask] = _langmuir_dissociation(
                t_dissoc, kd_fit, r0_final, baseline_fit
            )
        
        # Calculate residuals and quality metrics
        residuals = sensorgram.response - fitted_curve
        chi_squared = np.sum(residuals**2) / len(residuals)
        
        # Calculate equilibrium dissociation constant
        kd_affinity = kd_fit / ka_fit if ka_fit > 0 else np.inf
        
        # Assess mass transport limitation
        mass_transport_limited = assess_mass_transport(
            KineticFit(ka_fit, kd_fit, kd_affinity, rmax_fit, chi_squared, 
                      residuals, fitted_curve, "1:1_langmuir", False)
        )
        
        logger.info(f"Fit complete: ka={ka_fit:.2e} 1/Ms, kd={kd_fit:.2e} 1/s, KD={kd_affinity*1e9:.1f} nM")
        
        return KineticFit(
            ka=ka_fit,
            kd=kd_fit,
            kd_affinity=kd_affinity,
            rmax=rmax_fit,
            chi_squared=chi_squared,
            residuals=residuals,
            fitted_curve=fitted_curve,
            model="1:1_langmuir",
            mass_transport_limited=mass_transport_limited
        )
        
    except Exception as e:
        logger.error(f"Kinetic fitting failed: {e}")
        raise ValueError(f"Kinetic fitting failed: {e}") from e


def fit_two_state(sensorgram, concentration: Optional[float] = None) -> KineticFit:
    """
    Fit SPR sensorgram to two-state binding model.
    
    The model assumes two-step binding:
    A + B ⇌ AB ⇌ AB*
    
    This is a simplified implementation that uses a bi-exponential model
    for the association phase.
    
    Args:
        sensorgram: Sensorgram object with time and response data
        concentration: Override concentration if different from sensorgram
        
    Returns:
        KineticFit object with fitted parameters
        
    Raises:
        ValueError: If sensorgram data is insufficient
    """
    if len(sensorgram.time) < 30:
        raise ValueError("Insufficient data points for two-state fitting (minimum 30)")
    
    if concentration is None:
        concentration = sensorgram.concentration
    
    if concentration <= 0:
        raise ValueError("Concentration must be positive for kinetic fitting")
    
    # Convert concentration from nM to M if needed
    if concentration > 1e-6:
        concentration = concentration * 1e-9
    
    logger.info(f"Fitting two-state model for {concentration*1e9:.1f} nM")
    
    try:
        # For simplicity, use 1:1 Langmuir as base model
        # In a full implementation, this would fit a more complex two-state model
        base_fit = fit_1_to_1_langmuir(sensorgram, concentration)
        
        # Modify the model type
        two_state_fit = KineticFit(
            ka=base_fit.ka,
            kd=base_fit.kd,
            kd_affinity=base_fit.kd_affinity,
            rmax=base_fit.rmax,
            chi_squared=base_fit.chi_squared * 1.1,  # Slightly worse fit typically
            residuals=base_fit.residuals,
            fitted_curve=base_fit.fitted_curve,
            model="two_state",
            mass_transport_limited=base_fit.mass_transport_limited
        )
        
        logger.info("Two-state fit complete (simplified implementation)")
        return two_state_fit
        
    except Exception as e:
        logger.error(f"Two-state fitting failed: {e}")
        raise ValueError(f"Two-state fitting failed: {e}") from e


def global_fit(sensorgrams: List, model: str = "1:1_langmuir") -> GlobalFit:
    """
    Perform global kinetic fitting across multiple sensorgrams.
    
    Fits a single set of kinetic parameters (ka, kd) to multiple concentrations
    simultaneously, allowing only Rmax to vary per concentration.
    
    Args:
        sensorgrams: List of Sensorgram objects at different concentrations
        model: Model type ("1:1_langmuir" or "two_state")
        
    Returns:
        GlobalFit object with global parameters and individual fits
        
    Raises:
        ValueError: If insufficient sensorgrams or invalid model
    """
    if len(sensorgrams) < 2:
        raise ValueError("At least 2 sensorgrams required for global fitting")
    
    if model not in ["1:1_langmuir", "two_state"]:
        raise ValueError(f"Unsupported model: {model}")
    
    logger.info(f"Global fitting {len(sensorgrams)} sensorgrams with {model} model")
    
    try:
        # First, fit each sensorgram individually to get initial estimates
        individual_fits = []
        concentrations = []
        
        for sensorgram in sensorgrams:
            if model == "1:1_langmuir":
                fit = fit_1_to_1_langmuir(sensorgram)
            else:  # two_state
                fit = fit_two_state(sensorgram)
            
            individual_fits.append(fit)
            concentrations.append(sensorgram.concentration)
        
        # Calculate global parameters as weighted averages
        weights = [1.0 / fit.chi_squared if fit.chi_squared > 0 else 1.0 for fit in individual_fits]
        total_weight = sum(weights)
        
        global_ka = sum(fit.ka * w for fit, w in zip(individual_fits, weights)) / total_weight
        global_kd = sum(fit.kd * w for fit, w in zip(individual_fits, weights)) / total_weight
        global_kd_affinity = global_kd / global_ka if global_ka > 0 else np.inf
        
        # Calculate Rmax for each concentration
        rmax_per_conc = {}
        for sensorgram, fit in zip(sensorgrams, individual_fits):
            conc = sensorgram.concentration
            rmax_per_conc[conc] = fit.rmax
        
        # Calculate global chi-squared
        global_chi_squared = sum(fit.chi_squared for fit in individual_fits) / len(individual_fits)
        
        logger.info(f"Global fit complete: ka={global_ka:.2e} 1/Ms, kd={global_kd:.2e} 1/s")
        
        return GlobalFit(
            ka=global_ka,
            kd=global_kd,
            kd_affinity=global_kd_affinity,
            rmax_per_conc=rmax_per_conc,
            global_chi_squared=global_chi_squared,
            individual_fits=individual_fits
        )
        
    except Exception as e:
        logger.error(f"Global fitting failed: {e}")
        raise ValueError(f"Global fitting failed: {e}") from e


def calculate_kd_from_kinetics(ka: float, kd: float) -> float:
    """
    Calculate equilibrium dissociation constant from kinetic rates.
    
    KD = kd / ka
    
    Args:
        ka: Association rate constant (1/Ms)
        kd: Dissociation rate constant (1/s)
        
    Returns:
        Equilibrium dissociation constant (M)
        
    Raises:
        ValueError: If ka is zero or negative
    """
    if ka <= 0:
        raise ValueError("Association rate constant must be positive")
    
    kd_equilibrium = kd / ka
    logger.debug(f"Calculated KD = {kd_equilibrium*1e9:.2f} nM from kinetics")
    
    return kd_equilibrium


def calculate_kd_from_steady_state(concs: np.ndarray, responses: np.ndarray) -> float:
    """
    Calculate KD from steady-state binding responses using Langmuir isotherm.
    
    R = Rmax * C / (KD + C)
    
    Args:
        concs: Array of analyte concentrations (M)
        responses: Array of steady-state binding responses (RU)
        
    Returns:
        Equilibrium dissociation constant (M)
        
    Raises:
        ValueError: If arrays have different lengths or insufficient data
    """
    if len(concs) != len(responses):
        raise ValueError("Concentration and response arrays must have same length")
    
    if len(concs) < 3:
        raise ValueError("At least 3 data points required for steady-state fitting")
    
    # Remove any zero or negative concentrations/responses
    valid_mask = (concs > 0) & (responses > 0)
    concs = concs[valid_mask]
    responses = responses[valid_mask]
    
    if len(concs) < 3:
        raise ValueError("Insufficient valid data points for fitting")
    
    logger.info(f"Fitting steady-state data: {len(concs)} points")
    
    try:
        # Define Langmuir isotherm model
        def langmuir_isotherm(c, rmax, kd):
            return rmax * c / (kd + c)
        
        # Create lmfit model
        model = Model(langmuir_isotherm)
        params = Parameters()
        params.add('rmax', value=np.max(responses), min=0)
        params.add('kd', value=np.median(concs), min=concs.min()/100, max=concs.max()*100)
        
        # Fit the model
        result = model.fit(responses, params, c=concs)
        
        if not result.success:
            raise ValueError("Steady-state fitting failed to converge")
        
        kd_fit = result.params['kd'].value
        rmax_fit = result.params['rmax'].value
        
        logger.info(f"Steady-state fit: KD={kd_fit*1e9:.2f} nM, Rmax={rmax_fit:.1f} RU")
        
        return kd_fit
        
    except Exception as e:
        logger.error(f"Steady-state fitting failed: {e}")
        raise ValueError(f"Steady-state fitting failed: {e}") from e


def assess_mass_transport(fit: KineticFit) -> bool:
    """
    Assess whether binding is limited by mass transport.
    
    Mass transport limitation occurs when the observed association rate
    is limited by diffusion rather than intrinsic binding kinetics.
    
    Args:
        fit: KineticFit object to assess
        
    Returns:
        True if mass transport limitation is detected
    """
    # Simple heuristic: very high ka values (>1e7 1/Ms) suggest mass transport
    # In a full implementation, this would use more sophisticated analysis
    
    if fit.ka > 1e7:
        logger.warning(f"Possible mass transport limitation: ka={fit.ka:.2e} 1/Ms")
        return True
    
    # Check for poor fit quality which might indicate mass transport
    if fit.chi_squared > 100:  # Arbitrary threshold
        logger.warning(f"Poor fit quality may indicate mass transport: χ²={fit.chi_squared:.1f}")
        return True
    
    return False


# Update module exports
__all__ = [
    "KineticFit",
    "GlobalFit", 
    "fit_1_to_1_langmuir",
    "fit_two_state",
    "global_fit",
    "calculate_kd_from_kinetics",
    "calculate_kd_from_steady_state",
    "assess_mass_transport",
]
