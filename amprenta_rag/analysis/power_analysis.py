"""Statistical power analysis utilities."""
from __future__ import annotations

from typing import Literal

try:
    from scipy.stats import power as scipy_power
    from scipy.stats.power import TTestIndPower, FTestAnovaPower, TTestPower
    SCIPY_AVAILABLE = True
except ImportError:
    SCIPY_AVAILABLE = False
    TTestIndPower = None  # type: ignore
    FTestAnovaPower = None  # type: ignore
    TTestPower = None  # type: ignore

from amprenta_rag.logging_utils import get_logger

logger = get_logger(__name__)

TestType = Literal["t-test", "anova", "correlation", "chi-square"]


def get_effect_size_preset(size: str) -> float:
    """
    Get preset effect size (Cohen's d conventions).
    
    Args:
        size: "small", "medium", or "large"
        
    Returns:
        Effect size value
    """
    presets = {
        "small": 0.2,
        "medium": 0.5,
        "large": 0.8,
    }
    return presets.get(size.lower(), 0.5)


def calculate_sample_size(
    effect_size: float,
    alpha: float = 0.05,
    power: float = 0.80,
    test_type: TestType = "t-test",
) -> int:
    """
    Calculate required sample size for a given effect size and power.
    
    Args:
        effect_size: Effect size (Cohen's d for t-test, f for ANOVA, r for correlation)
        alpha: Significance level (default 0.05)
        power: Desired power (default 0.80)
        test_type: Type of test ("t-test", "anova", "correlation", "chi-square")
        
    Returns:
        Required sample size (per group for t-test, total for others)
    """
    if not SCIPY_AVAILABLE:
        logger.error("[POWER] scipy.stats not available")
        raise ImportError("scipy.stats is required for power analysis")
    
    try:
        if test_type == "t-test":
            power_analysis = TTestIndPower()
            n = power_analysis.solve_power(
                effect_size=effect_size,
                alpha=alpha,
                power=power,
                ratio=1.0,  # Equal group sizes
            )
            return int(n) if n else 0
        
        elif test_type == "anova":
            power_analysis = FTestAnovaPower()
            # For ANOVA, effect_size is f (Cohen's f)
            n = power_analysis.solve_power(
                effect_size=effect_size,
                alpha=alpha,
                power=power,
                nobs=None,
            )
            return int(n) if n else 0
        
        elif test_type == "correlation":
            # For correlation, use TTestPower with transformed effect size
            # r to t transformation: t = r * sqrt((n-2)/(1-r^2))
            # We'll use an approximation
            power_analysis = TTestPower()
            # Approximate: for correlation, we need to solve iteratively
            # Using a simplified approach: treat r as effect size
            n = power_analysis.solve_power(
                effect_size=effect_size,
                alpha=alpha,
                power=power,
                nobs=None,
            )
            return int(n) if n else 0
        
        elif test_type == "chi-square":
            # Chi-square power analysis is more complex
            # Using a simplified approximation based on effect size
            # For chi-square, effect_size is typically Cohen's w
            # This is a simplified calculation
            from scipy.stats import chi2
            # Approximate calculation for chi-square
            # This is a placeholder - full implementation would require more parameters
            logger.warning("[POWER] Chi-square power analysis is simplified")
            power_analysis = TTestPower()  # Using t-test as approximation
            n = power_analysis.solve_power(
                effect_size=effect_size,
                alpha=alpha,
                power=power,
                nobs=None,
            )
            return int(n) if n else 0
        
        else:
            raise ValueError(f"Unsupported test type: {test_type}")
    
    except Exception as e:
        logger.error("[POWER] Error calculating sample size: %r", e)
        raise


def calculate_power(
    n: int,
    effect_size: float,
    alpha: float = 0.05,
    test_type: TestType = "t-test",
) -> float:
    """
    Calculate statistical power for a given sample size and effect size.
    
    Args:
        n: Sample size (per group for t-test, total for others)
        effect_size: Effect size (Cohen's d for t-test, f for ANOVA, r for correlation)
        alpha: Significance level (default 0.05)
        test_type: Type of test ("t-test", "anova", "correlation", "chi-square")
        
    Returns:
        Statistical power (0.0 to 1.0)
    """
    if not SCIPY_AVAILABLE:
        logger.error("[POWER] scipy.stats not available")
        raise ImportError("scipy.stats is required for power analysis")
    
    try:
        if test_type == "t-test":
            power_analysis = TTestIndPower()
            power = power_analysis.power(
                effect_size=effect_size,
                nobs1=n,
                alpha=alpha,
                ratio=1.0,
            )
            return float(power)
        
        elif test_type == "anova":
            power_analysis = FTestAnovaPower()
            power = power_analysis.power(
                effect_size=effect_size,
                nobs=n,
                alpha=alpha,
            )
            return float(power)
        
        elif test_type == "correlation":
            power_analysis = TTestPower()
            # For correlation, approximate using t-test power
            power = power_analysis.power(
                effect_size=effect_size,
                nobs=n,
                alpha=alpha,
            )
            return float(power)
        
        elif test_type == "chi-square":
            # Simplified chi-square power calculation
            logger.warning("[POWER] Chi-square power calculation is simplified")
            power_analysis = TTestPower()
            power = power_analysis.power(
                effect_size=effect_size,
                nobs=n,
                alpha=alpha,
            )
            return float(power)
        
        else:
            raise ValueError(f"Unsupported test type: {test_type}")
    
    except Exception as e:
        logger.error("[POWER] Error calculating power: %r", e)
        raise
