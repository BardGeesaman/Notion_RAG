"""Statistical power analysis utilities."""
from __future__ import annotations

from typing import Literal

try:
    from statsmodels.stats.power import TTestIndPower, FTestAnovaPower, TTestPower
    STATSMODELS_AVAILABLE = True
except ImportError:
    STATSMODELS_AVAILABLE = False
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
    if not STATSMODELS_AVAILABLE:
        logger.error("[POWER] statsmodels not available")
        raise ImportError("statsmodels is required for power analysis")

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
    if not STATSMODELS_AVAILABLE:
        logger.error("[POWER] statsmodels not available")
        raise ImportError("statsmodels is required for power analysis")

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


def estimate_effect_size_from_data(group1: list, group2: list) -> float:
    """
    Calculate Cohen's d effect size from two groups of data.

    Args:
        group1: First group measurements
        group2: Second group measurements

    Returns:
        Cohen's d effect size (standardized mean difference)
    """
    import numpy as np
    
    if not group1 or not group2:
        return 0.0
    
    arr1 = np.array(group1)
    arr2 = np.array(group2)
    
    mean1 = np.mean(arr1)
    mean2 = np.mean(arr2)
    
    # Calculate pooled standard deviation
    var1 = np.var(arr1, ddof=1) if len(arr1) > 1 else 0.0
    var2 = np.var(arr2, ddof=1) if len(arr2) > 1 else 0.0
    
    pooled_std = np.sqrt((var1 + var2) / 2)
    
    if pooled_std == 0:
        return 0.0
    
    cohen_d = (mean1 - mean2) / pooled_std
    
    return float(cohen_d)


def calculate_plate_layout(n: int, plate_format: int = 96) -> dict:
    """
    Calculate plate layout requirements for n samples.

    Args:
        n: Number of samples
        plate_format: Wells per plate (96, 384, or 1536)

    Returns:
        Dictionary with plates_needed, wells_used, empty_wells
    """
    if n <= 0:
        return {"plates_needed": 0, "wells_used": 0, "empty_wells": 0}
    
    if plate_format not in [96, 384, 1536]:
        raise ValueError(f"Invalid plate format: {plate_format}. Use 96, 384, or 1536")
    
    plates_needed = (n + plate_format - 1) // plate_format
    wells_used = n
    total_wells = plates_needed * plate_format
    empty_wells = total_wells - wells_used
    
    return {
        "plates_needed": plates_needed,
        "wells_used": wells_used,
        "empty_wells": empty_wells,
    }


def estimate_experiment_cost(
    n: int,
    cost_per_sample: float,
    overhead_pct: float = 0.1,
) -> dict:
    """
    Estimate total experiment cost with overhead.

    Args:
        n: Number of samples
        cost_per_sample: Cost per sample in dollars
        overhead_pct: Overhead percentage (default 10%)

    Returns:
        Dictionary with sample_cost, overhead, total
    """
    sample_cost = n * cost_per_sample
    overhead = sample_cost * overhead_pct
    total = sample_cost + overhead
    
    return {
        "sample_cost": round(sample_cost, 2),
        "overhead": round(overhead, 2),
        "total": round(total, 2),
    }
