"""Statistical feature selection methods for biomarker discovery."""

from typing import List, Tuple, Optional
import numpy as np

try:
    from scipy import stats
    SCIPY_AVAILABLE = True
except ImportError:
    SCIPY_AVAILABLE = False

try:
    from statsmodels.stats.multitest import multipletests
    STATSMODELS_AVAILABLE = True
except ImportError:
    STATSMODELS_AVAILABLE = False


def t_test_selection(
    X: np.ndarray,
    y: np.ndarray,
    feature_names: Optional[List[str]] = None,
) -> List[Tuple[str, float, float]]:
    """
    Two-sample t-test for binary classification feature selection.
    
    Args:
        X: Feature matrix (n_samples, n_features)
        y: Binary labels (0/1)
        feature_names: Optional feature names
        
    Returns:
        List of (feature_name, t_statistic, p_value) sorted by p-value
    """
    if not SCIPY_AVAILABLE:
        raise ImportError("scipy is required for t_test_selection")
    
    X = np.asarray(X)
    y = np.asarray(y)
    
    if feature_names is None:
        feature_names = [f"feature_{i}" for i in range(X.shape[1])]
    
    unique_classes = np.unique(y)
    if len(unique_classes) != 2:
        raise ValueError(f"t_test_selection requires binary labels, got {len(unique_classes)} classes")
    
    class_0_mask = y == unique_classes[0]
    class_1_mask = y == unique_classes[1]
    
    results = []
    for i in range(X.shape[1]):
        group_0 = X[class_0_mask, i]
        group_1 = X[class_1_mask, i]
        
        # Handle constant features
        if np.std(group_0) == 0 and np.std(group_1) == 0:
            t_stat, p_val = 0.0, 1.0
        else:
            t_stat, p_val = stats.ttest_ind(group_0, group_1, equal_var=False)
        
        results.append((feature_names[i], float(t_stat), float(p_val)))
    
    # Sort by p-value
    results.sort(key=lambda x: x[2])
    return results


def mann_whitney_selection(
    X: np.ndarray,
    y: np.ndarray,
    feature_names: Optional[List[str]] = None,
) -> List[Tuple[str, float, float]]:
    """
    Mann-Whitney U test (non-parametric) for binary classification.
    
    Args:
        X: Feature matrix (n_samples, n_features)
        y: Binary labels (0/1)
        feature_names: Optional feature names
        
    Returns:
        List of (feature_name, u_statistic, p_value) sorted by p-value
    """
    if not SCIPY_AVAILABLE:
        raise ImportError("scipy is required for mann_whitney_selection")
    
    X = np.asarray(X)
    y = np.asarray(y)
    
    if feature_names is None:
        feature_names = [f"feature_{i}" for i in range(X.shape[1])]
    
    unique_classes = np.unique(y)
    if len(unique_classes) != 2:
        raise ValueError(f"mann_whitney_selection requires binary labels, got {len(unique_classes)} classes")
    
    class_0_mask = y == unique_classes[0]
    class_1_mask = y == unique_classes[1]
    
    results = []
    for i in range(X.shape[1]):
        group_0 = X[class_0_mask, i]
        group_1 = X[class_1_mask, i]
        
        try:
            u_stat, p_val = stats.mannwhitneyu(group_0, group_1, alternative='two-sided')
        except ValueError:
            # Handle edge cases (all identical values)
            u_stat, p_val = 0.0, 1.0
        
        results.append((feature_names[i], float(u_stat), float(p_val)))
    
    # Sort by p-value
    results.sort(key=lambda x: x[2])
    return results


def anova_selection(
    X: np.ndarray,
    y: np.ndarray,
    feature_names: Optional[List[str]] = None,
) -> List[Tuple[str, float, float]]:
    """
    One-way ANOVA for multi-class feature selection.
    
    Args:
        X: Feature matrix (n_samples, n_features)
        y: Class labels (can be multi-class)
        feature_names: Optional feature names
        
    Returns:
        List of (feature_name, f_statistic, p_value) sorted by p-value
    """
    if not SCIPY_AVAILABLE:
        raise ImportError("scipy is required for anova_selection")
    
    X = np.asarray(X)
    y = np.asarray(y)
    
    if feature_names is None:
        feature_names = [f"feature_{i}" for i in range(X.shape[1])]
    
    unique_classes = np.unique(y)
    if len(unique_classes) < 2:
        raise ValueError("anova_selection requires at least 2 classes")
    
    results = []
    for i in range(X.shape[1]):
        groups = [X[y == c, i] for c in unique_classes]
        
        # Filter out empty groups
        groups = [g for g in groups if len(g) > 0]
        
        if len(groups) < 2:
            f_stat, p_val = 0.0, 1.0
        else:
            try:
                f_stat, p_val = stats.f_oneway(*groups)
            except ValueError:
                f_stat, p_val = 0.0, 1.0
        
        # Handle NaN results
        if np.isnan(f_stat):
            f_stat = 0.0
        if np.isnan(p_val):
            p_val = 1.0
            
        results.append((feature_names[i], float(f_stat), float(p_val)))
    
    # Sort by p-value
    results.sort(key=lambda x: x[2])
    return results


def fdr_correction(
    p_values: List[float],
    alpha: float = 0.05,
    method: str = "fdr_bh",
) -> Tuple[List[bool], List[float]]:
    """
    Apply FDR correction (Benjamini-Hochberg) to p-values.
    
    Args:
        p_values: List of p-values
        alpha: Significance threshold (default 0.05)
        method: Correction method ('fdr_bh' for Benjamini-Hochberg)
        
    Returns:
        Tuple of (reject array, corrected p-values)
    """
    if not STATSMODELS_AVAILABLE:
        raise ImportError("statsmodels is required for fdr_correction")
    
    p_values = np.asarray(p_values)
    
    # Handle edge cases
    if len(p_values) == 0:
        return [], []
    
    reject, pvals_corrected, _, _ = multipletests(
        p_values, alpha=alpha, method=method
    )
    
    return list(reject), list(pvals_corrected)

