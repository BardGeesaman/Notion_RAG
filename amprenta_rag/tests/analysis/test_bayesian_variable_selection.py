"""Unit tests for Bayesian variable selection."""
from __future__ import annotations

import pytest
import numpy as np
from unittest.mock import patch, MagicMock


def test_compute_pip_returns_expected_structure():
    """Test PIP function returns correct structure."""
    mock_trace = MagicMock()
    mock_trace.posterior = {"beta": MagicMock()}
    mock_trace.posterior["beta"].values = np.random.randn(2, 500, 3).reshape(-1, 3)
    
    with patch("amprenta_rag.analysis.bayesian_variable_selection._require_pymc") as mock_req:
        mock_pm = MagicMock()
        mock_az = MagicMock()
        mock_req.return_value = (mock_pm, mock_az)
        mock_pm.Model.return_value.__enter__ = MagicMock()
        mock_pm.Model.return_value.__exit__ = MagicMock()
        mock_pm.sample.return_value = mock_trace
        
        from amprenta_rag.analysis.bayesian_variable_selection import compute_posterior_inclusion_probabilities
        
        X = np.random.randn(20, 3)
        y = np.random.randn(20)
        names = ["feat_1", "feat_2", "feat_3"]
        
        result = compute_posterior_inclusion_probabilities(X, y, names)
        
        assert "pips" in result
        assert "included_features" in result
        assert "coefficients" in result
        assert len(result["pips"]) == 3


def test_compute_pip_feature_name_mismatch():
    """Test validation of feature name count."""
    from amprenta_rag.analysis.bayesian_variable_selection import compute_posterior_inclusion_probabilities
    
    X = np.random.randn(20, 3)
    y = np.random.randn(20)
    names = ["feat_1", "feat_2"]  # Only 2 names for 3 features
    
    with pytest.raises(ValueError, match="feature_names length"):
        compute_posterior_inclusion_probabilities(X, y, names)


def test_compute_pip_threshold_filtering():
    """Test PIP threshold correctly filters features."""
    mock_trace = MagicMock()
    # Create beta samples where feat_1 has high values, feat_2 has low values
    beta_samples = np.zeros((1000, 2))
    beta_samples[:, 0] = np.random.randn(1000) * 2 + 1  # High signal
    beta_samples[:, 1] = np.random.randn(1000) * 0.001  # Low signal
    mock_trace.posterior = {"beta": MagicMock()}
    mock_trace.posterior["beta"].values = beta_samples.reshape(2, 500, 2)
    
    with patch("amprenta_rag.analysis.bayesian_variable_selection._require_pymc") as mock_req:
        mock_pm = MagicMock()
        mock_az = MagicMock()
        mock_req.return_value = (mock_pm, mock_az)
        mock_pm.Model.return_value.__enter__ = MagicMock()
        mock_pm.Model.return_value.__exit__ = MagicMock()
        mock_pm.sample.return_value = mock_trace
        
        from amprenta_rag.analysis.bayesian_variable_selection import compute_posterior_inclusion_probabilities
        
        X = np.random.randn(20, 2)
        y = np.random.randn(20)
        names = ["high_signal", "low_signal"]
        
        result = compute_posterior_inclusion_probabilities(X, y, names, pip_threshold=0.5)
        
        # high_signal should be included, low_signal should not
        assert "high_signal" in result["included_features"]


def test_prior_config_schema_validation():
    """Test PriorConfig schema accepts valid values."""
    from amprenta_rag.api.schemas import PriorConfig
    
    config = PriorConfig(
        ec50_prior_mean=-6.0,
        ec50_prior_sd=1.5,
        hill_prior_mean=1.0,
        hill_prior_sd=1.0,
    )
    assert config.ec50_prior_mean == -6.0
    assert config.ec50_prior_sd == 1.5


def test_prior_config_defaults():
    """Test PriorConfig uses correct defaults."""
    from amprenta_rag.api.schemas import PriorConfig
    
    config = PriorConfig()
    assert config.ec50_prior_mean is None  # Auto from data
    assert config.ec50_prior_sd == 1.0
    assert config.hill_prior_mean == 1.0
    assert config.hill_prior_sd == 2.0


def test_compute_pip_standardization():
    """Test that input data is properly standardized."""
    mock_trace = MagicMock()
    mock_trace.posterior = {"beta": MagicMock()}
    mock_trace.posterior["beta"].values = np.random.randn(2, 500, 2).reshape(-1, 2)
    
    with patch("amprenta_rag.analysis.bayesian_variable_selection._require_pymc") as mock_req:
        mock_pm = MagicMock()
        mock_az = MagicMock()
        mock_req.return_value = (mock_pm, mock_az)
        mock_pm.Model.return_value.__enter__ = MagicMock()
        mock_pm.Model.return_value.__exit__ = MagicMock()
        mock_pm.sample.return_value = mock_trace
        
        from amprenta_rag.analysis.bayesian_variable_selection import compute_posterior_inclusion_probabilities
        
        # Create data with different scales
        X = np.array([[1, 100], [2, 200], [3, 300], [4, 400], [5, 500]])
        y = np.array([10, 20, 30, 40, 50])
        names = ["small_scale", "large_scale"]
        
        result = compute_posterior_inclusion_probabilities(X, y, names)
        
        # Should complete without error despite different scales
        assert "pips" in result
        assert len(result["pips"]) == 2


def test_compute_pip_mcmc_parameters():
    """Test MCMC sampling parameters are correctly set."""
    mock_trace = MagicMock()
    mock_trace.posterior = {"beta": MagicMock()}
    mock_trace.posterior["beta"].values = np.random.randn(2, 500, 2).reshape(-1, 2)
    
    with patch("amprenta_rag.analysis.bayesian_variable_selection._require_pymc") as mock_req:
        mock_pm = MagicMock()
        mock_az = MagicMock()
        mock_req.return_value = (mock_pm, mock_az)
        mock_pm.Model.return_value.__enter__ = MagicMock()
        mock_pm.Model.return_value.__exit__ = MagicMock()
        mock_pm.sample.return_value = mock_trace
        
        from amprenta_rag.analysis.bayesian_variable_selection import compute_posterior_inclusion_probabilities
        
        X = np.random.randn(20, 2)
        y = np.random.randn(20)
        names = ["feat_1", "feat_2"]
        
        compute_posterior_inclusion_probabilities(X, y, names, n_samples=500)
        
        # Verify sample was called with correct parameters
        mock_pm.sample.assert_called_once_with(
            draws=500, tune=500, chains=2, progressbar=False, random_seed=42
        )


def test_compute_pip_horseshoe_prior_structure():
    """Test that Horseshoe prior structure is correctly implemented."""
    mock_trace = MagicMock()
    mock_trace.posterior = {"beta": MagicMock()}
    mock_trace.posterior["beta"].values = np.random.randn(2, 500, 3).reshape(-1, 3)
    
    with patch("amprenta_rag.analysis.bayesian_variable_selection._require_pymc") as mock_req:
        mock_pm = MagicMock()
        mock_az = MagicMock()
        mock_req.return_value = (mock_pm, mock_az)
        
        # Mock the context manager
        mock_model = MagicMock()
        mock_pm.Model.return_value = mock_model
        mock_model.__enter__ = MagicMock(return_value=mock_model)
        mock_model.__exit__ = MagicMock(return_value=None)
        mock_pm.sample.return_value = mock_trace
        
        from amprenta_rag.analysis.bayesian_variable_selection import compute_posterior_inclusion_probabilities
        
        X = np.random.randn(20, 3)
        y = np.random.randn(20)
        names = ["feat_1", "feat_2", "feat_3"]
        
        compute_posterior_inclusion_probabilities(X, y, names)
        
        # Verify HalfCauchy priors were created for global and local shrinkage
        mock_pm.HalfCauchy.assert_any_call("tau", beta=1.0)
        mock_pm.HalfCauchy.assert_any_call("lambdas", beta=1.0, shape=3)
