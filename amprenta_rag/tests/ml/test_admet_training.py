"""Tests for ADMET training pipeline components."""

from __future__ import annotations

import subprocess
import sys
from pathlib import Path
import numpy as np


def test_download_script_dry_run():
    """Test that download script exists and dry-run works."""
    script_path = Path(__file__).parents[3] / "scripts" / "download_chembl_sqlite.py"
    assert script_path.exists(), f"Download script not found: {script_path}"
    
    # Test dry-run functionality
    result = subprocess.run([
        sys.executable, str(script_path), "--dry-run"
    ], capture_output=True, text=True, timeout=30)
    
    assert result.returncode == 0, f"Dry-run failed: {result.stderr}"
    assert "DRY RUN: ChEMBL SQLite download" in result.stdout
    assert "chembl_34" in result.stdout


def test_train_script_dry_run():
    """Test that training script dry-run lists all 9 endpoints."""
    script_path = Path(__file__).parents[3] / "scripts" / "train_admet_models.py"
    assert script_path.exists(), f"Training script not found: {script_path}"
    
    # Test dry-run functionality
    result = subprocess.run([
        sys.executable, str(script_path), "--dry-run"
    ], capture_output=True, text=True, timeout=30)
    
    assert result.returncode == 0, f"Dry-run failed: {result.stderr}"
    assert "Training endpoints: ['herg', 'logs', 'logp', 'cyp3a4', 'cyp2d6', 'cyp2c9', 'bbb', 'caco2', 'clearance']" in result.stdout
    
    # Check all 9 endpoints are listed
    expected_endpoints = ["herg", "logs", "logp", "cyp3a4", "cyp2d6", "cyp2c9", "bbb", "caco2", "clearance"]
    for endpoint in expected_endpoints:
        assert endpoint in result.stdout, f"Endpoint {endpoint} not found in dry-run output"


def test_admet_models_dict_has_9_endpoints():
    """Test that ADMET_MODELS dict has all 9 endpoints."""
    from amprenta_rag.ml.admet.predictor import ADMET_MODELS
    
    expected_endpoints = {"herg", "logs", "logp", "cyp3a4", "cyp2d6", "cyp2c9", "bbb", "caco2", "clearance"}
    actual_endpoints = set(ADMET_MODELS.keys())
    
    assert len(ADMET_MODELS) == 9, f"Expected 9 endpoints, got {len(ADMET_MODELS)}"
    assert actual_endpoints == expected_endpoints, f"Endpoint mismatch. Expected: {expected_endpoints}, Got: {actual_endpoints}"
    
    # Check all model names follow ensemble pattern
    for endpoint, model_name in ADMET_MODELS.items():
        expected_name = f"admet_{endpoint}_ensemble"
        assert model_name == expected_name, f"Model name mismatch for {endpoint}: expected {expected_name}, got {model_name}"


def test_classification_endpoints_correct():
    """Test that CLASSIFICATION_ENDPOINTS contains correct endpoints."""
    from amprenta_rag.ml.admet.predictor import CLASSIFICATION_ENDPOINTS
    
    expected_classification = {"herg", "cyp3a4", "cyp2d6", "cyp2c9", "bbb"}
    actual_classification = set(CLASSIFICATION_ENDPOINTS)
    
    assert actual_classification == expected_classification, f"Classification endpoints mismatch. Expected: {expected_classification}, Got: {actual_classification}"


def test_bootstrap_ensemble_fit():
    """Test BootstrapEnsemble class can be instantiated and has expected attributes."""
    from amprenta_rag.ml.admet.ensemble import BootstrapEnsemble
    
    # Test ensemble initialization
    ens = BootstrapEnsemble(n_models=3)
    
    assert ens.n_models == 3
    assert ens.models == []
    assert ens.task is None
    assert ens.training_centroid is None
    assert ens.feature_dim is None
    
    # Test with base_params
    base_params = {"max_depth": 6, "n_estimators": 100}
    ens_with_params = BootstrapEnsemble(n_models=5, base_params=base_params)
    
    assert ens_with_params.n_models == 5
    assert ens_with_params.base_params == base_params


def test_calibration_wrapper_fit():
    """Test CalibrationWrapper with mock probabilities."""
    from amprenta_rag.ml.admet.calibration import CalibrationWrapper
    
    # Create mock calibration data
    n_samples = 100
    y_true = np.random.randint(0, 2, n_samples)
    y_prob = np.random.rand(n_samples)
    
    # Test isotonic calibration
    cal = CalibrationWrapper(method="isotonic")
    cal.fit(y_prob, y_true)
    
    assert cal.method == "isotonic"
    assert hasattr(cal, "calibrator")
    
    # Test calibration
    y_prob_test = np.random.rand(50)
    y_calibrated = cal.calibrate(y_prob_test)
    
    assert len(y_calibrated) == 50
    assert np.all((y_calibrated >= 0) & (y_calibrated <= 1)), "Calibrated probabilities should be in [0, 1]"


def test_applicability_checker_fit():
    """Test ApplicabilityChecker domain checking."""
    from amprenta_rag.ml.admet.applicability import ApplicabilityChecker
    
    # Create mock training data
    X_train = np.random.rand(100, 50)
    
    # Test applicability checker
    app = ApplicabilityChecker(threshold=0.3)
    app.fit(X_train)
    
    assert app.threshold == 0.3
    assert hasattr(app, "training_centroid")
    assert app.training_centroid.shape == (50,)
    
    # Test applicability checking
    X_test = np.random.rand(20, 50)
    is_applicable, similarities = app.check(X_test)
    
    assert len(is_applicable) == 20
    assert len(similarities) == 20
    assert all(isinstance(x, (bool, np.bool_)) for x in is_applicable)
    assert all(isinstance(x, (float, np.floating)) for x in similarities)


def test_registry_model_naming():
    """Test that model names follow admet_{endpoint}_ensemble pattern."""
    from amprenta_rag.ml.admet.predictor import ADMET_MODELS
    
    # Test naming convention
    for endpoint, model_name in ADMET_MODELS.items():
        # Should follow pattern: admet_{endpoint}_ensemble
        expected_pattern = f"admet_{endpoint}_ensemble"
        assert model_name == expected_pattern, f"Model name {model_name} doesn't follow pattern {expected_pattern}"
        
        # Should not be legacy XGB names
        legacy_pattern = f"admet_{endpoint}_xgb"
        assert model_name != legacy_pattern, f"Model name {model_name} is using legacy XGB pattern"
    
    # Test that we have both classification and regression endpoints
    from amprenta_rag.ml.admet.predictor import CLASSIFICATION_ENDPOINTS
    
    classification_models = {ep for ep in ADMET_MODELS.keys() if ep in CLASSIFICATION_ENDPOINTS}
    regression_models = {ep for ep in ADMET_MODELS.keys() if ep not in CLASSIFICATION_ENDPOINTS}
    
    assert len(classification_models) == 5, f"Expected 5 classification endpoints, got {len(classification_models)}"
    assert len(regression_models) == 4, f"Expected 4 regression endpoints, got {len(regression_models)}"
    
    # Verify specific endpoint types
    expected_classification = {"herg", "cyp3a4", "cyp2d6", "cyp2c9", "bbb"}
    expected_regression = {"logs", "logp", "caco2", "clearance"}
    
    assert classification_models == expected_classification
    assert regression_models == expected_regression
