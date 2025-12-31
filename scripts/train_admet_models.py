#!/usr/bin/env python3
"""Train ADMET ensemble models for all 9 endpoints using TDC datasets.

This script trains BootstrapEnsemble models with calibration and applicability checking
for comprehensive ADMET property prediction.
"""

from __future__ import annotations

import argparse
import sys
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Dict, List, Optional


def _ensure_repo_on_path() -> None:
    """Ensure repo root is importable when running as script."""
    repo_root = Path(__file__).resolve().parents[1]
    if str(repo_root) not in sys.path:
        sys.path.insert(0, str(repo_root))


def _require_training_deps():
    """Import all required training dependencies."""
    try:
        import numpy as np  # type: ignore
        from sklearn.metrics import (  # type: ignore
            brier_score_loss,
            mean_squared_error,
            r2_score,
            roc_auc_score,
        )
        from tdc.single_pred import ADME, Tox  # type: ignore
    except Exception as e:  # noqa: BLE001
        raise ImportError(
            "Training requires: tdc, scikit-learn, numpy (and rdkit for feature extraction). "
            f"Underlying error: {type(e).__name__}: {e}"
        ) from e
    return np, roc_auc_score, brier_score_loss, mean_squared_error, r2_score, ADME, Tox


def _require_ml_modules():
    """Import ML modules from amprenta_rag."""
    try:
        from amprenta_rag.ml.admet.ensemble import BootstrapEnsemble  # type: ignore
        from amprenta_rag.ml.admet.calibration import CalibrationWrapper  # type: ignore
        from amprenta_rag.ml.admet.applicability import ApplicabilityChecker  # type: ignore
        from amprenta_rag.ml.admet.predictor import ADMETPredictor  # type: ignore
        from amprenta_rag.ml.registry import get_registry  # type: ignore
    except Exception as e:  # noqa: BLE001
        raise ImportError(
            f"Could not import amprenta_rag ML modules: {type(e).__name__}: {e}"
        ) from e
    return BootstrapEnsemble, CalibrationWrapper, ApplicabilityChecker, ADMETPredictor, get_registry


def _now_iso() -> str:
    """Get current timestamp in ISO format."""
    return datetime.now(timezone.utc).isoformat()


def _featurize_smiles(predictor, smiles_list: List[str]):
    """Convert SMILES to feature vectors, handling failures gracefully."""
    import numpy as np

    feats = [predictor._get_features(s) for s in smiles_list]  # noqa: SLF001
    mask = [x is not None for x in feats]
    X = np.array([x for x, ok in zip(feats, mask) if ok])
    return X, mask


def _coverage_90(y_true, y_pred, y_std):
    """Compute 90% prediction interval coverage."""
    import numpy as np
    
    lower = y_pred - 1.645 * y_std
    upper = y_pred + 1.645 * y_std
    coverage = np.mean((y_true >= lower) & (y_true <= upper))
    return coverage


def _reliability_diagram(y_prob, y_true, n_bins=10):
    """Compute reliability diagram and ECE for calibration assessment."""
    import numpy as np
    
    bin_boundaries = np.linspace(0, 1, n_bins + 1)
    bin_lowers = bin_boundaries[:-1]
    bin_uppers = bin_boundaries[1:]
    
    ece = 0
    bin_accuracies = []
    bin_confidences = []
    bin_counts = []
    
    for bin_lower, bin_upper in zip(bin_lowers, bin_uppers):
        in_bin = (y_prob > bin_lower) & (y_prob <= bin_upper)
        prop_in_bin = in_bin.mean()
        
        if prop_in_bin > 0:
            accuracy_in_bin = y_true[in_bin].mean()
            avg_confidence_in_bin = y_prob[in_bin].mean()
            ece += np.abs(avg_confidence_in_bin - accuracy_in_bin) * prop_in_bin
            
            bin_accuracies.append(accuracy_in_bin)
            bin_confidences.append(avg_confidence_in_bin)
            bin_counts.append(in_bin.sum())
        else:
            bin_accuracies.append(0)
            bin_confidences.append(0)
            bin_counts.append(0)
    
    return {
        "ece": ece,
        "bin_accuracies": bin_accuracies,
        "bin_confidences": bin_confidences,
        "bin_counts": bin_counts,
    }


# ADMET endpoint configurations
ADMET_ENDPOINTS = {
    "herg": {
        "tdc_class": "Tox",
        "tdc_name": "hERG",
        "type": "classification",
        "description": "hERG channel inhibition classifier"
    },
    "logs": {
        "tdc_class": "ADME",
        "tdc_name": "Solubility_AqSolDB", 
        "type": "regression",
        "description": "Aqueous solubility regressor"
    },
    "logp": {
        "tdc_class": "ADME",
        "tdc_name": "Lipophilicity_AstraZeneca",
        "type": "regression", 
        "description": "Lipophilicity regressor"
    },
    "cyp3a4": {
        "tdc_class": "Tox",
        "tdc_name": "CYP3A4_Veith",
        "type": "classification",
        "description": "CYP3A4 inhibition classifier"
    },
    "cyp2d6": {
        "tdc_class": "Tox", 
        "tdc_name": "CYP2D6_Veith",
        "type": "classification",
        "description": "CYP2D6 inhibition classifier"
    },
    "cyp2c9": {
        "tdc_class": "Tox",
        "tdc_name": "CYP2C9_Veith", 
        "type": "classification",
        "description": "CYP2C9 inhibition classifier"
    },
    "bbb": {
        "tdc_class": "ADME",
        "tdc_name": "BBB_Martins",
        "type": "classification",
        "description": "Blood-brain barrier penetration classifier"
    },
    "caco2": {
        "tdc_class": "ADME",
        "tdc_name": "Caco2_Wang",
        "type": "regression",
        "description": "Caco-2 permeability regressor"
    },
    "clearance": {
        "tdc_class": "ADME",
        "tdc_name": "Clearance_Hepatocyte_AZ",
        "type": "regression", 
        "description": "Hepatic clearance regressor"
    }
}


def train_endpoint(endpoint: str, n_models: int = 5, dry_run: bool = False) -> Optional[Dict[str, Any]]:
    """Train a single ADMET endpoint with bootstrap ensemble."""
    if endpoint not in ADMET_ENDPOINTS:
        raise ValueError(f"Unknown endpoint: {endpoint}. Available: {list(ADMET_ENDPOINTS.keys())}")
    
    config = ADMET_ENDPOINTS[endpoint]
    is_classification = config["type"] == "classification"
    
    print(f"\n=== Training {endpoint} ({config['type']}) ===")
    print(f"TDC Dataset: {config['tdc_class']}('{config['tdc_name']}')")
    
    if dry_run:
        print("DRY RUN: Would load TDC dataset and train ensemble")
        return None
    
    # Import dependencies
    np, roc_auc_score, brier_score_loss, mean_squared_error, r2_score, ADME, Tox = _require_training_deps()
    BootstrapEnsemble, CalibrationWrapper, ApplicabilityChecker, ADMETPredictor, get_registry = _require_ml_modules()
    
    # Load TDC dataset
    tdc_class = Tox if config["tdc_class"] == "Tox" else ADME
    dataset = tdc_class(name=config["tdc_name"])
    train, val, test = dataset.get_split()
    
    print(f"Dataset splits: train={len(train)}, val={len(val)}, test={len(test)}")
    
    # Initialize predictor for feature extraction
    predictor = ADMETPredictor()
    
    # Featurize all splits
    X_train, mask_train = _featurize_smiles(predictor, train["Drug"].tolist())
    y_train = train["Y"].values[mask_train]
    if not is_classification:
        y_train = y_train.astype(float)
    
    X_val, mask_val = _featurize_smiles(predictor, val["Drug"].tolist()) 
    y_val = val["Y"].values[mask_val]
    if not is_classification:
        y_val = y_val.astype(float)
        
    X_test, mask_test = _featurize_smiles(predictor, test["Drug"].tolist())
    y_test = test["Y"].values[mask_test]
    if not is_classification:
        y_test = y_test.astype(float)
    
    print(f"Features: train={X_train.shape}, val={X_val.shape}, test={X_test.shape}")
    
    # Train bootstrap ensemble
    print(f"Training BootstrapEnsemble (n_models={n_models})...")
    ens = BootstrapEnsemble(n_models=n_models)
    ens.fit(X_train, y_train)
    
    # Calibration and metrics
    calibrator_obj: Optional[Any] = None
    
    if is_classification:
        # Classification: use validation set for calibration
        p_val, p_val_std = ens.predict_proba(X_val)
        p_test, p_test_std = ens.predict_proba(X_test)
        
        # Calibrate probabilities
        cal = CalibrationWrapper(method="isotonic")
        cal.fit(p_val, y_val)
        p_test_cal = cal.calibrate(p_test)
        calibrator_obj = cal
        
        # Compute metrics
        auc = float(roc_auc_score(y_test, p_test_cal))
        brier = float(brier_score_loss(y_test, p_test_cal))
        ece = float(_reliability_diagram(p_test_cal, y_test, n_bins=10)["ece"])
        coverage = float(_coverage_90(y_test, p_test_cal, p_test_std))
        
        metrics = {
            "auc": auc,
            "brier": brier, 
            "ece": ece,
            "coverage_90": coverage
        }
        
        print(f"Results: AUC={auc:.3f}, Brier={brier:.3f}, ECE={ece:.3f}, Coverage={coverage:.3f}")
        
    else:
        # Regression: direct predictions
        yhat_test, yhat_test_std = ens.predict(X_test)
        
        # Compute metrics
        rmse = float(np.sqrt(mean_squared_error(y_test, yhat_test)))
        r2 = float(r2_score(y_test, yhat_test))
        coverage = float(_coverage_90(y_test, yhat_test, yhat_test_std))
        
        metrics = {
            "rmse": rmse,
            "r2": r2,
            "coverage_90": coverage
        }
        
        print(f"Results: RMSE={rmse:.3f}, R²={r2:.3f}, Coverage={coverage:.3f}")
    
    # Applicability checker
    app = ApplicabilityChecker(threshold=0.3)
    app.fit(X_train)
    
    # Create artifact
    artifact = {
        "ensemble": ens.to_artifact(),
        "calibrator": calibrator_obj,
        "applicability": {
            "centroid": app.training_centroid, 
            "threshold": float(app.threshold)
        },
        "feature_schema": {"morgan_bits": 2048, "rdkit_descs": 6},
        "metrics": metrics,
        "endpoint": endpoint,
        "trained_at": _now_iso(),
    }
    
    # Register with ML registry
    registry = get_registry()
    model_type = "admet_classification" if is_classification else "admet_regression"
    
    registry.register_model(
        name=f"admet_{endpoint}_ensemble",
        version="1.0.0",
        model_type=model_type,
        framework="xgboost_ensemble",
        model_object=artifact,
        metrics=metrics,
        description=f"Bootstrap ensemble + calibration + applicability for {config['description']}",
    )
    
    print(f"✓ Registered model: admet_{endpoint}_ensemble")
    
    return artifact


def main() -> None:
    """Main CLI entry point."""
    _ensure_repo_on_path()
    
    parser = argparse.ArgumentParser(
        description="Train ADMET ensemble models for TDC endpoints",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    
    parser.add_argument(
        "--endpoints",
        default=",".join(ADMET_ENDPOINTS.keys()),
        help="Comma-separated list of endpoints to train"
    )
    parser.add_argument(
        "--n-models", 
        type=int,
        default=5,
        help="Number of models in bootstrap ensemble"
    )
    parser.add_argument(
        "--dry-run",
        action="store_true", 
        help="Print configuration without training"
    )
    
    args = parser.parse_args()
    
    # Parse endpoints
    if args.endpoints.lower() == "all":
        endpoints = list(ADMET_ENDPOINTS.keys())
    else:
        endpoints = [ep.strip() for ep in args.endpoints.split(",")]
    
    # Validate endpoints
    invalid = [ep for ep in endpoints if ep not in ADMET_ENDPOINTS]
    if invalid:
        print(f"Error: Unknown endpoints: {invalid}")
        print(f"Available endpoints: {list(ADMET_ENDPOINTS.keys())}")
        return 1
    
    print(f"Training endpoints: {endpoints}")
    print(f"Ensemble size: {args.n_models}")
    print(f"Dry run: {args.dry_run}")
    
    if args.dry_run:
        print("\nDRY RUN: Endpoint configurations:")
        for endpoint in endpoints:
            config = ADMET_ENDPOINTS[endpoint]
            print(f"  {endpoint}: {config['tdc_class']}('{config['tdc_name']}') - {config['type']}")
        return 0
    
    # Train each endpoint
    try:
        for endpoint in endpoints:
            train_endpoint(endpoint, n_models=args.n_models, dry_run=args.dry_run)
        
        print(f"\n✓ Successfully trained {len(endpoints)} ADMET endpoints!")
        
    except Exception as e:
        print(f"\n✗ Error during training: {e}")
        return 1
    
    return 0


if __name__ == "__main__":
    exit(main())