"""ADMET Predictor service."""

from __future__ import annotations

from typing import Any, Dict, List, Optional

import numpy as np

from amprenta_rag.logging_utils import get_logger
from amprenta_rag.ml.registry import get_registry

logger = get_logger(__name__)

# ADMET model names
ADMET_MODELS = {
    "herg": "admet_herg_ensemble",
    "logs": "admet_logs_ensemble",
    "logp": "admet_logp_ensemble",
    "cyp3a4": "admet_cyp3a4_ensemble",
    "cyp2d6": "admet_cyp2d6_ensemble",
    "cyp2c9": "admet_cyp2c9_ensemble",
    "bbb": "admet_bbb_ensemble",
    "caco2": "admet_caco2_ensemble",
    "clearance": "admet_clearance_ensemble",
}

# Backward-compatible legacy names (pre-ensemble)
_LEGACY_ADMET_MODELS = {
    "herg": "admet_herg_xgb",
    "logs": "admet_logs_xgb",
    "logp": "admet_logp_xgb",
    # New endpoints don't have legacy names
}

# Classification endpoints (for task type inference)
CLASSIFICATION_ENDPOINTS = {"herg", "cyp3a4", "cyp2d6", "cyp2c9", "bbb"}


def _require_rdkit():
    try:
        from rdkit import Chem
        from rdkit.Chem import AllChem, Descriptors
    except ImportError as e:
        raise ImportError("RDKit required for ADMET prediction") from e
    return Chem, Descriptors, AllChem


class ADMETPredictor:
    """ADMET prediction service using registered models."""

    def __init__(self):
        self.registry = get_registry()

    def _get_features(self, smiles: str) -> Optional[np.ndarray]:
        """Calculate molecular features from SMILES."""
        Chem, Descriptors, AllChem = _require_rdkit()

        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            return None

        # Use Morgan fingerprint (2048 bits) + basic descriptors
        fp = AllChem.GetMorganFingerprintAsBitVect(mol, radius=2, nBits=2048)
        fp_arr = np.array(fp, dtype=np.float32)

        # Add basic RDKit descriptors
        basic_descs = np.array(
            [
                Descriptors.MolWt(mol),
                Descriptors.MolLogP(mol),
                Descriptors.TPSA(mol),
                Descriptors.NumHDonors(mol),
                Descriptors.NumHAcceptors(mol),
                Descriptors.NumRotatableBonds(mol),
            ],
            dtype=np.float32,
        )

        return np.concatenate([fp_arr, basic_descs])

    def predict(
        self,
        smiles_list: List[str],
        endpoints: Optional[List[str]] = None,
        include_shap: bool = False,
    ) -> List[Dict[str, Any]]:
        """
        Predict ADMET properties for compounds.
        """
        endpoints = endpoints or list(ADMET_MODELS.keys())
        results: List[Dict[str, Any]] = []

        for smiles in smiles_list:
            try:
                features = self._get_features(smiles)
            except ImportError as e:
                results.append({"smiles": smiles, "error": str(e)})
                continue
            if features is None:
                results.append({"smiles": smiles, "error": "Invalid SMILES"})
                continue

            pred: Dict[str, Any] = {"smiles": smiles}
            for endpoint in endpoints:
                model_name = ADMET_MODELS.get(endpoint) or _LEGACY_ADMET_MODELS.get(endpoint)
                if not model_name:
                    continue

                try:
                    try:
                        ml_model = self.registry.get_active_model(model_name)
                    except Exception as e:  # noqa: BLE001
                        pred[endpoint] = {"error": f"Model registry unavailable: {type(e).__name__}"}
                        continue
                    if ml_model is None and endpoint in _LEGACY_ADMET_MODELS:
                        # Fallback to legacy model name if ensemble isn't registered yet.
                        try:
                            ml_model = self.registry.get_active_model(_LEGACY_ADMET_MODELS[endpoint])
                        except Exception as e:  # noqa: BLE001
                            pred[endpoint] = {"error": f"Model registry unavailable: {type(e).__name__}"}
                            continue
                    if ml_model is None:
                        pred[endpoint] = {"error": "Model not registered"}
                        continue

                    try:
                        model = self.registry.load_model(ml_model.id)
                    except Exception as e:  # noqa: BLE001
                        pred[endpoint] = {"error": f"Failed to load model: {type(e).__name__}"}
                        continue

                    # Predict
                    X = features.reshape(1, -1)
                    if hasattr(model, "predict_proba"):
                        prob = model.predict_proba(X)[0, 1]
                        pred[endpoint] = {"probability": float(prob), "class": int(prob > 0.5)}
                    else:
                        val = model.predict(X)[0]
                        pred[endpoint] = {"value": float(val)}

                    # SHAP if requested
                    if include_shap:
                        try:
                            import shap  # type: ignore

                            explainer = shap.TreeExplainer(model)
                            shap_values = explainer.shap_values(X)
                            # Top 5 contributing features
                            if isinstance(shap_values, list):
                                sv = shap_values[1][0]  # For binary classifiers
                            else:
                                sv = shap_values[0]
                            top_idx = np.argsort(np.abs(sv))[-5:][::-1]
                            pred[f"{endpoint}_shap"] = {
                                "top_features": top_idx.tolist(),
                                "top_values": sv[top_idx].tolist(),
                            }
                        except ImportError:
                            # SHAP is optional; skip explanations if not installed.
                            logger.info("SHAP not installed; skipping SHAP explanations for %s", endpoint)
                        except Exception as e:  # noqa: BLE001
                            logger.warning("SHAP failed for %s: %s", endpoint, e)

                except Exception as e:  # noqa: BLE001
                    pred[endpoint] = {"error": str(e)}

            results.append(pred)

        return results

    def predict_with_uncertainty(
        self,
        smiles_list: List[str],
        endpoints: Optional[List[str]] = None,
        include_shap: bool = False,
        shap_top_k: int = 10,
    ) -> List[Dict[str, Any]]:
        """Predict ADMET with uncertainty quantification (ensemble + calibration + applicability)."""

        from amprenta_rag.ml.admet.applicability import ApplicabilityChecker

        endpoints = endpoints or list(ADMET_MODELS.keys())
        results: List[Dict[str, Any]] = []

        # Cache loaded artifacts per endpoint for efficiency.
        artifact_cache: Dict[str, Dict[str, Any]] = {}

        def _load_artifact(endpoint: str) -> Optional[Dict[str, Any]]:
            if endpoint in artifact_cache:
                return artifact_cache[endpoint]

            model_name = ADMET_MODELS.get(endpoint)
            if not model_name:
                return None
            try:
                ml_model = self.registry.get_active_model(model_name)
            except Exception:
                # Registry/DB not available in some environments (e.g., E2E runs without Postgres).
                return None
            if ml_model is None:
                # No ensemble registered -> not available for uncertainty method.
                return None
            try:
                obj = self.registry.load_model(ml_model.id)
            except Exception:
                return None
            if not isinstance(obj, dict) or "ensemble" not in obj:
                # Not an ensemble artifact (unexpected); treat as unavailable.
                return None
            artifact_cache[endpoint] = obj
            return obj

        for smiles in smiles_list:
            try:
                features = self._get_features(smiles)
            except ImportError as e:
                results.append({"smiles": smiles, "predictions": {}, "error": str(e)})
                continue
            if features is None:
                results.append({"smiles": smiles, "predictions": {}, "error": "Invalid SMILES"})
                continue

            X = features.reshape(1, -1)
            out: Dict[str, Any] = {"smiles": smiles, "predictions": {}, "error": None}

            for endpoint in endpoints:
                artifact = _load_artifact(endpoint)
                if artifact is None:
                    # Keep schema-compatible shape: mean required; use NaN when ensemble not available.
                    out["predictions"][endpoint] = {
                        "mean": float("nan"),
                        "std": None,
                        "ci_low": None,
                        "ci_high": None,
                        "in_domain": None,
                        "similarity": None,
                        "calibrated": False,
                    }
                    continue

                ens_art = artifact.get("ensemble") or {}
                models = ens_art.get("models") or []
                task = ens_art.get("task") or ("classification" if endpoint in CLASSIFICATION_ENDPOINTS else "regression")

                # Mean/std across ensemble members
                try:
                    if task == "classification":
                        probs = []
                        for m in models:
                            p = m.predict_proba(X)[:, 1]
                            probs.append(np.asarray(p, dtype=float))
                        mat = np.vstack(probs)
                        mean = float(np.mean(mat))
                        std = float(np.std(mat))
                    else:
                        preds = []
                        for m in models:
                            p = m.predict(X)
                            preds.append(np.asarray(p, dtype=float))
                        mat = np.vstack(preds)
                        mean = float(np.mean(mat))
                        std = float(np.std(mat))
                except Exception as e:  # noqa: BLE001
                    out["predictions"][endpoint] = {"calibrated": False, "error": str(e)}
                    continue

                calibrated = False
                cal_obj = artifact.get("calibrator")
                if task == "classification" and cal_obj is not None and hasattr(cal_obj, "calibrate"):
                    try:
                        mean = float(cal_obj.calibrate(np.asarray([mean], dtype=float))[0])
                        calibrated = True
                    except Exception:  # noqa: BLE001
                        calibrated = False

                # Applicability
                app_cfg = artifact.get("applicability") or {}
                threshold = float(app_cfg.get("threshold", 0.3))
                centroid = app_cfg.get("centroid")
                if centroid is None:
                    centroid = ens_art.get("training_centroid")

                in_domain = None
                sim = None
                ci_low = mean - 1.96 * std
                ci_high = mean + 1.96 * std

                if centroid is not None:
                    try:
                        app = ApplicabilityChecker(threshold=threshold)
                        app.training_centroid = np.asarray(centroid, dtype=float)
                        in_domain_arr, sim_arr = app.check(X)
                        in_domain = bool(in_domain_arr[0])
                        sim = float(sim_arr[0])

                        # widen CI for OOD
                        ci_low_a, ci_high_a = app.widen_ci(
                            np.asarray([ci_low], dtype=float),
                            np.asarray([ci_high], dtype=float),
                            np.asarray([in_domain], dtype=bool),
                        )
                        ci_low = float(ci_low_a[0])
                        ci_high = float(ci_high_a[0])
                    except Exception:  # noqa: BLE001
                        in_domain = None
                        sim = None

                out["predictions"][endpoint] = {
                    "mean": float(mean),
                    "std": float(std),
                    "ci_low": float(ci_low),
                    "ci_high": float(ci_high),
                    "in_domain": in_domain,
                    "similarity": sim,
                    "calibrated": bool(calibrated),
                }

                if include_shap:
                    try:
                        from amprenta_rag.ml.admet.explainer import EnsembleSHAPExplainer

                        expl = EnsembleSHAPExplainer(artifact)
                        shap_out = expl.explain_prediction(X, top_k=int(shap_top_k))
                        out["predictions"][endpoint]["shap"] = {
                            "top_features": shap_out.get("top_features") or [],
                            "other_sum": float(shap_out.get("other_sum", 0.0)),
                            "base_value": float(shap_out.get("base_value", 0.0)),
                        }
                    except ImportError:
                        out["predictions"][endpoint]["shap"] = {"error": "shap not installed"}
                    except Exception as e:  # noqa: BLE001
                        out["predictions"][endpoint]["shap"] = {"error": str(e)}

            results.append(out)

        return results


_predictor: Optional[ADMETPredictor] = None


def get_admet_predictor() -> ADMETPredictor:
    global _predictor
    if _predictor is None:
        _predictor = ADMETPredictor()
    return _predictor


__all__ = ["ADMETPredictor", "get_admet_predictor", "ADMET_MODELS"]


