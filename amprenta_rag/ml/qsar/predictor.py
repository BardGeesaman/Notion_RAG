"""QSAR predictor for per-target ensemble models."""

from __future__ import annotations

from dataclasses import dataclass
from typing import Any, Dict, List, Optional, Tuple

import numpy as np

from amprenta_rag.logging_utils import get_logger
from amprenta_rag.ml.admet.applicability import ApplicabilityChecker
from amprenta_rag.ml.admet.predictor import ADMETPredictor
from amprenta_rag.ml.registry import get_registry


logger = get_logger(__name__)


def _ci_95(mean: np.ndarray, std: np.ndarray) -> Tuple[np.ndarray, np.ndarray]:
    mean = np.asarray(mean, dtype=float).reshape(-1)
    std = np.asarray(std, dtype=float).reshape(-1)
    z = 1.96
    lo = mean - z * std
    hi = mean + z * std
    return lo, hi


def _predict_ensemble_proba(ensemble_artifact: Dict[str, Any], X: np.ndarray) -> Tuple[np.ndarray, np.ndarray]:
    """Compute mean/std probability from a serialized ensemble artifact."""
    models = ensemble_artifact.get("models") or []
    task = (ensemble_artifact.get("task") or "classification").lower()
    if task != "classification":
        raise ValueError("QSAR predictor expects classification ensembles")
    if not models:
        raise ValueError("No models in ensemble artifact")

    X = np.asarray(X)
    probs = []
    for m in models:
        p = m.predict_proba(X)[:, 1]
        probs.append(np.asarray(p, dtype=float))
    mat = np.vstack(probs)
    mean = np.clip(mat.mean(axis=0), 0.0, 1.0)
    std = mat.std(axis=0)
    return mean, std


@dataclass(frozen=True)
class TargetModelInfo:
    target: str
    model_name: str
    version: str
    metrics: Dict[str, Any]

    def to_dict(self) -> Dict[str, Any]:
        return {
            "target": self.target,
            "model_name": self.model_name,
            "version": self.version,
            "metrics": self.metrics,
        }


class TargetQSARPredictor:
    def __init__(self):
        self.registry = get_registry()
        # cache: target -> {"ml_model": MLModel, "artifact": dict}
        self._model_cache: Dict[str, Dict[str, Any]] = {}

        # Avoid ADMETPredictor.__init__ (registry/DB); we only need _get_features.
        self._feature_extractor = object.__new__(ADMETPredictor)

    def list_available_targets(self) -> List[Dict[str, Any]]:
        """List registered QSAR target models (best-effort)."""
        try:
            models = self.registry.list_models()
        except Exception as e:  # noqa: BLE001
            logger.warning("[QSAR] Registry unavailable for list_available_targets: %s", e)
            return []

        out: List[TargetModelInfo] = []
        for m in models:
            name = getattr(m, "name", "") or ""
            if not name.startswith("qsar_") or not name.endswith("_ensemble"):
                continue
            # qsar_{target}_ensemble
            target = name[len("qsar_") : -len("_ensemble")] or name
            out.append(
                TargetModelInfo(
                    target=target,
                    model_name=name,
                    version=str(getattr(m, "version", "") or ""),
                    metrics=getattr(m, "metrics", None) or {},
                )
            )
        return [x.to_dict() for x in out]

    def _load_target_model(self, target: str) -> Optional[Dict[str, Any]]:
        t = str(target)
        if t in self._model_cache:
            return self._model_cache[t]

        model_name = f"qsar_{t}_ensemble"
        try:
            ml_model = self.registry.get_active_model(model_name)
        except Exception as e:  # noqa: BLE001
            logger.warning("[QSAR] Registry unavailable: %s", e)
            return None
        if ml_model is None:
            return None

        try:
            artifact = self.registry.load_model(ml_model.id)
        except Exception as e:  # noqa: BLE001
            logger.warning("[QSAR] Failed to load model %s: %s", model_name, e)
            return None

        cached = {"ml_model": ml_model, "artifact": artifact}
        self._model_cache[t] = cached
        return cached

    def get_model_info(self, target: str) -> Optional[Dict[str, Any]]:
        cached = self._load_target_model(target)
        if not cached:
            return None
        art = cached.get("artifact") or {}
        meta = art.get("metadata") if isinstance(art, dict) else None
        if isinstance(meta, dict):
            return meta
        return None

    def predict(self, smiles_list: List[str], targets: List[str]) -> List[Dict[str, Any]]:
        results: List[Dict[str, Any]] = []
        targets = [str(t) for t in (targets or []) if str(t).strip()]
        if not targets:
            return [{"smiles": str(s or ""), "predictions": {}, "error": "No targets provided"} for s in smiles_list]

        # Preload models (best-effort). If registry unavailable, return empty predictions.
        loaded: Dict[str, Dict[str, Any]] = {}
        for t in targets:
            m = self._load_target_model(t)
            if m:
                loaded[t] = m

        if not loaded:
            return [{"smiles": str(s or ""), "predictions": {}, "error": "No target models available"} for s in smiles_list]

        for smi in smiles_list:
            smi_s = str(smi or "").strip()
            if not smi_s:
                results.append({"smiles": smi_s, "predictions": {}, "error": "Empty SMILES"})
                continue

            try:
                feats = ADMETPredictor._get_features(self._feature_extractor, smi_s)
            except ImportError as e:
                results.append({"smiles": smi_s, "predictions": {}, "error": str(e)})
                continue

            if feats is None:
                results.append({"smiles": smi_s, "predictions": {}, "error": "Invalid SMILES"})
                continue

            X = np.asarray(feats, dtype=np.float32).reshape(1, -1)
            per_target: Dict[str, Any] = {}

            for t, m in loaded.items():
                art = m.get("artifact") or {}
                if not isinstance(art, dict):
                    continue

                try:
                    ens = art.get("ensemble") or {}
                    if not isinstance(ens, dict):
                        continue
                    mean, std = _predict_ensemble_proba(ens, X)
                    prob = mean[0]
                    std0 = float(std[0])

                    calibrated = False
                    calibrator = art.get("calibrator")
                    if calibrator is not None and hasattr(calibrator, "calibrate"):
                        try:
                            prob = float(np.asarray(calibrator.calibrate(np.asarray([prob])))[0])
                            calibrated = True
                        except Exception:  # noqa: BLE001
                            calibrated = False

                    app_cfg = art.get("applicability") or {}
                    checker = ApplicabilityChecker(threshold=float(app_cfg.get("threshold", 0.3)))
                    checker.training_centroid = app_cfg.get("centroid")
                    in_dom, sim = checker.check(X)

                    ci_low, ci_high = _ci_95(np.asarray([prob]), np.asarray([std0]))
                    ci_low2, ci_high2 = checker.widen_ci(ci_low, ci_high, in_dom)

                    per_target[t] = {
                        "probability": float(prob),
                        "std": float(std0),
                        "ci_low": float(ci_low2[0]),
                        "ci_high": float(ci_high2[0]),
                        "in_domain": bool(in_dom[0]),
                        "similarity": float(sim[0]),
                        "calibrated": bool(calibrated),
                        "active": bool(prob > 0.5),
                    }
                except Exception as e:  # noqa: BLE001
                    logger.warning("[QSAR] Prediction failed for %s/%s: %s", smi_s, t, e)
                    continue

            results.append({"smiles": smi_s, "predictions": per_target, "error": None})

        return results


__all__ = ["TargetQSARPredictor"]



