"""SHAP explainability utilities for ADMET ensemble models."""

from __future__ import annotations

from typing import Any, Dict, List

import numpy as np

from amprenta_rag.ml.admet.features import get_feature_name


class EnsembleSHAPExplainer:
    """Aggregate SHAP values across an ensemble of tree models."""

    def __init__(self, ensemble_artifact: Dict[str, Any]):
        try:
            import shap  # type: ignore
        except Exception as e:  # noqa: BLE001
            raise ImportError(
                "shap is required for EnsembleSHAPExplainer. Install with: pip install shap"
            ) from e

        self._shap = shap
        self._artifact = dict(ensemble_artifact or {})

        try:
            self.models = list(self._artifact["ensemble"]["models"])
        except Exception as e:  # noqa: BLE001
            raise ValueError("Invalid ensemble_artifact: expected artifact['ensemble']['models']") from e

        if not self.models:
            raise ValueError("Invalid ensemble_artifact: models list is empty")

    def _shap_for_model(self, model: Any, X: np.ndarray) -> tuple[np.ndarray, float]:
        expl = self._shap.TreeExplainer(model)
        sv = expl.shap_values(X)

        # shap can return list for multiclass; keep first output for MVP
        if isinstance(sv, list):
            sv_arr = np.asarray(sv[0])
        else:
            sv_arr = np.asarray(sv)

        base = getattr(expl, "expected_value", 0.0)
        if isinstance(base, (list, tuple, np.ndarray)):
            base_val = float(np.asarray(base).ravel()[0])
        else:
            base_val = float(base)

        return sv_arr, base_val

    def explain_prediction(self, X: np.ndarray, top_k: int = 10) -> Dict[str, Any]:
        """Explain a single prediction.

        Args:
            X: shape (1, 2054)
            top_k: number of top features by |SHAP| to return
        """
        X = np.asarray(X)
        if X.shape != (1, 2054):
            raise ValueError(f"X must have shape (1, 2054), got {X.shape}")
        top_k = int(top_k)
        if top_k <= 0:
            raise ValueError("top_k must be > 0")

        shap_vals: List[np.ndarray] = []
        base_vals: List[float] = []
        for m in self.models:
            sv, base = self._shap_for_model(m, X)
            sv = np.asarray(sv).reshape(1, -1)
            if sv.shape[1] != 2054:
                raise ValueError(f"Model SHAP output must have 2054 features, got {sv.shape}")
            shap_vals.append(sv[0])
            base_vals.append(float(base))

        mean_sv = np.mean(np.stack(shap_vals, axis=0), axis=0)  # (2054,)
        base_value = float(np.mean(np.asarray(base_vals)))

        k = min(top_k, mean_sv.shape[0])
        idx = np.argsort(np.abs(mean_sv))[::-1]
        top_idx = idx[:k]
        top_set = set(int(i) for i in top_idx.tolist())

        top_features: List[Dict[str, Any]] = []
        for rank, i in enumerate(top_idx.tolist(), start=1):
            i_int = int(i)
            top_features.append(
                {"name": get_feature_name(i_int), "value": float(mean_sv[i_int]), "rank": int(rank)}
            )

        other_sum = float(np.sum([mean_sv[i] for i in range(mean_sv.shape[0]) if i not in top_set]))

        return {
            "shap_values": mean_sv,
            "top_features": top_features,
            "base_value": base_value,
            "other_sum": other_sum,
        }

    def compute_global_importance(self, X: np.ndarray) -> List[Dict[str, Any]]:
        """Compute mean |SHAP| importance across all samples and ensemble models."""
        X = np.asarray(X)
        if X.ndim != 2 or X.shape[1] != 2054:
            raise ValueError(f"X must have shape (N, 2054), got {X.shape}")

        per_model: List[np.ndarray] = []
        for m in self.models:
            sv, _base = self._shap_for_model(m, X)
            sv_arr = np.asarray(sv)
            if sv_arr.ndim == 1:
                sv_arr = sv_arr.reshape(1, -1)
            if sv_arr.shape[1] != 2054:
                raise ValueError(f"Model SHAP output must have 2054 features, got {sv_arr.shape}")
            per_model.append(np.mean(np.abs(sv_arr), axis=0))

        mean_imp = np.mean(np.stack(per_model, axis=0), axis=0)  # (2054,)
        order = np.argsort(mean_imp)[::-1]

        out: List[Dict[str, Any]] = []
        for rank, i in enumerate(order.tolist(), start=1):
            i_int = int(i)
            out.append({"name": get_feature_name(i_int), "importance": float(mean_imp[i_int]), "rank": int(rank)})
        return out


__all__ = ["EnsembleSHAPExplainer"]


