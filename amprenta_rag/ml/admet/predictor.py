"""ADMET Predictor service."""

from __future__ import annotations

from typing import Any, Dict, List, Optional

import numpy as np

from amprenta_rag.logging_utils import get_logger
from amprenta_rag.ml.registry import get_registry

logger = get_logger(__name__)

# ADMET model names
ADMET_MODELS = {
    "herg": "admet_herg_xgb",
    "logs": "admet_logs_xgb",
    "logp": "admet_logp_xgb",
}


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
            features = self._get_features(smiles)
            if features is None:
                results.append({"smiles": smiles, "error": "Invalid SMILES"})
                continue

            pred: Dict[str, Any] = {"smiles": smiles}
            for endpoint in endpoints:
                model_name = ADMET_MODELS.get(endpoint)
                if not model_name:
                    continue

                try:
                    ml_model = self.registry.get_active_model(model_name)
                    if ml_model is None:
                        pred[endpoint] = {"error": "Model not registered"}
                        continue

                    model = self.registry.load_model(ml_model.id)

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


_predictor: Optional[ADMETPredictor] = None


def get_admet_predictor() -> ADMETPredictor:
    global _predictor
    if _predictor is None:
        _predictor = ADMETPredictor()
    return _predictor


__all__ = ["ADMETPredictor", "get_admet_predictor", "ADMET_MODELS"]


