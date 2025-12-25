"""QSAR dataset loaders (ChEMBL REST + local DB fallback)."""

from __future__ import annotations

import statistics
from typing import Any, Dict, List, Optional, Tuple

import numpy as np


COMMON_TARGETS: Dict[str, str] = {
    "EGFR": "CHEMBL203",
    "CDK2": "CHEMBL301",
    "BRAF": "CHEMBL5145",
    "JAK2": "CHEMBL2971",
    "ABL1": "CHEMBL1862",
}


def _normalize_to_nm(value: float, unit: str) -> Optional[float]:
    """Normalize to nM. Returns None for unsupported/unknown units."""
    if value is None or unit is None:
        return None
    try:
        v = float(value)
    except Exception:
        return None

    u = str(unit).strip().lower()
    u = u.replace("µ", "μ")  # normalize micro sign

    factors = {
        "nm": 1.0,
        "um": 1_000.0,
        "μm": 1_000.0,
        "mm": 1_000_000.0,
        "m": 1_000_000_000.0,
    }
    f = factors.get(u)
    if f is None:
        return None
    return v * f


def _aggregate_ic50(values: List[float]) -> float:
    """Aggregate replicate IC50 values per compound (median for robustness)."""
    if not values:
        raise ValueError("No values to aggregate")
    return float(statistics.median([float(v) for v in values]))


def _require_requests():
    try:
        import requests  # type: ignore
    except Exception as e:  # noqa: BLE001
        raise ImportError("requests is required for ChEMBL API integration") from e
    return requests


def _chembl_get_json(url: str, params: Dict[str, Any]) -> Dict[str, Any]:
    requests = _require_requests()
    r = requests.get(url, params=params, timeout=30)
    r.raise_for_status()
    out = r.json()
    return out if isinstance(out, dict) else {"raw": out}


def fetch_chembl_target_data(
    target_chembl_id: str,
    threshold_nm: float = 1000,
    limit: int = 5000,
) -> Tuple[np.ndarray, np.ndarray, List[str], Dict[str, Any]]:
    """Fetch IC50 activity data for a ChEMBL target and return features/labels.

    Returns: (X_features, y_labels, smiles_list, metadata)
    """
    base = "https://www.ebi.ac.uk/chembl/api/data"
    activity_url = f"{base}/activity.json"

    collected: Dict[str, List[float]] = {}
    unit_counts: Dict[str, int] = {}
    n_raw = 0

    # Basic pagination (offset/limit); stop on empty page or overall limit.
    page_size = min(1000, int(limit))
    offset = 0
    while len(collected) < int(limit):
        before_n = len(collected)
        params = {
            "target_chembl_id": target_chembl_id,
            "standard_type": "IC50",
            "limit": page_size,
            "offset": offset,
        }
        payload = _chembl_get_json(activity_url, params=params)
        acts = payload.get("activities") or []
        if not isinstance(acts, list) or not acts:
            break

        for a in acts:
            if not isinstance(a, dict):
                continue
            mol_id = a.get("molecule_chembl_id")
            std_val = a.get("standard_value")
            std_unit = a.get("standard_units")
            if not mol_id or std_val is None or not std_unit:
                continue
            nm = _normalize_to_nm(std_val, std_unit)
            if nm is None:
                continue
            unit_key = str(std_unit).strip()
            unit_counts[unit_key] = unit_counts.get(unit_key, 0) + 1
            collected.setdefault(str(mol_id), []).append(float(nm))
            n_raw += 1

            if len(collected) >= int(limit):
                break

        # Defensive: if pagination makes no progress (no new molecule IDs), stop.
        # This prevents infinite loops on repeated pages and makes mocks robust.
        if len(collected) == before_n:
            break

        offset += page_size

    # Fetch SMILES for each molecule.
    smiles_by_mol: Dict[str, str] = {}
    for mol_id in list(collected.keys()):
        mol_url = f"{base}/molecule/{mol_id}.json"
        payload = _chembl_get_json(mol_url, params={})
        ms = payload.get("molecule_structures") or {}
        smi = ms.get("canonical_smiles") if isinstance(ms, dict) else None
        if smi:
            smiles_by_mol[mol_id] = str(smi)

    # Aggregate IC50 per compound and build labels/features.
    from amprenta_rag.ml.admet.predictor import ADMETPredictor

    # Avoid ADMETPredictor.__init__ (which touches the model registry / DB).
    # We only need the feature extractor implementation in _get_features().
    pred = object.__new__(ADMETPredictor)
    features: List[np.ndarray] = []
    labels: List[int] = []
    smiles_list: List[str] = []
    ic50_median_by_mol: Dict[str, float] = {}

    for mol_id, vals in collected.items():
        smi = smiles_by_mol.get(mol_id)
        if not smi:
            continue
        med = _aggregate_ic50(vals)
        ic50_median_by_mol[mol_id] = med
        y = 1 if med < float(threshold_nm) else 0

        x = ADMETPredictor._get_features(pred, smi)
        if x is None:
            continue
        features.append(x.astype(np.float32))
        labels.append(int(y))
        smiles_list.append(smi)

    if features:
        X = np.stack(features, axis=0)
        y_arr = np.asarray(labels, dtype=np.int64)
    else:
        X = np.zeros((0, 2054), dtype=np.float32)
        y_arr = np.zeros((0,), dtype=np.int64)

    active_ratio = float(y_arr.mean()) if y_arr.size else 0.0
    meta: Dict[str, Any] = {
        "source": "chembl",
        "target_chembl_id": target_chembl_id,
        "threshold_nm": float(threshold_nm),
        "n_raw_activities": int(n_raw),
        "n_compounds_total": int(len(collected)),
        "n_compounds_with_smiles": int(len(smiles_by_mol)),
        "n_compounds_used": int(len(smiles_list)),
        "active_count": int(int(y_arr.sum())),
        "active_ratio": float(active_ratio),
        "unit_counts": unit_counts,
        "ic50_median_by_molecule": ic50_median_by_mol,
    }
    return X, y_arr, smiles_list, meta


def load_target_from_db(
    target_name: str,
    threshold_nm: float = 1000,
) -> Tuple[np.ndarray, np.ndarray, List[str], Dict[str, Any]]:
    """Load IC50 data for a target from local DB (BiochemicalResult)."""
    from amprenta_rag.database.session import db_session
    from amprenta_rag.models.chemistry import BiochemicalResult, Compound
    from amprenta_rag.ml.admet.predictor import ADMETPredictor

    pred = object.__new__(ADMETPredictor)
    per_cmp: Dict[str, List[float]] = {}
    smiles_by_cmp: Dict[str, str] = {}
    n_raw = 0

    with db_session() as db:
        rows = (
            db.query(
                BiochemicalResult.compound_id,
                Compound.smiles,
                BiochemicalResult.ic50,
                BiochemicalResult.units,
            )
            .join(Compound, Compound.id == BiochemicalResult.compound_id)
            .filter(BiochemicalResult.target == target_name)
            .filter(BiochemicalResult.ic50.isnot(None))
            .all()
        )

    for compound_id, smiles, ic50, units in rows:
        nm = _normalize_to_nm(ic50, units)
        if nm is None:
            continue
        cid = str(compound_id)
        per_cmp.setdefault(cid, []).append(float(nm))
        if smiles:
            smiles_by_cmp[cid] = str(smiles)
        n_raw += 1

    features: List[np.ndarray] = []
    labels: List[int] = []
    smiles_list: List[str] = []
    ic50_median_by_cmp: Dict[str, float] = {}

    for cid, vals in per_cmp.items():
        smi = smiles_by_cmp.get(cid)
        if not smi:
            continue
        med = _aggregate_ic50(vals)
        ic50_median_by_cmp[cid] = med
        y = 1 if med < float(threshold_nm) else 0
        x = ADMETPredictor._get_features(pred, smi)
        if x is None:
            continue
        features.append(x.astype(np.float32))
        labels.append(int(y))
        smiles_list.append(smi)

    if features:
        X = np.stack(features, axis=0)
        y_arr = np.asarray(labels, dtype=np.int64)
    else:
        X = np.zeros((0, 2054), dtype=np.float32)
        y_arr = np.zeros((0,), dtype=np.int64)

    active_ratio = float(y_arr.mean()) if y_arr.size else 0.0
    meta: Dict[str, Any] = {
        "source": "local",
        "target_name": target_name,
        "threshold_nm": float(threshold_nm),
        "n_raw_rows": int(n_raw),
        "n_compounds_total": int(len(per_cmp)),
        "n_compounds_used": int(len(smiles_list)),
        "active_count": int(int(y_arr.sum())),
        "active_ratio": float(active_ratio),
        "ic50_median_by_compound": ic50_median_by_cmp,
    }
    return X, y_arr, smiles_list, meta


class TargetDatasetLoader:
    """Unified loader for target datasets (ChEMBL or local DB)."""

    def load_target(
        self,
        target: str,
        source: str = "chembl",
        threshold_nm: float = 1000,
        min_compounds: int = 100,
        min_active_ratio: float = 0.2,
    ) -> Tuple[np.ndarray, np.ndarray, List[str], Dict[str, Any]]:
        source = (source or "chembl").strip().lower()
        if source not in {"chembl", "local"}:
            raise ValueError("source must be 'chembl' or 'local'")

        if source == "chembl":
            chembl_id = target
            if not chembl_id.upper().startswith("CHEMBL"):
                chembl_id = COMMON_TARGETS.get(target.upper(), target)
            X, y, smiles, meta = fetch_chembl_target_data(
                chembl_id, threshold_nm=float(threshold_nm), limit=5000
            )
        else:
            X, y, smiles, meta = load_target_from_db(target, threshold_nm=float(threshold_nm))

        n = int(X.shape[0])
        active_ratio = float(y.mean()) if y.size else 0.0

        if n < int(min_compounds):
            raise ValueError(f"Insufficient compounds for QSAR: {n} < {int(min_compounds)}")
        if active_ratio < float(min_active_ratio):
            raise ValueError(
                f"Insufficient active ratio for QSAR: {active_ratio:.3f} < {float(min_active_ratio):.3f}"
            )

        meta = dict(meta or {})
        meta.update(
            {
                "min_compounds": int(min_compounds),
                "min_active_ratio": float(min_active_ratio),
                "compound_count": n,
                "active_ratio": float(active_ratio),
            }
        )
        return X, y, smiles, meta


__all__ = [
    "COMMON_TARGETS",
    "TargetDatasetLoader",
    "fetch_chembl_target_data",
    "load_target_from_db",
    "_normalize_to_nm",
    "_aggregate_ic50",
]


