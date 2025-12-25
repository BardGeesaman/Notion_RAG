from __future__ import annotations

import json
from pathlib import Path
from typing import List, Tuple

import numpy as np

from amprenta_rag.ml.admet.predictor import ADMETPredictor


def _require_pandas():
    try:
        import pandas as pd  # type: ignore
    except Exception as e:  # noqa: BLE001
        raise ImportError("pandas is required for dataset loading (pip install pandas)") from e
    return pd


class ChEMBLDatasetLoader:
    def __init__(self, data_dir: str = "data/admet"):
        self.data_dir = Path(data_dir)

    def load_endpoint(self, endpoint: str, split: str = "train") -> Tuple[np.ndarray, np.ndarray, List[str]]:
        """Load endpoint parquet and featurize SMILES.

        Returns:
          (X, y, smiles)
        - Invalid or missing SMILES are skipped.
        - Rows whose features cannot be computed are skipped.
        """

        pd = _require_pandas()
        endpoint = (endpoint or "").strip().lower()
        split = (split or "").strip().lower()
        if not endpoint:
            raise ValueError("endpoint is required")
        if split not in ("train", "cal", "test"):
            raise ValueError("split must be one of: train, cal, test")

        path = self.data_dir / f"{endpoint}_{split}.parquet"
        if not path.exists():
            raise FileNotFoundError(f"Parquet not found: {path}")

        df = pd.read_parquet(path)
        if df.empty:
            return np.zeros((0, 0), dtype=np.float32), np.zeros((0,), dtype=np.float32), []

        if "canonical_smiles" not in df.columns or "y" not in df.columns:
            raise ValueError("Parquet must include columns: canonical_smiles, y")

        predictor = ADMETPredictor()

        feats: List[np.ndarray] = []
        ys: List[float] = []
        smiles_out: List[str] = []

        for _, row in df.iterrows():
            smi = row.get("canonical_smiles")
            if not isinstance(smi, str) or not smi.strip():
                continue
            f = predictor._get_features(smi)  # noqa: SLF001
            if f is None:
                continue
            try:
                yv = float(row.get("y"))
            except Exception:  # noqa: BLE001
                continue
            feats.append(np.asarray(f, dtype=np.float32))
            ys.append(yv)
            smiles_out.append(smi)

        if not feats:
            return np.zeros((0, 0), dtype=np.float32), np.zeros((0,), dtype=np.float32), []

        X = np.vstack([x.reshape(1, -1) for x in feats]).astype(np.float32)
        y = np.asarray(ys, dtype=np.float32)
        return X, y, smiles_out

    def get_summary(self) -> dict:
        path = self.data_dir / "summary.json"
        if not path.exists():
            raise FileNotFoundError(f"summary.json not found: {path}")
        return json.loads(path.read_text())

    def list_endpoints(self) -> List[str]:
        """List endpoints based on parquet files found in data_dir."""
        if not self.data_dir.exists():
            return []
        endpoints = set()
        for p in self.data_dir.glob("*.parquet"):
            name = p.name
            # expected: {endpoint}_{split}.parquet
            if not name.endswith(".parquet") or "_" not in name:
                continue
            base = name[: -len(".parquet")]
            parts = base.split("_")
            if len(parts) < 2:
                continue
            split = parts[-1]
            endpoint = "_".join(parts[:-1])
            if split in ("train", "cal", "test") and endpoint:
                endpoints.add(endpoint)
        return sorted(endpoints)


__all__ = ["ChEMBLDatasetLoader"]



