"""MOFA HDF5 output parser (tolerant to mofapy2 version differences)."""

from __future__ import annotations

from typing import Any, Dict, Iterator, Tuple

import pandas as pd


def parse_mofa_output(hdf5_path: str) -> Dict[str, Any]:
    """Parse a MOFA model HDF5 file into factors, loadings, scores, and variance metadata.

    Returns:
        {
          "factors": [0, 1, ...],
          "loadings": {view: DataFrame(features x factors)},
          "scores": DataFrame(samples x factors),
          "variance": {...}
        }
    """
    try:
        import h5py  # type: ignore
    except Exception as e:  # noqa: BLE001
        raise ImportError("h5py is required to parse MOFA output HDF5 files.") from e

    def _iter_datasets(g, prefix: str = "") -> Iterator[Tuple[str, Any]]:
        for k, v in g.items():
            p = f"{prefix}/{k}" if prefix else str(k)
            if hasattr(v, "shape") and not hasattr(v, "items"):
                yield p, v
            elif hasattr(v, "items"):
                yield from _iter_datasets(v, p)

    def _decode_list(x) -> list[str]:
        if x is None:
            return []
        out: list[str] = []
        for v in list(x):
            if isinstance(v, (bytes, bytearray)):
                out.append(v.decode("utf-8", errors="replace"))
            else:
                out.append(str(v))
        return out

    with h5py.File(hdf5_path, "r") as h5:
        datasets = list(_iter_datasets(h5))

        # Find Z (scores) dataset
        z_path = None
        for p, d in datasets:
            if p.split("/")[-1] == "Z" and getattr(d, "ndim", 0) == 2:
                z_path = p
                break
        if z_path is None:
            # fallback: first 2D dataset named like *factors* or *scores*
            for p, d in datasets:
                if getattr(d, "ndim", 0) == 2 and any(tok in p.lower() for tok in ("score", "factor")):
                    z_path = p
                    break
        if z_path is None:
            raise ValueError("Unable to locate factor scores (Z) matrix in HDF5.")

        Z = h5[z_path][:]
        n_samples, n_factors = int(Z.shape[0]), int(Z.shape[1])

        # Sample names
        sample_names: list[str] = []
        for p, d in datasets:
            name = p.split("/")[-1].lower()
            if name in ("samples", "sample_names") and getattr(d, "ndim", 0) == 1 and int(d.shape[0]) == n_samples:
                sample_names = _decode_list(d[:])
                break
        if not sample_names:
            sample_names = [f"S{i}" for i in range(n_samples)]

        factor_names = list(range(n_factors))
        scores_df = pd.DataFrame(Z, index=sample_names, columns=factor_names)

        # Loadings W per view
        loadings: Dict[str, pd.DataFrame] = {}

        # Common layout: expectations/W/<view>
        for p, d in datasets:
            parts = p.split("/")
            if len(parts) >= 3 and parts[-2] == "W" and getattr(d, "ndim", 0) == 2:
                view = parts[-1]
                W = h5[p][:]
                n_feat = int(W.shape[0])

                feat_names: list[str] = []
                # try to find feature names for this view
                for p2, d2 in datasets:
                    parts2 = p2.split("/")
                    if len(parts2) >= 2 and parts2[-1].lower() in ("features", "feature_names", "feature_ids", "genes"):
                        if view in parts2 and getattr(d2, "ndim", 0) == 1 and int(d2.shape[0]) == n_feat:
                            feat_names = _decode_list(d2[:])
                            break
                if not feat_names:
                    feat_names = [f"F{i}" for i in range(n_feat)]

                loadings[view] = pd.DataFrame(W, index=feat_names, columns=factor_names)

        # Fallback: single dataset named W (no views)
        if not loadings:
            for p, d in datasets:
                if p.split("/")[-1] == "W" and getattr(d, "ndim", 0) == 2:
                    W = h5[p][:]
                    loadings["view1"] = pd.DataFrame(W, index=[f"F{i}" for i in range(W.shape[0])], columns=factor_names)
                    break

        # Variance explained (best-effort)
        variance: Dict[str, Any] = {}
        for p, d in datasets:
            base = p.split("/")[-1].lower()
            if base in ("r2", "variance_explained") or "variance" in p.lower():
                try:
                    arr = h5[p][:]
                    variance[p] = arr.tolist() if hasattr(arr, "tolist") else arr
                except Exception:
                    continue

        return {
            "factors": list(range(n_factors)),
            "loadings": loadings,
            "scores": scores_df,
            "variance": variance,
            "paths": {"Z": z_path},
        }


__all__ = ["parse_mofa_output"]


