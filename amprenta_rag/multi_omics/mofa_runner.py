"""MOFA runner wrapper (mofapy2)."""

from __future__ import annotations

from pathlib import Path
from typing import Dict

import pandas as pd


def run_mofa(
    data_dict: Dict[str, pd.DataFrame],
    n_factors: int = 10,
    convergence_mode: str = "fast",
    output_dir: str = "data/multi_omics",
) -> str:
    """Run MOFA via `mofapy2` and save model to HDF5.

    Args:
        data_dict: Dict of omics_type -> DataFrame (features x samples).
        n_factors: Number of latent factors.
        convergence_mode: Training mode, e.g. "fast".
        output_dir: Output folder; writes {output_dir}/mofa_model.hdf5

    Returns:
        Path to the written HDF5 model file.

    Raises:
        ImportError if mofapy2 is not installed.
    """
    create_mofa_model = None
    entry_point = None
    try:
        # Preferred API (as requested by Architect spec)
        from mofapy2.run.entry_point import create_mofa_model as _create_mofa_model  # type: ignore

        create_mofa_model = _create_mofa_model
    except Exception:
        create_mofa_model = None

    if create_mofa_model is None:
        try:
            # Backward-compatible API
            from mofapy2.run.entry_point import entry_point as _entry_point  # type: ignore

            entry_point = _entry_point
        except Exception as e:  # noqa: BLE001
            raise ImportError("mofapy2 is not installed. Install mofapy2 to run MOFA.") from e

    if not data_dict:
        raise ValueError("data_dict is empty")

    # mofapy2 expects samples x features per view; we store features x samples internally.
    views = {k: v.T.copy() for k, v in data_dict.items()}
    for k, df in views.items():
        if df.empty:
            raise ValueError(f"View '{k}' is empty")

    out_dir = Path(output_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    out_path = out_dir / "mofa_model.hdf5"

    if create_mofa_model is not None:
        model = create_mofa_model()
        model.set_data_df(views)
        model.set_model_options(factors=int(n_factors))
        model.set_train_options(convergence_mode=str(convergence_mode))
        model.build()
        model.run()
        model.save(str(out_path))
    else:
        if entry_point is None:
            raise ImportError("mofapy2 entry_point is not available")
        ent = entry_point()
        ent.set_data_df(views)
        ent.set_model_options(factors=int(n_factors))
        ent.set_train_options(convergence_mode=str(convergence_mode))
        ent.build()
        ent.run()
        ent.save(str(out_path))

    return str(out_path)


__all__ = ["run_mofa"]


