"""Seed placeholder ADMET models into the ML registry.

This is intended for local/dev environments to exercise the registry + ADMET
prediction API without running full training.

Prefers XGBoost models if available; otherwise falls back to sklearn dummy models.
"""

from __future__ import annotations

from typing import Any, Dict, Tuple


def _make_models() -> Tuple[Any, Any, Any, Dict[str, float]]:
    import numpy as np

    X = np.random.default_rng(0).normal(size=(50, 2054))

    # classification labels
    y_cls = (np.random.default_rng(1).random(size=(50,)) > 0.5).astype(int)
    # regression labels
    y_reg = np.random.default_rng(2).normal(size=(50,))

    try:
        from xgboost import XGBClassifier, XGBRegressor  # type: ignore

        m_herg = XGBClassifier(n_estimators=20, max_depth=3, learning_rate=0.2)
        m_logs = XGBRegressor(n_estimators=20, max_depth=3, learning_rate=0.2)
        m_logp = XGBRegressor(n_estimators=20, max_depth=3, learning_rate=0.2)

        m_herg.fit(X, y_cls)
        m_logs.fit(X, y_reg)
        m_logp.fit(X, y_reg)
        return m_herg, m_logs, m_logp, {"seed": 1.0}
    except Exception:
        from sklearn.dummy import DummyClassifier, DummyRegressor  # type: ignore

        m_herg = DummyClassifier(strategy="most_frequent")
        m_logs = DummyRegressor(strategy="mean")
        m_logp = DummyRegressor(strategy="mean")

        m_herg.fit(X, y_cls)
        m_logs.fit(X, y_reg)
        m_logp.fit(X, y_reg)
        return m_herg, m_logs, m_logp, {"seed": 1.0}


def main() -> None:
    from amprenta_rag.ml.registry import get_registry

    registry = get_registry()
    m_herg, m_logs, m_logp, metrics = _make_models()

    registry.register_model(
        name="admet_herg_xgb",
        version="v0.0.0-seed",
        model_type="admet_classification",
        framework="xgboost",
        model_object=m_herg,
        metrics=metrics,
        description="Seed placeholder model (dev-only)",
    )
    registry.register_model(
        name="admet_logs_xgb",
        version="v0.0.0-seed",
        model_type="admet_regression",
        framework="xgboost",
        model_object=m_logs,
        metrics=metrics,
        description="Seed placeholder model (dev-only)",
    )
    registry.register_model(
        name="admet_logp_xgb",
        version="v0.0.0-seed",
        model_type="admet_regression",
        framework="xgboost",
        model_object=m_logp,
        metrics=metrics,
        description="Seed placeholder model (dev-only)",
    )

    print("Seeded 3 placeholder ADMET models.")


if __name__ == "__main__":
    main()


