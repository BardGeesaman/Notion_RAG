from __future__ import annotations

import numpy as np
import pytest

xgb = pytest.importorskip("xgboost")
sklearn = pytest.importorskip("sklearn")
_ = (xgb, sklearn)


def test_isotonic_calibration_monotonic():
    from amprenta_rag.ml.admet.calibration import CalibrationWrapper

    y_pred = np.linspace(0.0, 1.0, 200)
    y_true = (y_pred > 0.5).astype(int)

    cal = CalibrationWrapper(method="isotonic").fit(y_pred, y_true)
    y_cal = cal.calibrate(y_pred)

    # monotonic (non-decreasing) for sorted inputs
    assert np.all(np.diff(y_cal) >= -1e-9)
    assert np.all((y_cal >= 0.0) & (y_cal <= 1.0))


def test_platt_calibration_sigmoid():
    from amprenta_rag.ml.admet.calibration import CalibrationWrapper

    rng = np.random.default_rng(0)
    x = rng.normal(size=(500,))
    y_pred = 1.0 / (1.0 + np.exp(-x))
    y_true = (x > 0.0).astype(int)

    cal = CalibrationWrapper(method="platt").fit(y_pred, y_true)
    y_cal = cal.calibrate(y_pred)

    # should be bounded and roughly monotonic vs y_pred
    order = np.argsort(y_pred)
    assert np.all(np.diff(y_cal[order]) >= -1e-6)
    assert np.all((y_cal >= 0.0) & (y_cal <= 1.0))


def test_ece_perfect_calibration():
    from amprenta_rag.ml.admet.calibration import CalibrationWrapper

    y_pred = np.array([0.0, 1.0, 0.0, 1.0], dtype=float)
    y_true = np.array([0.0, 1.0, 0.0, 1.0], dtype=float)

    cal = CalibrationWrapper(method="isotonic")
    ece = cal.compute_ece(y_pred, y_true, n_bins=2)
    assert ece == 0.0


def test_reliability_diagram_structure():
    from amprenta_rag.ml.admet.calibration import reliability_diagram

    y_pred = np.array([0.1, 0.9, 0.2, 0.8], dtype=float)
    y_true = np.array([0.0, 1.0, 0.0, 1.0], dtype=float)

    out = reliability_diagram(y_pred, y_true, n_bins=4)
    assert set(out.keys()) == {"bin_edges", "bin_accs", "bin_confs", "ece"}
    assert len(out["bin_edges"]) == 5
    assert len(out["bin_accs"]) == 4
    assert len(out["bin_confs"]) == 4
    assert isinstance(out["ece"], float)


