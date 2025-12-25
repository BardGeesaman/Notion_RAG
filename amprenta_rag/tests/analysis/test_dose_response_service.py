from __future__ import annotations

import numpy as np
import pytest


def _synthetic_4pl(conc, *, bottom=0.1, top=1.2, ec50=10.0, hill=1.1, noise=0.03, seed=0):
    rng = np.random.default_rng(seed)
    x = np.asarray(conc, dtype=float)
    y = bottom + (top - bottom) / (1.0 + (x / ec50) ** (-hill))
    y = y + rng.normal(0.0, noise, size=y.shape[0])
    return y


def test_fit_3pl_returns_params():
    from amprenta_rag.analysis.dose_response_service import fit_3pl

    conc = np.logspace(-1, 2, 10)
    resp = _synthetic_4pl(conc, bottom=0.0, top=1.0, ec50=5.0, hill=1.2)
    out = fit_3pl(conc, resp)
    assert set(out.keys()) >= {"ec50", "hill_slope", "top", "r_squared"}
    assert out["ec50"] > 0


def test_fit_4pl_returns_params():
    from amprenta_rag.analysis.dose_response_service import fit_4pl

    conc = np.logspace(-1, 2, 12)
    resp = _synthetic_4pl(conc, bottom=0.05, top=1.1, ec50=12.0, hill=1.0)
    out = fit_4pl(conc, resp)
    assert set(out.keys()) >= {"ec50", "hill_slope", "top", "bottom", "r_squared"}
    assert out["ec50"] > 0


def test_fit_dose_response_dispatcher():
    from amprenta_rag.analysis.dose_response_service import fit_dose_response

    conc = list(map(float, np.logspace(-1, 2, 12)))
    resp = list(map(float, _synthetic_4pl(conc, bottom=0.1, top=1.0, ec50=8.0, hill=1.3)))

    out4 = fit_dose_response(conc, resp, model="4PL")
    assert out4.model == "4PL"
    assert float(out4.params["ec50"]) > 0

    out3 = fit_dose_response(conc, resp, model="3PL")
    assert out3.model == "3PL"


def test_fit_dose_response_dispatcher_bayesian_optional():
    pytest.importorskip("pymc")
    pytest.importorskip("arviz")

    from amprenta_rag.analysis.dose_response_service import fit_dose_response

    conc = list(map(float, np.logspace(-1, 2, 12)))
    resp = list(map(float, _synthetic_4pl(conc, bottom=0.1, top=1.0, ec50=8.0, hill=1.3)))

    bout = fit_dose_response(conc, resp, model="bayesian_4pl")
    assert bout.model == "bayesian_4pl"
    assert "ec50" in bout.params


def test_generate_fit_curve_creates_points():
    from amprenta_rag.analysis.dose_response_service import generate_fit_curve

    params = {"bottom": 0.1, "top": 1.0, "ec50": 10.0, "hill_slope": 1.0}
    xs, ys = generate_fit_curve(params, x_min=0.1, x_max=100.0, n_points=50)
    assert len(xs) == 50
    assert len(ys) == 50
    assert all(x > 0 for x in xs)


def test_bootstrap_ci_returns_intervals():
    from amprenta_rag.analysis.dose_response_service import bootstrap_ci, fit_4pl

    conc = list(map(float, np.logspace(-1, 2, 12)))
    resp = list(map(float, _synthetic_4pl(conc, bottom=0.1, top=1.0, ec50=10.0, hill=1.2, noise=0.02)))

    ci = bootstrap_ci(conc, resp, fit_4pl, n_bootstrap=50)
    assert set(ci.keys()) >= {"ec50_ci", "hill_ci", "r_squared_ci"}
    lo, hi = ci["ec50_ci"]
    assert lo < hi


def test_compare_curves_returns_table():
    from amprenta_rag.analysis.dose_response_service import compare_curves

    out = compare_curves(
        [
            {"label": "A", "ec50": 10.0, "hill_slope": 1.0, "r_squared": 0.9},
            {"label": "B", "ec50": 20.0, "hill_slope": 1.1, "r_squared": 0.85},
        ]
    )
    assert "parameter_table" in out
    assert "ec50_fold_changes" in out
    assert len(out["parameter_table"]) == 2
    assert out["ec50_fold_changes"][0] == pytest.approx(1.0)
    assert out["ec50_fold_changes"][1] == pytest.approx(2.0)


