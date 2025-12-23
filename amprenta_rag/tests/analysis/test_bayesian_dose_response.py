from __future__ import annotations

import pytest


def test_fit_bayesian_dose_response_synthetic():
    pytest.importorskip("pymc")
    _ = pytest.importorskip("arviz")

    from amprenta_rag.analysis.bayesian_dose_response import fit_bayesian_dose_response

    import numpy as np

    rng = np.random.default_rng(0)
    concentrations = np.logspace(-2, 2, 12)  # 0.01 .. 100

    bottom = 0.1
    top = 1.2
    ec50_true = 5.0
    hill = 1.3

    def four_pl(x):
        return bottom + (top - bottom) / (1.0 + (x / ec50_true) ** (-hill))

    y = four_pl(concentrations) + rng.normal(0, 0.05, size=len(concentrations))

    out = fit_bayesian_dose_response(
        concentrations=list(map(float, concentrations)),
        responses=list(map(float, y)),
        prior_ec50=ec50_true,
        likelihood="normal",
    )

    assert "ec50_mean" in out
    assert "ec50_ci" in out
    assert "hill_slope" in out
    assert "trace" in out
    assert "posterior_predictive" in out

    ec50_mean = float(out["ec50_mean"])
    lo, hi = out["ec50_ci"]
    assert lo < ec50_mean < hi
    # Loose check: within ~10x of truth (Bayesian fit can be noisy with small draws).
    assert ec50_mean > ec50_true / 10
    assert ec50_mean < ec50_true * 10


def test_bayesian_dose_response_diagnostics_api(monkeypatch):
    """Test include_diagnostics returns trace plot."""
    pytest.importorskip("arviz")
    pytest.importorskip("matplotlib")

    import numpy as np
    import arviz as az
    from fastapi.testclient import TestClient

    from amprenta_rag.api.main import app
    from amprenta_rag.api.routers import analysis as analysis_router

    # Fake trace with two chains so R-hat is defined.
    posterior = {
        "ec50": np.random.default_rng(0).lognormal(mean=1.0, sigma=0.2, size=(2, 50)),
        "hill_slope": np.random.default_rng(1).normal(loc=1.0, scale=0.2, size=(2, 50)),
    }
    fake_trace = az.from_dict(posterior=posterior)

    monkeypatch.setattr(
        analysis_router,
        "fit_bayesian_dose_response",
        lambda **kwargs: {
            "ec50_mean": 10.0,
            "ec50_ci": (5.0, 20.0),
            "hill_slope": 1.2,
            "trace": fake_trace,
            "posterior_predictive": None,
        },
    )

    client = TestClient(app)
    resp = client.post(
        "/api/analysis/dose-response/bayesian",
        json={
            "concentrations": [0.1, 1, 10, 100, 1000, 10000],
            "responses": [5, 15, 35, 65, 85, 92],
            "include_diagnostics": True,
        },
    )
    assert resp.status_code == 200
    body = resp.json()
    assert "diagnostics" in body
    diag = body["diagnostics"]
    assert isinstance(diag, dict)
    assert diag["trace_plot"].startswith("data:image/png;base64,")
    assert "rhat" in diag
    assert "ess_bulk" in diag


