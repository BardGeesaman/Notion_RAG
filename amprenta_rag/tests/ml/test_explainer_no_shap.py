from __future__ import annotations

import builtins

import pytest

from amprenta_rag.ml.admet.explainer import EnsembleSHAPExplainer


def test_shap_not_installed_error(monkeypatch):
    real_import = builtins.__import__

    def fake_import(name, globals=None, locals=None, fromlist=(), level=0):  # noqa: ANN001
        if name == "shap":
            raise ImportError("no shap")
        return real_import(name, globals, locals, fromlist, level)

    monkeypatch.setattr(builtins, "__import__", fake_import)

    with pytest.raises(ImportError) as e:
        EnsembleSHAPExplainer({"ensemble": {"models": [object()]}})

    assert "pip install shap" in str(e.value)


