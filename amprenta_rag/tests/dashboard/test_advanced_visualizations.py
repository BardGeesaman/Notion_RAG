from __future__ import annotations

import importlib.util
import sys
import types
from contextlib import contextmanager
from pathlib import Path

import pytest


FILES = [
    Path("scripts/dashboard/components/cytoscape.py"),
    Path("scripts/dashboard/pages/visualizations/cytoscape_network.py"),
    Path("scripts/dashboard/pages/visualizations/entity_explorer.py"),
    Path("scripts/dashboard/pages/visualizations/data_grid.py"),
    Path("scripts/dashboard/pages/visualizations/genome_browser.py"),
]


def _load_module_from_path(module_name: str, path: Path):
    spec = importlib.util.spec_from_file_location(module_name, str(path))
    assert spec and spec.loader
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)  # type: ignore[attr-defined]
    return mod


def _install_streamlit_stubs(monkeypatch):
    """Stub streamlit + components for import-only tests."""
    st_mod = types.ModuleType("streamlit")
    st_mod.header = lambda *a, **k: None
    st_mod.caption = lambda *a, **k: None
    st_mod.markdown = lambda *a, **k: None
    st_mod.info = lambda *a, **k: None
    st_mod.selectbox = lambda *a, **k: (a[1][0] if len(a) > 1 and a[1] else None)
    st_mod.slider = lambda *a, **k: (k.get("value") if "value" in k else None)
    st_mod.text_input = lambda *a, **k: (k.get("value") if "value" in k else "")
    st_mod.button = lambda *a, **k: False
    st_mod.sidebar = types.SimpleNamespace(markdown=lambda *a, **k: None, info=lambda *a, **k: None)

    comps_mod = types.ModuleType("streamlit.components")
    comps_v1_mod = types.ModuleType("streamlit.components.v1")
    comps_v1_mod.html = lambda *a, **k: None

    monkeypatch.setitem(sys.modules, "streamlit", st_mod)
    monkeypatch.setitem(sys.modules, "streamlit.components", comps_mod)
    monkeypatch.setitem(sys.modules, "streamlit.components.v1", comps_v1_mod)


def _install_aggrid_stub(monkeypatch):
    """Stub st_aggrid to validate import path."""
    ag_mod = types.ModuleType("st_aggrid")

    class _GridOptionsBuilder:
        @staticmethod
        def from_dataframe(_df):
            return _GridOptionsBuilder()

        def configure_pagination(self, *a, **k):
            return None

        def configure_default_column(self, *a, **k):
            return None

        def configure_column(self, *a, **k):
            return None

        def build(self):
            return {}

    class _GridUpdateMode:
        VALUE_CHANGED = "VALUE_CHANGED"

    ag_mod.AgGrid = lambda *a, **k: {"data": None}
    ag_mod.GridOptionsBuilder = _GridOptionsBuilder
    ag_mod.GridUpdateMode = _GridUpdateMode
    monkeypatch.setitem(sys.modules, "st_aggrid", ag_mod)


def _install_scripts_db_session_stub(monkeypatch):
    """Stub scripts.dashboard.db_session.db_session used by visualization pages."""
    # Ensure package-ish hierarchy exists
    scripts_mod = sys.modules.get("scripts") or types.ModuleType("scripts")
    dash_mod = sys.modules.get("scripts.dashboard") or types.ModuleType("scripts.dashboard")
    comps_mod = sys.modules.get("scripts.dashboard.components") or types.ModuleType("scripts.dashboard.components")

    @contextmanager
    def _db_session():
        # yields a db object with minimal .query shape
        yield types.SimpleNamespace(query=lambda *a, **k: types.SimpleNamespace(
            filter=lambda *aa, **kk: types.SimpleNamespace(order_by=lambda *aaa, **kkk: types.SimpleNamespace(all=lambda: [])),
            order_by=lambda *aa, **kk: types.SimpleNamespace(all=lambda: []),
        ))

    dbs_mod = types.ModuleType("scripts.dashboard.db_session")
    dbs_mod.db_session = _db_session

    monkeypatch.setitem(sys.modules, "scripts", scripts_mod)
    monkeypatch.setitem(sys.modules, "scripts.dashboard", dash_mod)
    monkeypatch.setitem(sys.modules, "scripts.dashboard.components", comps_mod)
    monkeypatch.setitem(sys.modules, "scripts.dashboard.db_session", dbs_mod)


def _install_cytoscape_component_stub(monkeypatch):
    """Ensure imports like from scripts.dashboard.components.cytoscape import render_cytoscape succeed."""
    # Package-ish hierarchy
    scripts_mod = sys.modules.get("scripts") or types.ModuleType("scripts")
    dash_mod = sys.modules.get("scripts.dashboard") or types.ModuleType("scripts.dashboard")
    comps_mod = sys.modules.get("scripts.dashboard.components") or types.ModuleType("scripts.dashboard.components")
    cyt_mod = types.ModuleType("scripts.dashboard.components.cytoscape")
    cyt_mod.render_cytoscape = lambda *a, **k: None

    monkeypatch.setitem(sys.modules, "scripts", scripts_mod)
    monkeypatch.setitem(sys.modules, "scripts.dashboard", dash_mod)
    monkeypatch.setitem(sys.modules, "scripts.dashboard.components", comps_mod)
    monkeypatch.setitem(sys.modules, "scripts.dashboard.components.cytoscape", cyt_mod)


def _install_amprenta_models_stub(monkeypatch):
    """Stub amprenta_rag.database.models for import-only checks."""
    models_mod = types.ModuleType("amprenta_rag.database.models")

    class _M:
        id = "id"
        name = "name"
        created_at = types.SimpleNamespace(desc=lambda: None)

    models_mod.Program = _M
    models_mod.Experiment = _M
    models_mod.Dataset = _M
    models_mod.Signature = _M
    models_mod.Feature = _M

    # Ensure package hierarchy exists
    db_pkg = sys.modules.get("amprenta_rag.database") or types.ModuleType("amprenta_rag.database")
    monkeypatch.setitem(sys.modules, "amprenta_rag.database", db_pkg)
    monkeypatch.setitem(sys.modules, "amprenta_rag.database.models", models_mod)


@pytest.mark.unit
def test_no_merge_artifacts_or_duplicate_future_imports():
    for rel in FILES:
        path = Path(__file__).resolve().parents[3] / rel
        text = path.read_text(encoding="utf-8")
        assert "<<<<<<<" not in text
        assert ">>>>>>>" not in text
        assert "=======" not in text
        # Common accidental duplication signal:
        assert text.count("from __future__ import annotations") <= 1


@pytest.mark.unit
def test_genome_browser_uses_json_dumps_for_locus():
    path = Path(__file__).resolve().parents[3] / "scripts/dashboard/pages/visualizations/genome_browser.py"
    text = path.read_text(encoding="utf-8")
    assert "json.dumps(locus)" in text


@pytest.mark.unit
def test_aggrid_imports_streamlit_aggrid_module_name():
    path = Path(__file__).resolve().parents[3] / "scripts/dashboard/pages/visualizations/data_grid.py"
    text = path.read_text(encoding="utf-8")
    assert "from st_aggrid import" in text


@pytest.mark.unit
def test_visualization_modules_import_without_syntaxerror(monkeypatch):
    _install_streamlit_stubs(monkeypatch)
    _install_aggrid_stub(monkeypatch)
    _install_scripts_db_session_stub(monkeypatch)
    _install_cytoscape_component_stub(monkeypatch)
    _install_amprenta_models_stub(monkeypatch)

    for rel in FILES:
        path = Path(__file__).resolve().parents[3] / rel
        # Load directly; any SyntaxError will fail the test.
        _load_module_from_path(f"viz_{rel.stem}", path)


