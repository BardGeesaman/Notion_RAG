import builtins
import importlib.util
import logging
import types
from contextlib import contextmanager
from pathlib import Path
from types import SimpleNamespace
from unittest.mock import Mock

import pytest


def _load_module_from_path(module_name: str, path: Path):
    spec = importlib.util.spec_from_file_location(module_name, str(path))
    assert spec and spec.loader
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)  # type: ignore[attr-defined]
    return module


def _install_fake_rdkit_for_sar(monkeypatch, *, mol_from_smiles_impl=None):
    """
    Install a minimal fake RDKit module tree that satisfies imports used by:
    - amprenta_rag.api.services.sar_data
    - deploy/jupyterhub/templates/notebook_utils.py
    """
    rdkit_mod = types.ModuleType("rdkit")
    chem_mod = types.ModuleType("rdkit.Chem")
    allchem_mod = types.ModuleType("rdkit.Chem.AllChem")
    datastructs_mod = types.ModuleType("rdkit.DataStructs")

    class _Flags:
        SANITIZE_ALL = 0b11111111
        SANITIZE_KEKULIZE = 0b00000010

    def default_mol_from_smiles(smiles: str, **kwargs):
        # Default: always "parses"
        return {"smiles": smiles, "kwargs": kwargs}

    chem_mod.MolFromSmiles = mol_from_smiles_impl or default_mol_from_smiles
    chem_mod.MolFromSmarts = lambda smarts: {"smarts": smarts}
    chem_mod.SanitizeMol = Mock()
    chem_mod.SanitizeFlags = _Flags

    allchem_mod.GetMorganFingerprintAsBitVect = Mock(side_effect=lambda m, r, nBits=2048: ("fp", m, r, nBits))
    datastructs_mod.TanimotoSimilarity = Mock(return_value=0.75)

    rdkit_mod.Chem = chem_mod
    rdkit_mod.DataStructs = datastructs_mod
    rdlogger_mod = types.ModuleType("rdkit.RDLogger")
    rdlogger_mod.DisableLog = Mock()
    rdkit_mod.RDLogger = rdlogger_mod

    chem_mod.AllChem = allchem_mod

    monkeypatch.setitem(__import__("sys").modules, "rdkit", rdkit_mod)
    monkeypatch.setitem(__import__("sys").modules, "rdkit.Chem", chem_mod)
    monkeypatch.setitem(__import__("sys").modules, "rdkit.Chem.AllChem", allchem_mod)
    monkeypatch.setitem(__import__("sys").modules, "rdkit.DataStructs", datastructs_mod)
    monkeypatch.setitem(__import__("sys").modules, "rdkit.RDLogger", rdlogger_mod)

    return SimpleNamespace(
        rdkit=rdkit_mod,
        Chem=chem_mod,
        AllChem=allchem_mod,
        DataStructs=datastructs_mod,
        RDLogger=rdlogger_mod,
    )


@pytest.mark.unit
def test_mol_from_smiles_best_effort_empty_returns_none():
    from amprenta_rag.api.services import sar_data

    assert sar_data._mol_from_smiles_best_effort("") is None
    assert sar_data._mol_from_smiles_best_effort(None) is None  # type: ignore[arg-type]


@pytest.mark.unit
def test_mol_from_smiles_best_effort_rdkit_missing_returns_none(monkeypatch):
    from amprenta_rag.api.services import sar_data

    orig_import = builtins.__import__

    def _import(name, *args, **kwargs):
        if name.startswith("rdkit"):
            raise ImportError("rdkit not installed")
        return orig_import(name, *args, **kwargs)

    monkeypatch.setattr(builtins, "__import__", _import)
    assert sar_data._mol_from_smiles_best_effort("CCO") is None


@pytest.mark.unit
def test_mol_from_smiles_best_effort_fallback_sanitizes_without_kekulize(monkeypatch):
    from amprenta_rag.api.services import sar_data

    calls = {"n": 0}

    def mol_from_smiles(smiles: str, **kwargs):
        calls["n"] += 1
        # First attempt fails, fallback returns a mol-like object
        if calls["n"] == 1 and kwargs.get("sanitize") is None:
            return None
        return {"mol": smiles, "kwargs": kwargs}

    fake = _install_fake_rdkit_for_sar(monkeypatch, mol_from_smiles_impl=mol_from_smiles)
    m = sar_data._mol_from_smiles_best_effort("c1ccccc1")
    assert m is not None
    # Ensure fallback sanitize happened
    assert fake.Chem.SanitizeMol.called


@pytest.mark.unit
def test_list_targets_transforms_rows_and_coerces_counts(monkeypatch):
    from amprenta_rag.api.services import sar_data

    db = Mock()
    q = Mock()
    db.query.return_value = q
    q.filter.return_value = q
    q.group_by.return_value = q
    q.order_by.return_value = q
    q.limit.return_value = q
    q.all.return_value = [("CDK2", 3), ("EGFR", None)]

    @contextmanager
    def _fake_db_session():
        yield db

    monkeypatch.setattr(sar_data, "db_session", _fake_db_session)
    out = sar_data.list_targets(limit=123)
    assert out == [
        {"target": "CDK2", "compound_count": 3},
        {"target": "EGFR", "compound_count": 0},
    ]
    q.limit.assert_called_once_with(123)


@pytest.mark.unit
def test_get_compounds_by_target_filters_missing_essentials(monkeypatch):
    from amprenta_rag.api.services import sar_data

    br1 = SimpleNamespace(ic50=10.0, units="nM", assay_name="A", result_id="R1")
    c1 = SimpleNamespace(compound_id="C1", smiles="CCO")
    br2 = SimpleNamespace(ic50=5.0, units="nM", assay_name="A", result_id="R2")
    c2 = SimpleNamespace(compound_id=None, smiles="CCN")  # filtered (missing compound_id)
    br3 = SimpleNamespace(ic50=1.0, units="nM", assay_name="A", result_id="R3")
    c3 = SimpleNamespace(compound_id="C3", smiles=None)  # filtered (missing smiles)

    db = Mock()
    q = Mock()
    db.query.return_value = q
    q.join.return_value = q
    q.filter.return_value = q
    q.order_by.return_value = q
    q.limit.return_value = q
    q.all.return_value = [(br1, c1), (br2, c2), (br3, c3)]

    @contextmanager
    def _fake_db_session():
        yield db

    monkeypatch.setattr(sar_data, "db_session", _fake_db_session)
    out = sar_data.get_compounds_by_target("CDK2", limit=10)
    assert out == [
        {
            "compound_id": "C1",
            "smiles": "CCO",
            "ic50": 10.0,
            "units": "nM",
            "assay_name": "A",
            "result_id": "R1",
        }
    ]


@pytest.mark.unit
def test_get_activity_cliffs_returns_empty_when_rdkit_missing(monkeypatch):
    from amprenta_rag.api.services import sar_data

    orig_import = builtins.__import__

    def _import(name, *args, **kwargs):
        if name.startswith("rdkit"):
            raise ImportError("rdkit not installed")
        return orig_import(name, *args, **kwargs)

    monkeypatch.setattr(builtins, "__import__", _import)
    assert sar_data.get_activity_cliffs_for_target("CDK2") == []


@pytest.mark.unit
def test_get_activity_cliffs_detects_fold_change_and_sorts(monkeypatch):
    from amprenta_rag.api.services import sar_data

    _install_fake_rdkit_for_sar(monkeypatch)

    # Input rows include invalid/non-numeric/<=0 activity and invalid smiles; those should be dropped
    monkeypatch.setattr(
        sar_data,
        "get_compounds_by_target",
        Mock(
            return_value=[
                {"compound_id": "C1", "smiles": "S1", "ic50": 100.0},
                {"compound_id": "C2", "smiles": "S2", "ic50": 5.0},
                {"compound_id": "C3", "smiles": "S3", "ic50": 40.0},  # below fold-change threshold vs C2
                {"compound_id": "C4", "smiles": "BAD", "ic50": 1.0},  # invalid mol
                {"compound_id": "C5", "smiles": "S5", "ic50": 0},  # invalid activity
                {"compound_id": "C6", "smiles": "S6", "ic50": "nope"},  # invalid activity
            ]
        ),
    )
    monkeypatch.setattr(
        sar_data,
        "_mol_from_smiles_best_effort",
        lambda smi: None if smi == "BAD" else {"mol": smi},
    )

    cliffs = sar_data.get_activity_cliffs_for_target(
        "CDK2",
        similarity_threshold=0.6,
        fold_change=10.0,
        limit=50,
    )
    assert cliffs == [
        {
            "compound_1": "C1",
            "smiles_1": "S1",
            "activity_1": 100.0,
            "compound_2": "C2",
            "smiles_2": "S2",
            "activity_2": 5.0,
            "similarity": 0.75,
            "fold_change": 20.0,
            "assay_id": None,
        }
    ]


@pytest.mark.unit
def test_notebook_utils_suppress_rdkit_warnings_sets_levels_and_disables(monkeypatch):
    notebook_utils = _load_module_from_path(
        "notebook_utils",
        Path(__file__).resolve().parents[3] / "deploy/jupyterhub/templates/notebook_utils.py",
    )

    fake = _install_fake_rdkit_for_sar(monkeypatch)
    lg = Mock()
    monkeypatch.setattr(logging, "getLogger", Mock(return_value=lg))

    notebook_utils.suppress_rdkit_warnings()
    fake.RDLogger.DisableLog.assert_called_once_with("rdApp.*")
    lg.setLevel.assert_called_once_with(logging.WARNING)


@pytest.mark.unit
def test_notebook_utils_get_api_client_uses_provided_then_env_then_defaults(monkeypatch):
    notebook_utils = _load_module_from_path(
        "notebook_utils",
        Path(__file__).resolve().parents[3] / "deploy/jupyterhub/templates/notebook_utils.py",
    )

    attempted = []

    class _Http:
        def __init__(self, ok: bool):
            self._ok = ok

        def get(self, path: str):
            if not self._ok:
                raise RuntimeError("health check failed")
            return {"status": "healthy"}

    class FakeRAGClient:
        def __init__(self, api_url: str):
            attempted.append(api_url)
            # Only env URL works in this test
            self.http = _Http(ok=(api_url == "http://env-ok:8000"))

    import os
    import amprenta_rag.client as client_mod

    monkeypatch.setattr(client_mod, "RAGClient", FakeRAGClient)
    monkeypatch.setenv("AMPRENTA_API_URL", "http://env-ok:8000")

    client, is_demo = notebook_utils.get_api_client(api_url="http://provided-bad:8000")
    assert is_demo is False
    assert client is not None
    assert attempted[:2] == ["http://provided-bad:8000", "http://env-ok:8000"]


@pytest.mark.unit
def test_notebook_utils_get_api_client_demo_mode_when_all_fail(monkeypatch):
    notebook_utils = _load_module_from_path(
        "notebook_utils",
        Path(__file__).resolve().parents[3] / "deploy/jupyterhub/templates/notebook_utils.py",
    )

    class _Http:
        def get(self, path: str):
            raise RuntimeError("nope")

    class FakeRAGClient:
        def __init__(self, api_url: str):
            self.http = _Http()

    import amprenta_rag.client as client_mod

    monkeypatch.setattr(client_mod, "RAGClient", FakeRAGClient)
    client, is_demo = notebook_utils.get_api_client(api_url="http://also-bad:8000")
    assert client is None
    assert is_demo is True


@pytest.mark.unit
def test_notebook_utils_mol_from_smiles_safe_fallback_sanitizes(monkeypatch):
    notebook_utils = _load_module_from_path(
        "notebook_utils",
        Path(__file__).resolve().parents[3] / "deploy/jupyterhub/templates/notebook_utils.py",
    )

    calls = {"n": 0}

    def mol_from_smiles(smiles: str, **kwargs):
        calls["n"] += 1
        if calls["n"] == 1 and kwargs.get("sanitize") is None:
            return None
        return {"mol": smiles, "kwargs": kwargs}

    fake = _install_fake_rdkit_for_sar(monkeypatch, mol_from_smiles_impl=mol_from_smiles)
    m = notebook_utils.mol_from_smiles_safe("c1ccccc1")
    assert m is not None
    assert fake.Chem.SanitizeMol.called


