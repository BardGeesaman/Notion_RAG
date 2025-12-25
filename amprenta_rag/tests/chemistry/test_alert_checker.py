from __future__ import annotations

import pytest

pytest.importorskip("rdkit")

from amprenta_rag.chemistry.alert_checker import StructuralAlertChecker, get_traffic_light  # noqa: E402
from amprenta_rag.chemistry.alerts import AlertResult  # noqa: E402


def test_check_compound_returns_structure() -> None:
    chk = StructuralAlertChecker()
    out = chk.check_compound("CCO")
    assert out["smiles"] == "CCO"
    assert "alerts" in out and isinstance(out["alerts"], list)
    assert "summary" in out and isinstance(out["summary"], dict)
    assert "traffic_light" in out


def test_traffic_light_red_on_pains() -> None:
    alerts = [AlertResult("PAINS", "x", "x", "high", "")]
    assert get_traffic_light(alerts) == "RED"


def test_traffic_light_yellow_on_brenk_only() -> None:
    alerts = [AlertResult("BRENK", "x", "x", "medium", "")]
    assert get_traffic_light(alerts) == "YELLOW"


def test_traffic_light_green_on_clean() -> None:
    assert get_traffic_light([]) == "GREEN"


def test_check_batch_parallel() -> None:
    chk = StructuralAlertChecker()
    out = chk.check_batch(["CCO", "CC(=O)Cl", "O=[N+]([O-])c1ccccc1"], n_jobs=2)
    assert len(out) == 3
    assert all("smiles" in r for r in out)


def test_invalid_smiles_handled() -> None:
    chk = StructuralAlertChecker()
    out = chk.check_compound("NOT_A_SMILES")
    assert out.get("error") == "Invalid SMILES"
    assert out["alerts"] == []
    assert out["traffic_light"] == "RED"



