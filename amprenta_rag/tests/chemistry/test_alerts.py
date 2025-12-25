from __future__ import annotations

import pytest

pytest.importorskip("rdkit")

from amprenta_rag.chemistry.alerts import BRENKFilter, LillyFilter, PAINSFilter  # noqa: E402


def test_pains_filter_detects_known_pains() -> None:
    # Rhodanine-like motif (commonly flagged as PAINS)
    smiles = "O=C1NC(=S)SC1=Cc1ccccc1"
    f = PAINSFilter()
    out = f.check(smiles)
    assert len(out) >= 1
    assert all(x.alert_type == "PAINS" for x in out)


def test_brenk_filter_detects_unwanted() -> None:
    # Nitro aromatic (commonly flagged as undesirable)
    smiles = "O=[N+]([O-])c1ccccc1"
    f = BRENKFilter()
    out = f.check(smiles)
    assert len(out) >= 1
    assert all(x.alert_type == "BRENK" for x in out)


def test_lilly_filter_detects_reactive() -> None:
    # Acyl chloride
    smiles = "CC(=O)Cl"
    f = LillyFilter()
    out = f.check(smiles)
    assert any(x.pattern_name == "acyl_halide" for x in out)
    assert all(x.alert_type == "LILLY" for x in out)


def test_clean_compound_no_alerts() -> None:
    # "Clean" compound: ibuprofen (note: aspirin is flagged by RDKit Brenk as phenol_ester)
    smiles = "CC(C)Cc1ccc(cc1)[C@@H](C)C(=O)O"
    pains = PAINSFilter().check(smiles)
    brenk = BRENKFilter().check(smiles)
    lilly = LillyFilter().check(smiles)
    assert pains == []
    assert brenk == []
    assert lilly == []


