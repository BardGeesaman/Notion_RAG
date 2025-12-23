from __future__ import annotations

from contextlib import contextmanager
from types import SimpleNamespace
from uuid import uuid4

import pytest


@pytest.mark.unit
def test_sphingolipid_imbalance_enzyme_mode(monkeypatch):
    from amprenta_rag.analysis.sphingolipid import scoring as sc

    ds_ids = [uuid4()]

    # Fake DB session returns features containing enzymes
    fake_feats = [
        SimpleNamespace(name="CERS2", normalized_name=None),
        SimpleNamespace(name="UGCG", normalized_name=None),
        SimpleNamespace(name="SMPD1", normalized_name=None),
    ]

    fake_db = SimpleNamespace(
        query=lambda model: SimpleNamespace(
            join=lambda *args, **kwargs: SimpleNamespace(
                filter=lambda *a, **k: SimpleNamespace(all=lambda: fake_feats)
            )
        )
    )

    @contextmanager
    def fake_db_session():
        yield fake_db

    monkeypatch.setattr(sc, "db_session", fake_db_session)
    out = sc.compute_pathway_imbalance(ds_ids, pathway="ceramide")
    assert out["method"] == "enzymes"
    assert out["score"] > 0
    assert "CERS2" in out["matched"]["enzymes"]


@pytest.mark.unit
def test_sphingolipid_imbalance_ratio_mode(monkeypatch):
    from amprenta_rag.analysis.sphingolipid import scoring as sc

    ds_ids = [uuid4()]

    fake_feats = [
        SimpleNamespace(name="Cer 16:0", normalized_name=None),
        SimpleNamespace(name="SM 24:1", normalized_name=None),
        SimpleNamespace(name="hex_cer_24_0", normalized_name=None),
    ]

    fake_db = SimpleNamespace(
        query=lambda model: SimpleNamespace(
            join=lambda *args, **kwargs: SimpleNamespace(
                filter=lambda *a, **k: SimpleNamespace(all=lambda: fake_feats)
            )
        )
    )

    @contextmanager
    def fake_db_session():
        yield fake_db

    monkeypatch.setattr(sc, "db_session", fake_db_session)
    out = sc.compute_pathway_imbalance(ds_ids, pathway="ceramide")
    assert out["method"] == "ratios"
    assert "lipid_class_counts" in out["matched"]
    counts = out["matched"]["lipid_class_counts"]
    assert counts["CER"] >= 1
    assert counts["SM"] >= 1
    assert counts["HEXCER"] >= 1


