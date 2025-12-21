from __future__ import annotations

from amprenta_rag.ingestion.features import linking as lk


def test_linking_reexports():
    assert "_find_or_create_feature_page" in lk.__all__
    assert lk.link_feature is not None

