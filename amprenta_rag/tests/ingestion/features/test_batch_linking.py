from __future__ import annotations

from amprenta_rag.ingestion.features import batch_linking as bl


def test_batch_find_or_create_feature_pages_uses_cache(monkeypatch):
    bl.clear_feature_cache()
    bl._feature_page_cache[("gene", "TP53")] = "page1"
    result = bl.batch_find_or_create_feature_pages([("gene", "TP53")], max_workers=1)
    assert result[("gene", "TP53")] == "page1"


def test_batch_find_or_create_feature_pages_calls_underlying(monkeypatch):
    bl.clear_feature_cache()

    def fake_find(ftype, fname):
        return f"{ftype}:{fname}"

    monkeypatch.setattr(bl, "_feature_page_cache", {})
    monkeypatch.setattr(bl, "_relation_property_cache", {})
    monkeypatch.setattr(
        bl,
        "ThreadPoolExecutor",
        lambda max_workers=1: __import__("contextlib").nullcontext(
            type("E", (), {"submit": lambda self, fn, *a: type("F", (), {"result": lambda self2: fn(*a)})()})()
        ),
    )
    monkeypatch.setattr(
        bl,
        "as_completed",
        lambda futures: futures,
    )
    monkeypatch.setattr(bl, "_feature_page_cache", {})
    monkeypatch.setattr(
        bl,
        "ThreadPoolExecutor",
        lambda max_workers=1: __import__("contextlib").nullcontext(
            type("E", (), {"submit": lambda self, fn, *a: type("F", (), {"result": lambda self2: fn(*a)})()})()
        ),
    )
    monkeypatch.setattr(bl, "as_completed", lambda futures: futures)
    monkeypatch.setattr(bl, "_feature_page_cache", {})
    monkeypatch.setattr(
        __import__("amprenta_rag.ingestion.features.general_linking", fromlist=["_find_or_create_feature_page"]),
        "_find_or_create_feature_page",
        fake_find,
    )

    res = bl.batch_find_or_create_feature_pages([("gene", "TP53")], max_workers=1)
    assert res[("gene", "TP53")] == "gene:TP53"


def test_batch_link_features_skips_when_disabled(monkeypatch):
    bl.clear_feature_cache()

    monkeypatch.setattr(bl, "batch_find_or_create_feature_pages", lambda feats, max_workers=1: {("g", "t"): "p1"})
    monkeypatch.setattr(bl, "batch_add_dataset_relations", lambda *a, **k: (_ for _ in ()).throw(AssertionError("should not be called")))

    res = bl.batch_link_features([("g", "t")], "ds1", enable_linking=False)
    assert res[("g", "t")] == "p1"


def test_clear_feature_cache_resets():
    bl._feature_page_cache[("x", "y")] = "p"
    bl._relation_property_cache["p"] = "rel"
    bl.clear_feature_cache()
    assert bl._feature_page_cache == {}
    assert bl._relation_property_cache == {}

