from __future__ import annotations

from uuid import uuid4

from amprenta_rag.ingestion.features import postgres_linking as pl


class FakeColumn:
    def __init__(self, attr: str):
        self.attr = attr

    def __eq__(self, other):
        return lambda obj: getattr(obj, self.attr) == other

    def ilike(self, other):
        return lambda obj: getattr(obj, self.attr, "").lower() == other.lower()

    def in_(self, values):
        values_set = set(values)
        return lambda obj: getattr(obj, self.attr) in values_set


class FakeFeature:
    name = FakeColumn("name")
    feature_type = FakeColumn("feature_type")

    def __init__(self, name: str, feature_type: str):
        self.name = name
        self.feature_type = feature_type
        self.id = uuid4()


class FakeDataset:
    id = FakeColumn("id")

    def __init__(self, id_):
        self.id = id_
        self.features = []


class FakeItem:
    id = FakeColumn("id")

    def __init__(self, id_, semantic_metadata=None):
        self.id = id_
        self.semantic_metadata = semantic_metadata or {}


class FakeQuery:
    def __init__(self, data):
        self._data = list(data)

    def filter(self, *conds):
        def matches(obj):
            return all(cond(obj) for cond in conds)

        self._data = [obj for obj in self._data if matches(obj)]
        return self

    def first(self):
        return self._data[0] if self._data else None

    def all(self):
        return list(self._data)


class FakeSession:
    def __init__(self, datasets=None, features=None, literature=None, emails=None):
        self.datasets = datasets or []
        self.features = features or []
        self.literature = literature or []
        self.emails = emails or []
        self.commits = 0
        self.rollbacks = 0

    def query(self, model):
        if model is pl.DatasetModel:
            return FakeQuery(self.datasets)
        if model is pl.FeatureModel:
            return FakeQuery(self.features)
        if model is pl.Literature:
            return FakeQuery(self.literature)
        if model is pl.Email:
            return FakeQuery(self.emails)
        return FakeQuery([])

    def add(self, obj):
        if isinstance(obj, FakeFeature):
            self.features.append(obj)

    def commit(self):
        self.commits += 1

    def rollback(self):
        self.rollbacks += 1

    def refresh(self, obj):
        return obj

    def close(self):
        pass


def test_link_features_to_dataset_uses_append(monkeypatch):
    ds_id = uuid4()
    dataset = FakeDataset(ds_id)
    session = FakeSession(datasets=[dataset])

    monkeypatch.setattr(pl, "_find_or_create_feature_impl", lambda name, t, db: FakeFeature(name, t))

    pl.link_features_to_postgres_items(["A", "B"], ds_id, "dataset", db=session)

    assert len(dataset.features) == 2


def test_link_features_to_literature_updates_metadata(monkeypatch):
    lit_id = uuid4()
    literature = FakeItem(lit_id, semantic_metadata={"feature_mentions": ["old"]})
    session = FakeSession(literature=[literature])

    monkeypatch.setattr(pl, "_find_or_create_feature_impl", lambda name, t, db: FakeFeature(name, t))

    pl.link_features_to_postgres_items(["X", "Y"], lit_id, "literature", db=session)

    assert "feature_mentions" in literature.semantic_metadata
    assert set(literature.semantic_metadata["feature_mentions"]) == {"old", "x", "y"}


def test_find_or_create_feature_uses_external_db(monkeypatch):
    session = FakeSession(features=[FakeFeature("a", "metabolite")])

    monkeypatch.setattr(pl, "_find_or_create_feature_impl", lambda name, t, db: "created")
    monkeypatch.setattr(pl, "get_db", lambda: iter([session]))

    result = pl.find_or_create_feature_in_postgres("A", "metabolite")
    assert result == "created"


def test_batch_link_features_creates_and_links(monkeypatch):
    ds_id = uuid4()
    existing = FakeFeature("existing", "metabolite")
    dataset = FakeDataset(ds_id)
    session = FakeSession(datasets=[dataset], features=[existing])

    # Monkeypatch model classes to our fake columns
    monkeypatch.setattr(pl, "FeatureModel", FakeFeature)
    monkeypatch.setattr(pl, "DatasetModel", FakeDataset)
    monkeypatch.setattr(pl, "get_db", lambda: iter([session]))

    features = [("existing", "metabolite"), ("new", "metabolite")]
    results = pl.batch_link_features_to_dataset_in_postgres(features, ds_id, db=session, max_workers=1)

    assert results["existing"] is True
    assert results["new"] is True
    # Two features linked to dataset
    assert len(dataset.features) == 2


def test_batch_link_returns_false_when_dataset_missing(monkeypatch):
    session = FakeSession(datasets=[])
    monkeypatch.setattr(pl, "get_db", lambda: iter([session]))

    features = [("f1", "metabolite")]
    results = pl.batch_link_features_to_dataset_in_postgres(features, uuid4(), db=session)
    assert results == {"f1": False}

