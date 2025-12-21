from __future__ import annotations

from datetime import datetime
from uuid import uuid4

from amprenta_rag.api.services import screening


class Field:
    def __init__(self, name: str):
        self.name = name

    def __eq__(self, other):
        return ("eq", self.name, other)

    def desc(self):
        return ("desc", self.name)

    def is_(self, val):
        return ("is_true" if val is True else "is_false", self.name)


class FakeCampaign:
    run_date = Field("run_date")
    campaign_id = Field("campaign_id")
    id = Field("id")

    def __init__(self, **kwargs):
        self.campaign_id = kwargs.get("campaign_id", "cid")
        self.campaign_name = kwargs.get("campaign_name")
        self.description = kwargs.get("description")
        self.assay_type = kwargs.get("assay_type")
        self.target = kwargs.get("target")
        self.library_id = kwargs.get("library_id")
        self.total_wells = kwargs.get("total_wells")
        self.hit_count = kwargs.get("hit_count")
        self.run_date = kwargs.get("run_date", datetime.utcnow())
        self.id = kwargs.get("id", uuid4())


class FakeResult:
    campaign_id_field = Field("campaign_id")
    hit_flag_field = Field("hit_flag")
    normalized_value_field = Field("normalized_value")
    campaign_id = campaign_id_field
    hit_flag = hit_flag_field
    normalized_value = normalized_value_field

    def __init__(self, **kwargs):
        self.id = kwargs.get("id", uuid4())
        self.result_id = kwargs.get("result_id")
        self.compound_id = kwargs.get("compound_id")
        self.campaign_id = kwargs.get("campaign_id")
        self.hit_flag = kwargs.get("hit_flag")
        self.normalized_value = kwargs.get("normalized_value")
        self.well_position = kwargs.get("well_position")
        self.raw_value = kwargs.get("raw_value")
        self.z_score = kwargs.get("z_score")
        self.hit_category = kwargs.get("hit_category")


class FakeQuery:
    def __init__(self, data):
        self._data = list(data)

    def order_by(self, *args):
        if args and isinstance(args[0], tuple) and args[0][0] == "desc":
            field = args[0][1]
            self._data.sort(key=lambda x: getattr(x, field, None), reverse=True)
        return self

    def filter(self, *predicates):
        filtered = []
        for obj in self._data:
            ok = True
            for pred in predicates:
                if isinstance(pred, tuple):
                    if len(pred) == 3 and pred[0] == "eq":
                        field, val = pred[1], pred[2]
                        if getattr(obj, field, None) != val:
                            ok = False
                            break
                    elif len(pred) == 2 and pred[0] in {"is_true", "is_false"}:
                        field = pred[1]
                        desired = pred[0] == "is_true"
                        if bool(getattr(obj, field, None)) != desired:
                            ok = False
                            break
            if ok:
                filtered.append(obj)
        return FakeQuery(filtered)

    def join(self, *args, **kwargs):
        return self

    def all(self):
        return list(self._data)

    def first(self):
        return self._data[0] if self._data else None


class FakeSession:
    def __init__(self, campaigns=None, results=None):
        self.campaigns = campaigns or []
        self.results = results or []

    def query(self, model):
        if model is FakeCampaign:
            return FakeQuery(self.campaigns)
        if model is FakeResult:
            return FakeQuery(self.results)
        return FakeQuery([])

    def __enter__(self):
        return self

    def __exit__(self, exc_type, exc, tb):
        return False


def test_list_campaigns(monkeypatch):
    c1 = FakeCampaign(campaign_id="C1")
    c2 = FakeCampaign(campaign_id="C2")
    monkeypatch.setattr(screening, "db_session", lambda: FakeSession([c1, c2], []))
    monkeypatch.setattr(screening, "HTSCampaign", FakeCampaign)
    rows = screening.list_campaigns()
    assert len(rows) == 2
    assert rows[0]["campaign_id"] in {"C1", "C2"}


def test_get_campaign(monkeypatch):
    c1 = FakeCampaign(campaign_id="C1")
    monkeypatch.setattr(screening, "db_session", lambda: FakeSession([c1], []))
    monkeypatch.setattr(screening, "HTSCampaign", FakeCampaign)
    row = screening.get_campaign("C1")
    assert row["campaign_id"] == "C1"
    assert screening.get_campaign("missing") is None


def test_get_campaign_hits(monkeypatch):
    c1 = FakeCampaign(id="C1", campaign_id="C1")
    r1 = FakeResult(campaign_id="C1", hit_flag=True, normalized_value=2.0, compound_id="cmp")
    r2 = FakeResult(campaign_id="C1", hit_flag=False, normalized_value=1.0)
    monkeypatch.setattr(screening, "db_session", lambda: FakeSession([c1], [r1, r2]))
    monkeypatch.setattr(screening, "HTSCampaign", FakeCampaign)
    monkeypatch.setattr(screening, "HTSResult", FakeResult)

    hits = screening.get_campaign_hits("C1")
    assert len(hits) == 1
    assert hits[0]["compound_id"] == "cmp"

