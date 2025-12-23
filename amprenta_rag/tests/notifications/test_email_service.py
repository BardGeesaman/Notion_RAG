from __future__ import annotations

from uuid import uuid4
from amprenta_rag.notifications import email_service as es

# Define FakeDB at module level so it's available to test_send_share_notification_unknown
class FakeQuery:
    def __init__(self, obj=None):
        self.obj = obj

    def filter(self, *a, **k):
        return self

    def first(self):
        return self.obj

class FakeDB:
    def __init__(self, obj=None):
        self.obj = obj

    def query(self, model):
        return FakeQuery(self.obj)

    def __enter__(self):
        return self

    def __exit__(self, *exc):
        return False

class _FakeSMTP:
    def __init__(self):
        self.started_tls = False
        self.logged_in = False
        self.sent = []

    def starttls(self):
        self.started_tls = True

    def login(self, user, password):
        self.logged_in = True
        self.user = user
        self.password = password

    def send_message(self, msg):
        self.sent.append(msg)

    def __enter__(self):
        return self

    def __exit__(self, exc_type, exc, tb):
        return False


def test_is_email_configured(monkeypatch):
    monkeypatch.setattr(es, "SMTP_USER", "u")
    monkeypatch.setattr(es, "SMTP_PASSWORD", "p")
    assert es.is_email_configured() is True
    monkeypatch.setattr(es, "SMTP_PASSWORD", "")
    assert es.is_email_configured() is False


def test_send_email_not_configured(monkeypatch, caplog):
    monkeypatch.setattr(es, "SMTP_USER", None)
    monkeypatch.setattr(es, "SMTP_PASSWORD", None)
    assert es.send_email("a@b.com", "subj", "body") is False
    assert any("not configured" in rec.message for rec in caplog.records)


def test_send_email_success_with_attachment(monkeypatch):
    smtp = _FakeSMTP()
    monkeypatch.setattr(es, "SMTP_USER", "u")
    monkeypatch.setattr(es, "SMTP_PASSWORD", "p")
    monkeypatch.setattr(es, "FROM_EMAIL", "from@example.com")
    monkeypatch.setattr(es, "smtplib", type("S", (), {"SMTP": lambda *a, **k: smtp}))

    ok = es.send_email(
        to="to@example.com",
        subject="hello",
        body="plain",
        html_body="<b>hi</b>",
        attachments=[("file.txt", b"data", "text/plain")],
    )
    assert ok is True
    assert smtp.started_tls and smtp.logged_in
    assert smtp.sent and smtp.sent[0]["To"] == "to@example.com"


def test_send_experiment_summary_missing(monkeypatch):
    monkeypatch.setattr(es, "send_email", lambda *a, **k: True)
    # Inject db_session into module globals for the service to pick up
    monkeypatch.setitem(es.__dict__, "db_session", lambda: FakeDB())
    assert es.send_experiment_summary(uuid4(), "to@example.com", None) is False


def test_send_experiment_summary_happy(monkeypatch):
    calls = {}
    class FakeExperiment:
        def __init__(self):
            self.id = uuid4()
            self.name = "Exp"
            self.description = "desc"
            self.design_type = "design"
            self.organism = "human"
            self.datasets = []

    fake_exp = FakeExperiment()
    monkeypatch.setitem(es.__dict__, "db_session", lambda: FakeDB(fake_exp))
    monkeypatch.setattr(es, "get_experiment_summary_html", lambda exp, ds: "html")
    monkeypatch.setattr(es, "send_email", lambda to, subj, body, html_body=None: calls.setdefault("sent", True))

    assert es.send_experiment_summary(fake_exp.id, "to@example.com", None) is True
    assert calls.get("sent") is True


def test_send_share_notification_compound(monkeypatch):
    calls = {}
    class FakeCompound:
        def __init__(self):
            self.id = uuid4()
            self.compound_id = "CID"
            self.name = "Name"

    monkeypatch.setitem(es.__dict__, "db_session", lambda: FakeDB(FakeCompound()))
    monkeypatch.setattr(es, "get_share_email_html", lambda *a, **k: "html")
    monkeypatch.setattr(es, "send_email", lambda to, subj, body, html_body=None: calls.setdefault("sent", True))

    assert es.send_share_notification("compound", uuid4(), "to@example.com", "user", "msg", None) is True
    assert calls.get("sent") is True


def test_send_share_notification_unknown(monkeypatch):
    monkeypatch.setitem(es.__dict__, "db_session", lambda: FakeDB())
    monkeypatch.setattr(es, "get_share_email_html", lambda *a, **k: "html")
    monkeypatch.setattr(es, "send_email", lambda *a, **k: True)
    assert es.send_share_notification("unknown", uuid4(), "to@example.com", "user", None, None) is True
