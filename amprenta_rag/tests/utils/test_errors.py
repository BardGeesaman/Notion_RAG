from __future__ import annotations

from amprenta_rag.utils import errors


def test_format_error_known_category_includes_suggestions():
    msg = errors.format_error("db_connection", details="timeout")
    assert "Database connection failed" in msg
    assert "Suggestions" in msg


def test_format_error_unknown_category():
    msg = errors.format_error("unknown", details="x")
    assert msg.startswith("[unknown]")


def test_render_cli_error(monkeypatch):
    captured = []

    class FakeLogger:
        def error(self, text):
            captured.append(text)

    monkeypatch.setattr(errors, "logger", FakeLogger())
    errors.render_cli_error("network_error", details="down")
    assert captured, "expected logger.error call"

