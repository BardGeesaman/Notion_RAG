from __future__ import annotations

import json

from amprenta_rag.analysis import literature_critique as lc


class _FakeChoice:
    def __init__(self, content: str):
        self.message = type("msg", (), {"content": content})


class _FakeClient:
    def __init__(self, payload: dict | None = None, should_raise: bool = False):
        self._payload = payload or {}
        self._raise = should_raise
        self.chat = type("chat", (), {"completions": self})

    def create(self, **_kwargs):
        if self._raise:
            raise RuntimeError("boom")
        return type("Resp", (), {"choices": [_FakeChoice(json.dumps(self._payload))]})


def test_generate_critique_success(monkeypatch) -> None:
    payload = {
        "strengths": ["good"],
        "weaknesses": ["bad"],
        "limitations": ["few"],
        "methodology_score": 80,
    }
    monkeypatch.setattr(lc, "get_openai_client", lambda: _FakeClient(payload))
    monkeypatch.setattr(lc, "get_default_models", lambda: ("gpt", None))

    result = lc.generate_critique("text")
    assert result["methodology_score"] == 80
    assert result["strengths"] == ["good"]


def test_generate_critique_error(monkeypatch) -> None:
    monkeypatch.setattr(lc, "get_openai_client", lambda: _FakeClient(should_raise=True))
    monkeypatch.setattr(lc, "get_default_models", lambda: ("gpt", None))

    result = lc.generate_critique("text")
    assert result["methodology_score"] == 50
    assert result["strengths"] == []


def test_extract_unanswered_questions(monkeypatch) -> None:
    payload = {"questions": ["Q1", "Q2"]}
    monkeypatch.setattr(lc, "get_openai_client", lambda: _FakeClient(payload))
    monkeypatch.setattr(lc, "get_default_models", lambda: ("gpt", None))

    questions = lc.extract_unanswered_questions("text")
    assert questions == ["Q1", "Q2"]


def test_detect_contradictions(monkeypatch) -> None:
    payload = {
        "contradictions": [
            {"claim_a": "A", "claim_b": "B", "topic": "T", "severity": "low"},
            {"claim_a": "A2", "claim_b": "B2", "topic": "T2", "severity": "high"},
        ]
    }
    monkeypatch.setattr(lc, "get_openai_client", lambda: _FakeClient(payload))
    monkeypatch.setattr(lc, "get_default_models", lambda: ("gpt", None))

    contradictions = lc.detect_contradictions("t1", "t2")
    assert len(contradictions) == 2
    assert contradictions[0]["topic"] == "T"


def test_detect_contradictions_error(monkeypatch) -> None:
    monkeypatch.setattr(lc, "get_openai_client", lambda: _FakeClient(should_raise=True))
    monkeypatch.setattr(lc, "get_default_models", lambda: ("gpt", None))

    contradictions = lc.detect_contradictions("t1", "t2")
    assert contradictions == []

