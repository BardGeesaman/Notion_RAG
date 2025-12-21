from __future__ import annotations

import json

import pytest

from amprenta_rag.analysis import design_engine as de


class _FakeChoice:
    def __init__(self, content: str):
        self.message = type("msg", (), {"content": content})


class _FakeChat:
    def __init__(self, response: dict):
        self._response = response

        class _Completions:
            def __init__(self, outer):
                self._outer = outer

            def create(self, **_kwargs):
                return type("Resp", (), {"choices": [_FakeChoice(json.dumps(self._outer._response))]})

        class _Chat:
            def __init__(self, outer):
                self.completions = _Completions(outer)

        self.chat = _Chat(self)


class _ErrorCompletions:
    def create(self, **_kwargs):
        raise RuntimeError("boom")


class _ErrorClient:
    def __init__(self):
        self.chat = type("chat", (), {"completions": _ErrorCompletions()})


def test_recommend_design_success(monkeypatch: pytest.MonkeyPatch) -> None:
    resp = {
        "design_type": "case_control",
        "rationale": "Test rationale",
        "min_samples": 12,
        "considerations": ["c1"],
    }

    monkeypatch.setattr(de, "get_openai_client", lambda: _FakeChat(resp))
    monkeypatch.setattr(de, "get_default_models", lambda: ("gpt-test", None))

    result = de.recommend_design("RQ", 10, ["var1"])

    assert result["design_type"] == "case_control"
    assert result["rationale"] == "Test rationale"
    assert result["min_samples"] == 12
    assert result["considerations"] == ["c1"]


def test_recommend_design_error(monkeypatch: pytest.MonkeyPatch) -> None:
    monkeypatch.setattr(de, "get_openai_client", lambda: _ErrorClient())
    monkeypatch.setattr(de, "get_default_models", lambda: ("gpt-test", None))

    result = de.recommend_design("RQ", 5, [])

    assert result["design_type"] == "observational"
    assert "Error" in result["considerations"][0]


def test_get_design_requirements_default_fallback() -> None:
    reqs = de.get_design_requirements("unknown")
    assert reqs["description"]
    assert "min_samples" in reqs


def test_validate_design_warnings() -> None:
    warnings = de.validate_design("case_control", sample_count=10, group_count=1)
    assert any("Insufficient samples" in w for w in warnings)
    assert any("Insufficient groups" in w for w in warnings)

