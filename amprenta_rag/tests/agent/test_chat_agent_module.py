from __future__ import annotations

from uuid import uuid4

from amprenta_rag.agent import chat_agent
from amprenta_rag.agent.chat_types import ChatSessionState


class FakeRagResult:
    def __init__(self, answer):
        self.answer = answer


def make_session():
    return ChatSessionState(turns=[])


def test_run_chat_turn_freeform(monkeypatch):
    monkeypatch.setattr(chat_agent, "route_intent", lambda text: "freeform_rag")
    monkeypatch.setattr(chat_agent, "query_rag", lambda **k: FakeRagResult("hi"))
    session, answer = chat_agent.run_chat_turn(make_session(), "hello")
    assert answer == "hi"
    assert session.turns


def test_dataset_summary_requires_id(monkeypatch):
    monkeypatch.setattr(chat_agent, "route_intent", lambda text: "dataset_summary")
    session, answer = chat_agent.run_chat_turn(make_session(), "dataset: not-a-uuid")
    assert "Please specify a dataset ID" in answer

    ds_id = str(uuid4())
    monkeypatch.setattr(chat_agent, "cross_omics_dataset_summary_postgres", lambda uid: f"ds-{uid}")
    session, answer = chat_agent.run_chat_turn(make_session(), f"dataset: {ds_id}")
    assert f"ds-{ds_id}" in answer


def test_program_signature_feature(monkeypatch):
    monkeypatch.setattr(chat_agent, "route_intent", lambda text: "program_summary")
    monkeypatch.setattr(chat_agent, "cross_omics_program_summary_postgres", lambda uid: "program-summary")
    _, ans = chat_agent.run_chat_turn(make_session(), f"program: {uuid4()}")
    assert "program-summary" in ans

    monkeypatch.setattr(chat_agent, "route_intent", lambda text: "signature_summary")
    monkeypatch.setattr(chat_agent, "cross_omics_signature_summary_postgres", lambda uid: "sig-summary")
    _, ans = chat_agent.run_chat_turn(make_session(), f"signature: {uuid4()}")
    assert "sig-summary" in ans

    monkeypatch.setattr(chat_agent, "route_intent", lambda text: "feature_summary")
    monkeypatch.setattr(chat_agent, "cross_omics_feature_summary_postgres", lambda uid, ft: "feat-summary")
    _, ans = chat_agent.run_chat_turn(make_session(), f"feature: {uuid4()}")
    assert "feat-summary" in ans


def test_help_and_unknown(monkeypatch):
    monkeypatch.setattr(chat_agent, "route_intent", lambda text: "help")
    _, ans = chat_agent.run_chat_turn(make_session(), "help")
    assert "I can help you" in ans

    monkeypatch.setattr(chat_agent, "route_intent", lambda text: "unknown")
    _, ans = chat_agent.run_chat_turn(make_session(), "??")
    assert "I'm not sure" in ans


def test_similar_datasets_branch(monkeypatch):
    monkeypatch.setattr(chat_agent, "route_intent", lambda text: "similar_datasets")
    _, ans = chat_agent.run_chat_turn(make_session(), "find similar")
    assert "Similarity" in ans


def test_error_handling(monkeypatch):
    monkeypatch.setattr(chat_agent, "route_intent", lambda text: "freeform_rag")

    def boom(**k):
        raise RuntimeError("bad")

    monkeypatch.setattr(chat_agent, "query_rag", boom)
    _, ans = chat_agent.run_chat_turn(make_session(), "hi")
    assert ans.startswith("Error:")

