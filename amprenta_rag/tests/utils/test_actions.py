from __future__ import annotations

import pytest

from amprenta_rag.utils import actions


def test_search_actions_empty_returns_all():
    result = actions.search_actions("")
    assert result == actions.QUICK_ACTIONS


def test_search_actions_matches_name_case_insensitive():
    result = actions.search_actions("experiment")
    names = [r["name"] for r in result]
    assert any("Experiment" in n for n in names)


def test_search_actions_matches_description_partial():
    result = actions.search_actions("validation")
    assert any(r["page"] == "Data Quality" for r in result)


def test_search_actions_no_match_returns_empty():
    assert actions.search_actions("nonexistentquery") == []

