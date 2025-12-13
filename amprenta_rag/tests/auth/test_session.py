from unittest.mock import patch

import pytest

from amprenta_rag.auth import session as session_mod


@pytest.fixture
def mock_session_state():
    with patch("amprenta_rag.auth.session.st") as mock_st:
        mock_st.session_state = {}
        yield mock_st.session_state


def test_set_and_get_current_user(mock_session_state):
    user = {"id": "u1", "username": "dev"}
    session_mod.set_current_user(user)
    assert session_mod.get_current_user() == user


def test_clear_session(mock_session_state):
    session_mod.set_current_user({"id": "u1"})
    assert session_mod.is_authenticated() is True
    session_mod.clear_session()
    assert session_mod.is_authenticated() is False
    assert session_mod.get_current_user() is None


def test_is_authenticated(mock_session_state):
    assert session_mod.is_authenticated() is False
    session_mod.set_current_user({"id": "u1"})
    assert session_mod.is_authenticated() is True


def test_check_session_timeout(mock_session_state):
    # Authenticate user
    session_mod.set_current_user({"id": "u1"})

    # If last_activity missing, it should initialize and not time out
    with patch("amprenta_rag.auth.session.time.time", return_value=1000.0):
        assert session_mod.check_session_timeout(timeout_minutes=30) is False
        assert "last_activity" in mock_session_state

    # If last_activity is old enough, should time out
    mock_session_state["last_activity"] = 0.0
    with patch("amprenta_rag.auth.session.time.time", return_value=(31 * 60) + 1.0):
        assert session_mod.check_session_timeout(timeout_minutes=30) is True


