from __future__ import annotations



from amprenta_rag.utils import config_check


def test_validate_config_all_set(monkeypatch, capsys):
    for var in config_check.REQUIRED_VARS:
        monkeypatch.setenv(var, "value")

    assert config_check.validate_config() is True
    out = capsys.readouterr().out
    for var in config_check.REQUIRED_VARS:
        assert f"{var} is set" in out


def test_validate_config_missing(monkeypatch, capsys):
    # Ensure clean env
    for var in config_check.REQUIRED_VARS:
        monkeypatch.delenv(var, raising=False)

    monkeypatch.setenv("DATABASE_URL", "db")
    result = config_check.validate_config()
    out = capsys.readouterr().out
    assert result is False
    assert "DATABASE_URL is set" in out
    for var in config_check.REQUIRED_VARS:
        if var != "DATABASE_URL":
            assert f"{var} is NOT set" in out

