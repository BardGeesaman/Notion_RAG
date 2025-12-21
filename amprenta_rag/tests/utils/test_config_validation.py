from __future__ import annotations

from amprenta_rag.utils import config_validation


class PipelineCfg:
    def __init__(
        self,
        use_postgres_as_sot=False,
        enable_feature_linking=False,
        feature_linking_max_workers=10,
        enable_notion_sync=False,
        enable_dual_write=False,
    ):
        self.use_postgres_as_sot = use_postgres_as_sot
        self.enable_feature_linking = enable_feature_linking
        self.feature_linking_max_workers = feature_linking_max_workers
        self.enable_notion_sync = enable_notion_sync
        self.enable_dual_write = enable_dual_write


class DatabaseCfg:
    def __init__(self, host="", db="", user=""):
        self.postgres_host = host
        self.postgres_db = db
        self.postgres_user = user


class Cfg:
    def __init__(self):
        self.pipeline = PipelineCfg()
        self.database = DatabaseCfg()


def test_validate_postgres_config_flags(monkeypatch):
    cfg = Cfg()
    cfg.pipeline.use_postgres_as_sot = True
    cfg.database = DatabaseCfg(host="", db="", user="")
    monkeypatch.setattr(config_validation, "get_config", lambda: cfg)

    issues = config_validation.validate_postgres_config()
    assert len(issues) == 3


def test_validate_feature_linking_workers(monkeypatch):
    cfg = Cfg()
    cfg.pipeline.enable_feature_linking = True
    cfg.pipeline.feature_linking_max_workers = 0
    monkeypatch.setattr(config_validation, "get_config", lambda: cfg)

    issues = config_validation.validate_feature_linking_config()
    assert issues and "FEATURE_LINKING_MAX_WORKERS" in issues[0]


def test_validate_notion_sync_missing_key(monkeypatch):
    cfg = Cfg()
    cfg.pipeline.enable_notion_sync = True
    monkeypatch.setattr(config_validation, "get_config", lambda: cfg)
    monkeypatch.setattr(config_validation, "NOTION_API_KEY", "", raising=False)

    issues = config_validation.validate_notion_sync_config()
    assert issues and "NOTION_API_KEY" in issues[0]


def test_validate_configuration_result(monkeypatch):
    cfg = Cfg()
    cfg.pipeline.use_postgres_as_sot = True
    monkeypatch.setattr(config_validation, "get_config", lambda: cfg)
    monkeypatch.setattr(config_validation, "NOTION_API_KEY", "", raising=False)

    result = config_validation.validate_configuration()
    assert result.success is False
    assert result.errors

