from __future__ import annotations

import types
from unittest.mock import Mock

import pytest
from fastapi.testclient import TestClient


@pytest.fixture()
def api_client():
    from amprenta_rag.api.main import app

    return TestClient(app)


@pytest.mark.api
def test_post_generate_report_returns_200_with_valid_entity(api_client, monkeypatch):
    from amprenta_rag.api.routers import reports as reports_router

    fn = Mock(return_value="/tmp/report.html")
    monkeypatch.setattr(reports_router, "generate_report", fn)

    payload = {"entity_type": "dataset", "entity_id": "123", "format": "html"}
    resp = api_client.post("/api/v1/reports/generate", json=payload)
    assert resp.status_code == 200
    assert resp.json() == {"file_path": "/tmp/report.html", "download_url": None}

    fn.assert_called_once_with(
        template_name="narrative_report.ipynb",
        params={"entity_type": "dataset", "entity_id": "123", "format": "html"},
        format="html",
    )


@pytest.mark.api
def test_post_generate_report_returns_404_for_invalid_entity_id(api_client, monkeypatch):
    from amprenta_rag.api.routers import reports as reports_router

    monkeypatch.setattr(
        reports_router, "generate_report", Mock(side_effect=FileNotFoundError("Entity not found"))
    )

    payload = {"entity_type": "dataset", "entity_id": "does-not-exist", "format": "html"}
    resp = api_client.post("/api/v1/reports/generate", json=payload)
    assert resp.status_code == 404
    assert resp.json()["detail"] == "Entity not found"


@pytest.mark.api
def test_post_generate_report_returns_400_for_invalid_entity_type(api_client, monkeypatch):
    from amprenta_rag.api.routers import reports as reports_router

    monkeypatch.setattr(
        reports_router, "generate_report", Mock(side_effect=ValueError("Invalid entity_type"))
    )

    payload = {"entity_type": "not-a-real-type", "entity_id": "123", "format": "html"}
    resp = api_client.post("/api/v1/reports/generate", json=payload)
    assert resp.status_code == 400
    assert resp.json()["detail"] == "Invalid entity_type"


@pytest.mark.unit
def test_generate_report_calls_execute_notebook_and_convert_html(monkeypatch, tmp_path):
    """
    Validate generate_report() orchestration:
    - resolves template path
    - creates output dir
    - calls execute_notebook(template, executed_notebook, params)
    - calls converter for the requested format
    """
    from amprenta_rag.reports import generator

    templates_dir = str(tmp_path / "templates")
    template_name = "narrative_report.ipynb"
    template_path = str(tmp_path / "templates" / template_name)

    # Pretend the template exists regardless of filesystem
    monkeypatch.setattr(generator.os.path, "exists", lambda p: p == template_path)

    makedirs = Mock()
    monkeypatch.setattr(generator.os, "makedirs", makedirs)

    exec_nb = Mock(return_value=str(tmp_path / "templates" / "output" / "executed_narrative_report.ipynb"))
    monkeypatch.setattr(generator, "execute_notebook", exec_nb)

    convert_html = Mock(return_value=str(tmp_path / "report.html"))
    monkeypatch.setattr(generator, "convert_to_html", convert_html)

    out = generator.generate_report(
        template_name=template_name,
        params={"entity_type": "dataset", "entity_id": "123", "format": "html"},
        format="html",
        templates_dir=templates_dir,
    )

    output_dir = str(tmp_path / "templates" / "output")
    executed = str(tmp_path / "templates" / "output" / "executed_narrative_report.ipynb")

    makedirs.assert_called_once_with(output_dir, exist_ok=True)
    exec_nb.assert_called_once_with(template_path, executed, {"entity_type": "dataset", "entity_id": "123", "format": "html"})
    convert_html.assert_called_once_with(executed, executed)
    assert out == str(tmp_path / "report.html")


@pytest.mark.unit
def test_execute_notebook_calls_papermill_correctly(monkeypatch, tmp_path):
    from amprenta_rag.reports import generator

    pm = types.ModuleType("papermill")
    pm.execute_notebook = Mock()
    monkeypatch.setitem(__import__("sys").modules, "papermill", pm)

    template = str(tmp_path / "t.ipynb")
    out = str(tmp_path / "o.ipynb")
    params = {"entity_type": "dataset", "entity_id": "123"}

    assert generator.execute_notebook(template, out, params) == out
    pm.execute_notebook.assert_called_once_with(
        template,
        out,
        parameters=params,
        kernel_name="python3",
    )


