from __future__ import annotations

import sys
from pathlib import Path
from types import ModuleType

import pytest

from amprenta_rag.reports import generator as gen


def test_execute_notebook_requires_papermill(monkeypatch, tmp_path):
    # Force the import inside execute_notebook() to fail even if papermill is installed.
    import builtins

    real_import = builtins.__import__

    def fake_import(name, globals=None, locals=None, fromlist=(), level=0):
        if name == "papermill":
            raise ModuleNotFoundError("No module named 'papermill'")
        return real_import(name, globals, locals, fromlist, level)

    monkeypatch.setattr(builtins, "__import__", fake_import)

    with pytest.raises(ImportError, match="papermill is required"):
        gen.execute_notebook(
            template_path=str(tmp_path / "t.ipynb"),
            output_path=str(tmp_path / "out.ipynb"),
            parameters={"a": 1},
        )


def test_execute_notebook_calls_papermill(monkeypatch, tmp_path):
    called: dict[str, object] = {}

    fake_pm = ModuleType("papermill")

    def execute_notebook(template_path, output_path, parameters=None, kernel_name=None):
        called["template_path"] = template_path
        called["output_path"] = output_path
        called["parameters"] = parameters
        called["kernel_name"] = kernel_name

    fake_pm.execute_notebook = execute_notebook  # type: ignore[attr-defined]
    monkeypatch.setitem(sys.modules, "papermill", fake_pm)

    tpl = tmp_path / "t.ipynb"
    out = tmp_path / "out.ipynb"
    tpl.write_text("{}")

    res = gen.execute_notebook(str(tpl), str(out), parameters={"x": "y"})
    assert res == str(out)
    assert called["template_path"] == str(tpl)
    assert called["output_path"] == str(out)
    assert called["parameters"] == {"x": "y"}
    assert called["kernel_name"] == "python3"


def test_convert_to_pdf_writes_bytes(monkeypatch, tmp_path):
    nb = tmp_path / "n.ipynb"
    nb.write_text("{}")

    monkeypatch.setattr(gen, "read", lambda f, as_version=4: {"nb": True})

    class FakePDFExporter:
        exclude_input_prompt = False
        exclude_output_prompt = False

        def from_notebook_node(self, notebook):
            return (b"PDFDATA", {})

    monkeypatch.setattr(gen, "PDFExporter", FakePDFExporter)

    out = tmp_path / "report"
    pdf_path = gen.convert_to_pdf(str(nb), str(out))
    assert pdf_path.endswith(".pdf")
    assert Path(pdf_path).read_bytes() == b"PDFDATA"


def test_convert_to_html_writes_text(monkeypatch, tmp_path):
    nb = tmp_path / "n.ipynb"
    nb.write_text("{}")

    monkeypatch.setattr(gen, "read", lambda f, as_version=4: {"nb": True})

    class FakeHTMLExporter:
        exclude_input_prompt = False
        exclude_output_prompt = False

        def from_notebook_node(self, notebook):
            return ("<html>OK</html>", {})

    monkeypatch.setattr(gen, "HTMLExporter", FakeHTMLExporter)

    out = tmp_path / "report"
    html_path = gen.convert_to_html(str(nb), str(out))
    assert html_path.endswith(".html")
    assert Path(html_path).read_text() == "<html>OK</html>"


@pytest.mark.parametrize("bad_name", ["../x.ipynb", "a/../b.ipynb", "a\\b.ipynb", "a/b.ipynb"])
def test_generate_report_rejects_path_traversal(bad_name, tmp_path):
    with pytest.raises(ValueError, match="Invalid template_name"):
        gen.generate_report(template_name=bad_name, params={}, templates_dir=str(tmp_path))


def test_generate_report_missing_template(tmp_path):
    with pytest.raises(FileNotFoundError, match="Template not found"):
        gen.generate_report(template_name="missing.ipynb", params={}, templates_dir=str(tmp_path))


def test_generate_report_html_and_pdf(monkeypatch, tmp_path):
    templates = tmp_path / "templates"
    templates.mkdir()
    tpl = templates / "t.ipynb"
    tpl.write_text("{}")

    # Avoid actually running papermill/nbconvert
    def fake_execute(template_path: str, output_path: str, parameters=None):
        Path(output_path).write_text("{}")
        return output_path

    monkeypatch.setattr(gen, "execute_notebook", fake_execute)
    monkeypatch.setattr(gen, "convert_to_html", lambda notebook_path, output_path: str(Path(output_path).with_suffix(".html")))
    monkeypatch.setattr(gen, "convert_to_pdf", lambda notebook_path, output_path: str(Path(output_path).with_suffix(".pdf")))

    html = gen.generate_report("t.ipynb", params={"a": 1}, format="html", templates_dir=str(templates))
    assert html.endswith(".html")

    pdf = gen.generate_report("t.ipynb", params={"a": 1}, format="pdf", templates_dir=str(templates))
    assert pdf.endswith(".pdf")


def test_generate_report_unsupported_format(monkeypatch, tmp_path):
    templates = tmp_path / "templates"
    templates.mkdir()
    (templates / "t.ipynb").write_text("{}")

    monkeypatch.setattr(gen, "execute_notebook", lambda *a, **k: str(tmp_path / "x.ipynb"))

    with pytest.raises(ValueError, match="Unsupported format"):
        gen.generate_report("t.ipynb", params={}, format="docx", templates_dir=str(templates))


