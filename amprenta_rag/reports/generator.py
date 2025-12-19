"""Report generator using Papermill and nbconvert."""

import os
from pathlib import Path
from typing import Dict, Any, Optional

from nbconvert import HTMLExporter, PDFExporter
from nbformat import read


def execute_notebook(
    template_path: str,
    output_path: str,
    parameters: Optional[Dict[str, Any]] = None
) -> str:
    """Execute a notebook template with parameters using Papermill.

    Args:
        template_path: Path to the notebook template (.ipynb)
        output_path: Path where executed notebook will be saved
        parameters: Dictionary of parameters to inject into notebook

    Returns:
        Path to the executed notebook
    """
    if parameters is None:
        parameters = {}

    try:
        import papermill as pm  # type: ignore
    except ModuleNotFoundError as exc:
        raise ImportError(
            "papermill is required for report generation. Install papermill."
        ) from exc

    pm.execute_notebook(
        template_path,
        output_path,
        parameters=parameters,
        kernel_name="python3"
    )

    return output_path


def convert_to_pdf(notebook_path: str, output_path: str) -> str:
    """Convert executed notebook to PDF using nbconvert.

    Args:
        notebook_path: Path to executed notebook (.ipynb)
        output_path: Path where PDF will be saved

    Returns:
        Path to the PDF file
    """
    exporter = PDFExporter()
    exporter.exclude_input_prompt = True
    exporter.exclude_output_prompt = True

    with open(notebook_path, "r") as f:
        notebook = read(f, as_version=4)

    (body, resources) = exporter.from_notebook_node(notebook)

    pdf_path = Path(output_path).with_suffix(".pdf")
    with open(pdf_path, "wb") as f:
        f.write(body)

    return str(pdf_path)


def convert_to_html(notebook_path: str, output_path: str) -> str:
    """Convert executed notebook to HTML using nbconvert.

    Args:
        notebook_path: Path to executed notebook (.ipynb)
        output_path: Path where HTML will be saved

    Returns:
        Path to the HTML file
    """
    exporter = HTMLExporter()
    exporter.exclude_input_prompt = True
    exporter.exclude_output_prompt = True

    with open(notebook_path, "r") as f:
        notebook = read(f, as_version=4)

    (body, resources) = exporter.from_notebook_node(notebook)

    html_path = Path(output_path).with_suffix(".html")
    with open(html_path, "w") as f:
        f.write(body)

    return str(html_path)


def generate_report(
    template_name: str,
    params: Dict[str, Any],
    format: str = "html",
    templates_dir: Optional[str] = None
) -> str:
    """Generate a report from a notebook template.

    Args:
        template_name: Name of template notebook (e.g., "getting_started.ipynb")
        params: Parameters to inject into notebook
        format: Output format ("html" or "pdf")
        templates_dir: Directory containing templates (default: deploy/jupyterhub/templates)

    Returns:
        Path to generated report
    """
    # Security: Reject path separators to prevent path traversal
    if "/" in template_name or "\\" in template_name or ".." in template_name:
        raise ValueError(
            f"Invalid template_name: '{template_name}'. "
            "Path separators and '..' are not allowed."
        )

    if templates_dir is None:
        # Default to templates directory relative to project root
        templates_dir = os.path.join(
            os.path.dirname(os.path.dirname(os.path.dirname(__file__))),
            "deploy",
            "jupyterhub",
            "templates"
        )

    template_path = os.path.join(templates_dir, template_name)
    if not os.path.exists(template_path):
        raise FileNotFoundError(f"Template not found: {template_path}")

    # Execute notebook
    output_dir = os.path.join(os.path.dirname(template_path), "output")
    os.makedirs(output_dir, exist_ok=True)

    executed_notebook = os.path.join(
        output_dir,
        f"executed_{Path(template_name).stem}.ipynb"
    )

    execute_notebook(template_path, executed_notebook, params)

    # Convert to requested format
    if format == "pdf":
        report_path = convert_to_pdf(executed_notebook, executed_notebook)
    elif format == "html":
        report_path = convert_to_html(executed_notebook, executed_notebook)
    else:
        raise ValueError(f"Unsupported format: {format}. Use 'html' or 'pdf'")

    return report_path

