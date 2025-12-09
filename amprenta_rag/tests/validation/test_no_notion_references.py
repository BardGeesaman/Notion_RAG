import pathlib


def test_no_notion_references():
    """Ensure deprecated Notion symbols are absent from code and docs."""
    root = pathlib.Path(__file__).resolve().parents[3]  # repo root
    patterns = [
        "notion_headers(",
        "from amprenta_rag.clients.notion_client",
        "import notion_client",
        "api.notion.com",
    ]
    excluded = {
        ".git",
        ".pytest_cache",
        "__pycache__",
        "md_archive",
        "docs",
        "context",
        "ChatGPT context",
        "tests",
    }

    for path in root.rglob("*"):
        if path.is_dir():
            if path.name in excluded:
                # skip excluded dirs entirely
                dirs = []
            continue
        if not path.is_file():
            continue
        if any(part in excluded for part in path.parts):
            continue
        try:
            text = path.read_text(encoding="utf-8", errors="ignore")
        except Exception:
            continue
        lower = text.lower()
        for pat in patterns:
            if pat.lower() in lower:
                raise AssertionError(f"Found deprecated Notion reference '{pat}' in {path}")

