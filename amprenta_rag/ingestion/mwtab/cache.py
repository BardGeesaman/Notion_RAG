"""Simple file-based cache for mwTab JSON responses."""

from __future__ import annotations

import json
import os
from pathlib import Path
from typing import Any, Dict, Optional
from datetime import datetime, timedelta
import re
import tempfile

from amprenta_rag.logging_utils import get_logger

logger = get_logger(__name__)

CACHE_DIR = Path(".cache/mwtab")
TTL_DAYS = 7


def _ensure_cache_dir() -> bool:
    try:
        CACHE_DIR.mkdir(parents=True, exist_ok=True)
        return True
    except Exception:
        return False


def _cache_path(study_id: str) -> Path:
    return CACHE_DIR / f"{_sanitize_key(study_id)}.json"


def _sanitize_key(key: str) -> str:
    cleaned = re.sub(r"[^A-Za-z0-9_-]", "", key or "")
    if not cleaned:
        cleaned = "unknown"
    return cleaned[:100]


def is_expired(path: Path, ttl_days: int = TTL_DAYS) -> bool:
    try:
        stat = path.stat()
        mtime = datetime.utcfromtimestamp(stat.st_mtime)
        return datetime.utcnow() - mtime > timedelta(days=ttl_days)
    except Exception:
        return True


def get_cached(study_id: str) -> Optional[Dict[str, Any]]:
    path = _cache_path(study_id)
    if not path.exists():
        return None
    if is_expired(path, TTL_DAYS):
        return None
    try:
        with path.open("r", encoding="utf-8") as f:
            return json.load(f)
    except Exception:
        return None


def set_cached(study_id: str, data: Dict[str, Any]) -> None:
    if not _ensure_cache_dir():
        return
    path = _cache_path(study_id)
    try:
        with tempfile.NamedTemporaryFile("w", delete=False, encoding="utf-8", dir=str(CACHE_DIR)) as tmp:
            json.dump(data, tmp)
            tmp_path = Path(tmp.name)
        os.replace(tmp_path, path)
    except Exception as exc:
        logger.debug("Failed to write cache for %s: %r", study_id, exc)
        try:
            if tmp_path and tmp_path.exists():
                tmp_path.unlink()
        except Exception:
            pass

