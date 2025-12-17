#!/usr/bin/env python3
"""
Harvest mwTab studies from Metabolomics Workbench with optional parallel fetching.

Provides `fetch_mw_mwtab(study_id)` for reuse (e.g., ingestion pipelines) and a CLI
to fetch one or more study IDs, optionally writing results to an output directory.
"""

from __future__ import annotations

import argparse
import json
import time
import threading
import sys
from pathlib import Path
from typing import Dict, List, Optional
from concurrent.futures import ThreadPoolExecutor, as_completed

import requests

# Ensure repo root on sys.path for script execution
REPO_ROOT = Path(__file__).resolve().parent.parent
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

from amprenta_rag.logging_utils import get_logger

logger = get_logger(__name__)

MW_BASE_URL = "https://www.metabolomicsworkbench.org/rest"
_rate_lock = threading.Lock()
_last_request_time = 0.0
_RATE_SECONDS = 0.5


def _rate_limited_get(url: str, timeout: int = 30):
    global _last_request_time
    with _rate_lock:
        now = time.time()
        elapsed = now - _last_request_time
        if elapsed < _RATE_SECONDS:
            time.sleep(_RATE_SECONDS - elapsed)
        _last_request_time = time.time()
    return requests.get(url, timeout=timeout)


def fetch_mw_mwtab(study_id: str) -> str:
    """
    Fetch mwTab data for a study.

    Args:
        study_id: MW study ID (e.g., "ST001111")

    Returns:
        Raw mwTab text (or JSON string). Empty string on failure.
    """
    url = f"{MW_BASE_URL}/study/study_id/{study_id}/mwtab"
    try:
        resp = _rate_limited_get(url, timeout=30)
        resp.raise_for_status()
        text_content = resp.text
        content_type = resp.headers.get("content-type", "")
        if "application/json" in content_type:
            try:
                data = resp.json()
                return json.dumps(data, indent=2)
            except Exception:
                logger.debug("[MW] Response looks JSON but parse failed; returning raw text")
        return text_content
    except Exception as e:
        logger.warning("[MW] Error fetching mwTab for %s: %r", study_id, e)
        try:
            alt_url = f"{MW_BASE_URL}/study/study_id/{study_id}/data"
            resp = _rate_limited_get(alt_url, timeout=30)
            resp.raise_for_status()
            return resp.text
        except Exception as e2:
            logger.warning("[MW] Alternative fetch failed for %s: %r", study_id, e2)
            return ""


def _write_output(study_id: str, content: str, out_dir: Optional[Path]) -> None:
    if not out_dir:
        return
    try:
        out_dir.mkdir(parents=True, exist_ok=True)
        path = out_dir / f"{study_id}.mwtab"
        path.write_text(content, encoding="utf-8")
    except Exception as exc:
        logger.warning("[MW] Failed to write output for %s: %r", study_id, exc)


def _process_study(study_id: str, out_dir: Optional[Path]) -> Dict[str, str]:
    content = fetch_mw_mwtab(study_id)
    if content:
        _write_output(study_id, content, out_dir)
        return {"study_id": study_id, "status": "ok"}
    return {"study_id": study_id, "status": "failed"}


def _load_ids_from_file(path: Path) -> List[str]:
    if not path.exists():
        logger.warning("[MW] study-id-file not found: %s", path)
        return []
    return [line.strip() for line in path.read_text().splitlines() if line.strip()]


def main() -> None:
    parser = argparse.ArgumentParser(description="Fetch mwTab studies from Metabolomics Workbench.")
    parser.add_argument("--study-id", action="append", help="Study ID to fetch (can repeat).")
    parser.add_argument("--study-id-file", type=Path, help="File with study IDs (one per line).")
    parser.add_argument("--output-dir", type=Path, help="Directory to write mwTab files.")
    parser.add_argument("--workers", type=int, default=4, help="Number of parallel workers (default: 4). Use 1 for sequential.")
    args = parser.parse_args()

    study_ids: List[str] = []
    if args.study_id:
        study_ids.extend(args.study_id)
    if args.study_id_file:
        study_ids.extend(_load_ids_from_file(args.study_id_file))

    study_ids = sorted(set([sid.strip() for sid in study_ids if sid and sid.strip()]))

    if not study_ids:
        parser.error("No study IDs provided. Use --study-id or --study-id-file.")

    logger.info("[MW] Fetching %d studies with %d worker(s)", len(study_ids), args.workers)

    results: List[Dict[str, str]] = []
    if args.workers <= 1:
        for sid in study_ids:
            results.append(_process_study(sid, args.output_dir))
    else:
        with ThreadPoolExecutor(max_workers=args.workers) as executor:
            future_map = {executor.submit(_process_study, sid, args.output_dir): sid for sid in study_ids}
            for future in as_completed(future_map):
                sid = future_map[future]
                try:
                    res = future.result()
                except Exception as exc:
                    logger.warning("[MW] Worker failed for %s: %r", sid, exc)
                    res = {"study_id": sid, "status": "failed"}
                results.append(res)

    ok = sum(1 for r in results if r["status"] == "ok")
    fail = len(results) - ok
    logger.info("[MW] Done. Success=%d, Failed=%d", ok, fail)


if __name__ == "__main__":
    main()

