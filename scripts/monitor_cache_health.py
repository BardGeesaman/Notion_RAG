#!/usr/bin/env python3
"""
Monitor dataset feature cache health.

Runs `manage_feature_cache.py --stats`, logs the statistics, and enforces
simple SLO-style thresholds on hit rate and eviction rate.

Intended to be run on a schedule (e.g., via cron) roughly hourly.
"""

from __future__ import annotations

import argparse
import re
import subprocess
import sys
from datetime import datetime
from pathlib import Path
from typing import Dict, Tuple

from amprenta_rag.logging_utils import get_logger


logger = get_logger(__name__)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Check dataset feature cache health using manage_feature_cache.py --stats",
    )
    parser.add_argument(
        "--min-hit-rate",
        type=float,
        default=0.50,
        help="Minimum acceptable cache hit rate (0.0–1.0, default: 0.50).",
    )
    parser.add_argument(
        "--max-eviction-rate",
        type=float,
        default=0.10,
        help="Maximum acceptable eviction rate (0.0–1.0, default: 0.10).",
    )
    return parser.parse_args()


def _run_manage_stats(repo_root: Path) -> Tuple[str, int]:
    """Invoke manage_feature_cache.py --stats and return (output, returncode)."""
    script_path = repo_root / "scripts" / "manage_feature_cache.py"

    cmd = [sys.executable, str(script_path), "--stats"]
    proc = subprocess.run(
        cmd,
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
        text=True,
    )
    output = proc.stdout or ""
    return output, proc.returncode


def _parse_stats_output(output: str) -> Dict[str, int]:
    """
    Parse the textual stats output from manage_feature_cache.py.

    Expects lines like:
      cached_datasets   12
      hits              123
      misses            45
      evictions         3
    """
    stats: Dict[str, int] = {}
    pattern = re.compile(r"^\s*(cached_datasets|hits|misses|evictions)\s+(\d+)\s*$")

    for line in output.splitlines():
        m = pattern.match(line)
        if not m:
            continue
        key, value = m.group(1), m.group(2)
        stats[key] = int(value)

    return stats


def _compute_rates(stats: Dict[str, int]) -> Tuple[float, float]:
    """
    Compute (hit_rate, eviction_rate) from raw stats.

    - hit_rate = hits / (hits + misses)             (0 if no traffic)
    - eviction_rate = evictions / (hits+misses+evictions)  (0 if no events)
    """
    hits = stats.get("hits", 0)
    misses = stats.get("misses", 0)
    evictions = stats.get("evictions", 0)

    traffic = hits + misses
    events = hits + misses + evictions

    hit_rate = (hits / traffic) if traffic > 0 else 0.0
    eviction_rate = (evictions / events) if events > 0 else 0.0
    return hit_rate, eviction_rate


def _append_log(repo_root: Path, text: str) -> None:
    logs_dir = repo_root / "logs"
    logs_dir.mkdir(parents=True, exist_ok=True)
    log_path = logs_dir / "cache_health.log"

    ts = datetime.utcnow().strftime("%Y-%m-%d %H:%M:%S")
    with log_path.open("a", encoding="utf-8") as f:
        for line in text.splitlines():
            f.write(f"[{ts}] {line}\n")


def main() -> None:
    args = parse_args()

    repo_root = Path(__file__).resolve().parents[1]

    output, rc = _run_manage_stats(repo_root)

    # Always log raw stats output for auditability
    if output:
        _append_log(repo_root, "Raw manage_feature_cache --stats output:")
        _append_log(repo_root, output)

    if rc != 0:
        msg = (
            f"manage_feature_cache.py --stats exited with code {rc}. "
            "Cache health check cannot proceed."
        )
        logger.error(msg)
        _append_log(repo_root, f"ERROR: {msg}")
        print(f"ALERT: {msg}")
        raise SystemExit(1)

    stats = _parse_stats_output(output)
    hit_rate, eviction_rate = _compute_rates(stats)

    cached = stats.get("cached_datasets", 0)
    hits = stats.get("hits", 0)
    misses = stats.get("misses", 0)
    evictions = stats.get("evictions", 0)

    summary = (
        "Cache health summary:\n"
        f"  cached_datasets   = {cached}\n"
        f"  hits              = {hits}\n"
        f"  misses            = {misses}\n"
        f"  evictions         = {evictions}\n"
        f"  hit_rate          = {hit_rate:.1%}\n"
        f"  eviction_rate     = {eviction_rate:.1%}\n"
        f"  min_hit_rate      = {args.min_hit_rate:.1%}\n"
        f"  max_eviction_rate = {args.max_eviction_rate:.1%}"
    )
    print(summary)
    _append_log(repo_root, summary)

    unhealthy_reasons = []
    if hit_rate < args.min_hit_rate:
        unhealthy_reasons.append(
            f"hit_rate {hit_rate:.1%} is below threshold {args.min_hit_rate:.1%}"
        )
    if eviction_rate > args.max_eviction_rate:
        unhealthy_reasons.append(
            f"eviction_rate {eviction_rate:.1%} exceeds threshold {args.max_eviction_rate:.1%}"
        )

    if unhealthy_reasons:
        alert_msg = "Cache health check FAILED: " + "; ".join(unhealthy_reasons)
        logger.warning(alert_msg)
        _append_log(repo_root, f"ALERT: {alert_msg}")
        print(f"\nALERT: {alert_msg}")
        raise SystemExit(1)

    ok_msg = "Cache health check OK - within configured thresholds."
    logger.info(ok_msg)
    _append_log(repo_root, ok_msg)
    print(f"\n{ok_msg}")


if __name__ == "__main__":
    main()


