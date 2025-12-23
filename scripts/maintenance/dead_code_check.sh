#!/usr/bin/env bash
set -euo pipefail

ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
REPORTS_DIR="${ROOT_DIR}/reports"
WHITELIST="${ROOT_DIR}/scripts/maintenance/vulture_whitelist.py"

mkdir -p "${REPORTS_DIR}"

RAW_OUT="${REPORTS_DIR}/dead_code_vulture.txt"
SUMMARY_OUT="${REPORTS_DIR}/dead_code_summary.md"

NOW_UTC="$(date -u +"%Y-%m-%dT%H:%M:%SZ")"

if ! command -v vulture >/dev/null 2>&1; then
  {
    echo "vulture not found on PATH."
    echo "Install requirements and rerun (analysis-only; no changes made)."
  } > "${RAW_OUT}"

  {
    echo "# Dead Code Summary"
    echo
    echo "- Generated: ${NOW_UTC}"
    echo "- Working dir: ${ROOT_DIR}"
    echo
    echo "**Status:** FAILED (vulture not installed)"
  } > "${SUMMARY_OUT}"

  echo "Wrote:"
  echo "  - ${RAW_OUT}"
  echo "  - ${SUMMARY_OUT}"
  exit 0
fi

# Raw vulture output (min-confidence 60)
vulture amprenta_rag/ scripts/ "${WHITELIST}" --min-confidence 60 > "${RAW_OUT}" || true

# Build a filtered "high confidence" report + include unused-import count (ruff F401)
ROOT_DIR="${ROOT_DIR}" NOW_UTC="${NOW_UTC}" python - <<'PY'
import os
import re
import subprocess
from pathlib import Path

root = Path(os.environ.get("ROOT_DIR", ".")).resolve()
reports = root / "reports"
raw_path = reports / "dead_code_vulture.txt"
summary_path = reports / "dead_code_summary.md"

raw = raw_path.read_text(encoding="utf-8", errors="replace").splitlines()

# Vulture formats vary:
# - "(confidence: 60%)"
# - "(60% confidence)"
conf_re = re.compile(r"\((?:confidence:\s*)?(\d+)%\s*(?:confidence)?\)\s*$")

entries = []
for ln in raw:
    m = conf_re.search(ln)
    if not m:
        continue
    conf = int(m.group(1))
    entries.append((conf, ln))

entries_sorted = sorted(entries, key=lambda t: (-t[0], t[1]))
hi = [ln for conf, ln in entries_sorted if conf >= 80]

def run(cmd):
    try:
        p = subprocess.run(cmd, cwd=str(root), text=True, capture_output=True)
        return p.returncode, (p.stdout or "") + (p.stderr or "")
    except Exception as e:  # noqa: BLE001
        return 1, str(e)

ruff_code, ruff_out = run(["ruff", "check", "--select", "F401", "amprenta_rag", "scripts"])
unused_imports_count = 0
if ruff_code in (0, 1):  # 1 = findings
    unused_imports_count = sum(1 for ln in ruff_out.splitlines() if "F401" in ln)

lines = []
lines.append("# Dead Code Summary")
lines.append("")
lines.append(f"- Generated: {os.environ.get('NOW_UTC', '')}")
lines.append(f"- Working dir: {root}")
lines.append("")
lines.append("## Counts")
lines.append("")
lines.append(f"- Vulture candidates (confidence >= 60): **{len(entries)}**")
lines.append(f"- High-confidence candidates (confidence >= 80): **{len(hi)}**")
lines.append(f"- Ruff unused imports (F401): **{unused_imports_count}**")
lines.append("")
lines.append("## High-Confidence Candidates (>= 80%)")
lines.append("")
if hi:
    lines.append("```")
    lines.extend(hi[:200])
    if len(hi) > 200:
        lines.append(f"... truncated ({len(hi) - 200} more)")
    lines.append("```")
else:
    lines.append("_None._")
    lines.append("")

summary_path.write_text("\n".join([ln for ln in lines if ln != ""]) + "\n", encoding="utf-8")
PY

echo "Wrote:"
echo "  - ${RAW_OUT}"
echo "  - ${SUMMARY_OUT}"


