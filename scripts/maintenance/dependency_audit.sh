#!/usr/bin/env bash
set -euo pipefail

ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
REPORTS_DIR="${ROOT_DIR}/reports"

mkdir -p "${REPORTS_DIR}"

SEC_JSON="${REPORTS_DIR}/security_audit.json"
SEC_MD="${REPORTS_DIR}/security_audit.md"
OUTDATED_JSON="${REPORTS_DIR}/outdated.json"
OUTDATED_SUMMARY_MD="${REPORTS_DIR}/outdated_summary.md"

NOW_UTC="$(date -u +"%Y-%m-%dT%H:%M:%SZ")"

{
  echo "# Dependency Audit Report"
  echo
  echo "- Generated: ${NOW_UTC}"
  echo "- Working dir: ${ROOT_DIR}"
  echo
} > "${OUTDATED_SUMMARY_MD}"

run_or_report() {
  local title="$1"
  shift
  set +e
  local out
  out="$("$@" 2>&1)"
  local code=$?
  set -e

  if [[ $code -ne 0 ]]; then
    {
      echo "## ${title}"
      echo
      echo "**Status:** FAILED (exit ${code})"
      echo
      echo '```'
      echo "${out}"
      echo '```'
      echo
    } >> "${OUTDATED_SUMMARY_MD}"
    return $code
  fi
  return 0
}

# Run a command, capture output, and write it to a file. Optionally allow exit code 1.
run_to_file() {
  local title="$1"
  local outfile="$2"
  local allow_code1="${3:-false}"
  shift 3 || true

  set +e
  local out
  out="$("$@" 2>&1)"
  local code=$?
  set -e

  printf "%s\n" "${out}" > "${outfile}"

  if [[ "${allow_code1}" == "true" ]]; then
    if [[ $code -eq 0 || $code -eq 1 ]]; then
      return 0
    fi
  fi

  if [[ $code -ne 0 ]]; then
    {
      echo "## ${title}"
      echo
      echo "**Status:** FAILED (exit ${code})"
      echo
      echo '```'
      echo "${out}"
      echo '```'
      echo
    } >> "${OUTDATED_SUMMARY_MD}"
    return $code
  fi
  return 0
}

# Security audit (pip-audit)
if command -v pip-audit >/dev/null 2>&1; then
  # pip-audit exits 1 if vulnerabilities are found; treat that as a successful run.
  set +e
  pa_out="$(pip-audit -f json -o "${SEC_JSON}" 2>&1)"
  pa_code=$?
  set -e
  if [[ $pa_code -ne 0 && $pa_code -ne 1 ]]; then
    {
      echo "## pip-audit (json)"
      echo
      echo "**Status:** FAILED (exit ${pa_code})"
      echo
      echo '```'
      echo "${pa_out}"
      echo '```'
      echo
    } >> "${OUTDATED_SUMMARY_MD}"
    echo "{}" > "${SEC_JSON}"
  fi

  set +e
  pa_md_out="$(pip-audit -f markdown -o "${SEC_MD}" 2>&1)"
  pa_md_code=$?
  set -e
  if [[ $pa_md_code -ne 0 && $pa_md_code -ne 1 ]]; then
    {
      echo "## pip-audit (markdown)"
      echo
      echo "**Status:** FAILED (exit ${pa_md_code})"
      echo
      echo '```'
      echo "${pa_md_out}"
      echo '```'
      echo
    } >> "${OUTDATED_SUMMARY_MD}"
    {
      echo "# Security Audit"
      echo
      echo "_pip-audit failed to generate markdown output; see summary for details._"
    } > "${SEC_MD}"
  fi
else
  {
    echo "## pip-audit"
    echo
    echo "**Status:** SKIPPED (pip-audit not found on PATH)"
    echo
  } >> "${OUTDATED_SUMMARY_MD}"
  echo "{}" > "${SEC_JSON}"
  {
    echo "# Security Audit"
    echo
    echo "_pip-audit not found on PATH; install dependencies then rerun._"
  } > "${SEC_MD}"
fi

# Outdated packages
run_to_file "pip list --outdated (json)" "${OUTDATED_JSON}" false pip list --outdated --format json || true

# Summaries (best-effort; never fail the script on parsing)
ROOT_DIR="${ROOT_DIR}" python - <<'PY' || true
import json
import os
from pathlib import Path

root = Path(os.environ.get("ROOT_DIR", ".")).resolve()
reports = root / "reports"

sec_json = reports / "security_audit.json"
out_json = reports / "outdated.json"
summary_md = reports / "outdated_summary.md"

def _read_json(p: Path):
    try:
        return json.loads(p.read_text())
    except Exception:
        return None

sec = _read_json(sec_json)
out = _read_json(out_json)

lines = []
lines.append("## Summary")
lines.append("")

vuln_count = None
if isinstance(sec, list):
    vuln_count = sum(len(item.get("vulns") or []) for item in sec)
elif isinstance(sec, dict) and "dependencies" in sec:
    deps = sec.get("dependencies") or []
    vuln_count = sum(len((d.get("vulns") or [])) for d in deps)

if vuln_count is None:
    lines.append("- Security audit: unable to parse JSON output.")
else:
    lines.append(f"- Security audit: **{vuln_count}** vulnerabilities reported (see `reports/security_audit.md`).")

if isinstance(out, list):
    lines.append(f"- Outdated packages: **{len(out)}** (see details below).")
else:
    lines.append("- Outdated packages: unable to parse JSON output.")

lines.append("")

if isinstance(out, list) and out:
    lines.append("## Outdated Packages (top 50)")
    lines.append("")
    lines.append("| Package | Current | Latest | Type |")
    lines.append("|---|---:|---:|---|")
    for item in out[:50]:
        name = item.get("name", "")
        cur = item.get("version", "")
        latest = item.get("latest_version", "")
        typ = item.get("latest_filetype", "")
        lines.append(f"| {name} | {cur} | {latest} | {typ} |")
    lines.append("")

existing = summary_md.read_text().splitlines() if summary_md.exists() else []
with summary_md.open("w", encoding="utf-8") as f:
    for ln in existing:
        f.write(ln + "\n")
    f.write("\n")
    for ln in lines:
        f.write(ln + "\n")
PY

echo "Wrote:"
echo "  - ${SEC_MD}"
echo "  - ${OUTDATED_SUMMARY_MD}"


