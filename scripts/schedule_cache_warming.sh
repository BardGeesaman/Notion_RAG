#!/usr/bin/env bash

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]:-$0}")" && pwd)"
REPO_ROOT="$(cd "${SCRIPT_DIR}/.." && pwd)"
LOG_FILE="${REPO_ROOT}/logs/cache_warming.log"
PYTHON_BIN="${PYTHON_BIN:-python3}"

mkdir -p "$(dirname "${LOG_FILE}")"

log() {
  local level="$1"
  shift
  local ts
  ts="$(date +"%Y-%m-%d %H:%M:%S")"
  echo "[${ts}] [${level}] $*" | tee -a "${LOG_FILE}"
}

log INFO "Starting cache warming for all datasets"

if ! "${PYTHON_BIN}" "${REPO_ROOT}/scripts/warm_feature_cache.py" --all-datasets >>"${LOG_FILE}" 2>&1; then
  status=$?
  log ERROR "warm_feature_cache.py failed with exit code ${status}"
  echo "ALERT: Cache warming failed (see ${LOG_FILE})" >&2
  exit "${status}"
fi

log INFO "Cache warming complete; running cache health check"

if ! "${PYTHON_BIN}" "${REPO_ROOT}/scripts/monitor_cache_health.py" >>"${LOG_FILE}" 2>&1; then
  status=$?
  log ERROR "Cache health check reported unhealthy state (exit=${status})"
  echo "ALERT: Cache health below threshold or eviction rate too high (see ${LOG_FILE})" >&2
  exit "${status}"
fi

log INFO "Cache warming and health check completed successfully"

# Example cron (run daily at 03:00):
# 0 3 * * * /usr/bin/env bash /path/to/repo/scripts/schedule_cache_warming.sh


