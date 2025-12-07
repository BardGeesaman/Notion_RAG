#!/usr/bin/env bash

# Start the Streamlit dashboard for the Amprenta Multi-Omics Platform.
#
# Defaults:
#   - Host: 0.0.0.0
#   - Port: 8501 (override with DASHBOARD_PORT env)
#   - Logs: appended to logs/dashboard.log
#
# Usage:
#   bash scripts/start_dashboard.sh
#   DASHBOARD_PORT=8601 bash scripts/start_dashboard.sh

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]:-$0}")" && pwd)"
REPO_ROOT="$(cd "${SCRIPT_DIR}/.." && pwd)"
LOG_DIR="${REPO_ROOT}/logs"
LOG_FILE="${LOG_DIR}/dashboard.log"

HOST="${DASHBOARD_HOST:-0.0.0.0}"
PORT="${DASHBOARD_PORT:-8501}"

PYTHON_BIN="${PYTHON_BIN:-python3}"

mkdir -p "${LOG_DIR}"

timestamp() {
  date +"%Y-%m-%d %H:%M:%S"
}

echo "[$(timestamp)] [INFO] Starting Streamlit dashboard on ${HOST}:${PORT}" | tee -a "${LOG_FILE}"
cd "${REPO_ROOT}"

STREAMLIT_CMD=(
  "${PYTHON_BIN}" -m streamlit run
  "scripts/run_dashboard.py"
  "--server.address" "${HOST}"
  "--server.port" "${PORT}"
)

echo "[$(timestamp)] [INFO] Command: ${STREAMLIT_CMD[*]}" | tee -a "${LOG_FILE}"

"${STREAMLIT_CMD[@]}" >> "${LOG_FILE}" 2>&1


