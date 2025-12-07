#!/usr/bin/env bash

# Start the FastAPI server (uvicorn) for the Amprenta API.
#
# Defaults:
#   - Host: 0.0.0.0
#   - Port: 8000 (override with API_PORT env)
#   - Reload: enabled by default for development (disable with API_RELOAD=0)
#   - Logs: appended to logs/api.log
#
# Usage:
#   bash scripts/start_api.sh
#   API_PORT=9000 bash scripts/start_api.sh
#   API_RELOAD=0 bash scripts/start_api.sh  # disable auto-reload (closer to prod)

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]:-$0}")" && pwd)"
REPO_ROOT="$(cd "${SCRIPT_DIR}/.." && pwd)"
LOG_DIR="${REPO_ROOT}/logs"
LOG_FILE="${LOG_DIR}/api.log"

HOST="${API_HOST:-0.0.0.0}"
PORT="${API_PORT:-8000}"
RELOAD_FLAG="${API_RELOAD:-1}"

PYTHON_BIN="${PYTHON_BIN:-python3}"

mkdir -p "${LOG_DIR}"

timestamp() {
  date +"%Y-%m-%d %H:%M:%S"
}

echo "[$(timestamp)] [INFO] Starting FastAPI (uvicorn) on ${HOST}:${PORT} (reload=${RELOAD_FLAG})" | tee -a "${LOG_FILE}"
cd "${REPO_ROOT}"

UVICORN_CMD=(
  "${PYTHON_BIN}" -m uvicorn
  "amprenta_rag.api.main:app"
  "--host" "${HOST}"
  "--port" "${PORT}"
)

if [[ "${RELOAD_FLAG}" == "1" ]]; then
  UVICORN_CMD+=("--reload")
fi

echo "[$(timestamp)] [INFO] Command: ${UVICORN_CMD[*]}" | tee -a "${LOG_FILE}"

# Run in the foreground; caller can manage process (systemd, tmux, etc.).
"${UVICORN_CMD[@]}" >> "${LOG_FILE}" 2>&1


