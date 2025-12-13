#!/usr/bin/env bash
set -euo pipefail

# -----------------------------------------------------------------------------
# Amprenta - Lightsail bootstrap (Ubuntu)
#
# What this script does:
# - Installs basic system dependencies
# - Installs Miniconda
# - Clones the repo (if needed)
# - Creates and activates a conda environment
# - Installs Python dependencies from requirements.txt
#
# What this script does NOT do:
# - It does NOT set secrets or production env vars
# - It does NOT open AWS firewall ports
# - It does NOT create systemd services (left to README)
# -----------------------------------------------------------------------------

REPO_DIR="${REPO_DIR:-$HOME/amprenta}"
REPO_URL="${REPO_URL:-}"
CONDA_DIR="${CONDA_DIR:-$HOME/miniconda3}"
CONDA_ENV_NAME="${CONDA_ENV_NAME:-amprenta}"
PYTHON_VERSION="${PYTHON_VERSION:-3.12}"

echo "[1/6] Installing system dependencies..."
sudo apt-get update -y
sudo apt-get install -y \
  git curl wget ca-certificates \
  build-essential \
  unzip \
  python3 python3-pip

echo "[2/6] Installing Miniconda (if needed)..."
if [ ! -d "$CONDA_DIR" ]; then
  MINICONDA_SH="/tmp/miniconda.sh"
  wget -qO "$MINICONDA_SH" "https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh"
  bash "$MINICONDA_SH" -b -p "$CONDA_DIR"
  rm -f "$MINICONDA_SH"
fi

# shellcheck disable=SC1090
source "$CONDA_DIR/etc/profile.d/conda.sh"

echo "[3/6] Cloning repository (if needed)..."
if [ ! -d "$REPO_DIR/.git" ]; then
  if [ -z "$REPO_URL" ]; then
    echo "ERROR: REPO_URL is not set. Run:"
    echo "  REPO_URL=<your_repo_url> bash deploy/aws/setup-lightsail.sh"
    exit 1
  fi
  git clone "$REPO_URL" "$REPO_DIR"
fi

cd "$REPO_DIR"

echo "[4/6] Creating conda environment (if needed)..."
if ! conda env list | awk '{print $1}' | grep -qx "$CONDA_ENV_NAME"; then
  conda create -y -n "$CONDA_ENV_NAME" "python=$PYTHON_VERSION"
fi

echo "[5/6] Installing Python dependencies..."
conda activate "$CONDA_ENV_NAME"
python -m pip install --upgrade pip
python -m pip install -r requirements.txt

echo "[6/6] Done."
echo
echo "NEXT STEPS:"
echo "1) Create a .env file in the repo root with your secrets and RDS settings."
echo "2) Run migrations:"
echo "   conda activate $CONDA_ENV_NAME && alembic upgrade head"
echo "3) Start Streamlit:"
echo "   conda activate $CONDA_ENV_NAME && streamlit run streamlit_app/app.py --server.port 8501 --server.address 0.0.0.0"
echo
echo "Placeholders you must fill (do NOT commit secrets):"
echo "  - POSTGRES_HOST / POSTGRES_USER / POSTGRES_PASSWORD / POSTGRES_DB"
echo "  - OPENAI_API_KEY / PINECONE_API_KEY"


