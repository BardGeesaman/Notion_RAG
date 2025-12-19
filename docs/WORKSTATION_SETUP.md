# Workstation Setup Guide

Complete guide for setting up the Amprenta RAG development environment on a new workstation.

## Prerequisites

### macOS (Homebrew)

```bash
# Install Homebrew if not present
/bin/bash -c "$(curl -fsSL https://raw.githubusercontent.com/Homebrew/install/HEAD/install.sh)"

# Install required tools
brew install git
brew install --cask docker
brew install --cask miniconda

# Start Docker Desktop
open -a Docker
```

### Cursor IDE

1. Download from https://cursor.sh
2. Install and open
3. Sign in to sync settings (optional)

---

## 1. Install Cursor Extensions

Install these extensions via Command Palette (`Cmd+Shift+P` → "Install Extension"):

| Extension | Extension ID |
|-----------|--------------|
| Python | ms-python.python |
| Python Debugger | ms-python.debugpy |
| Pylance (Type Checking) | ms-python.vscode-pylance |
| Jupyter Notebook | ms-toolsai.jupyter |
| Docker | ms-azuretools.vscode-docker |
| Dev Containers | ms-vscode-remote.remote-containers |
| GitLens | eamodio.gitlens |
| GitHub Pull Requests | github.vscode-pull-request-github |
| GitHub Actions | github.vscode-github-actions |
| Git History | donjayamanne.githistory |
| markdownlint | davidanson.vscode-markdownlint |
| Terraform | hashicorp.terraform |
| Kubernetes | ms-kubernetes-tools.vscode-kubernetes-tools |
| Vim | vscodevim.vim |

Or install via terminal:
```bash
cursor --install-extension ms-python.python
cursor --install-extension ms-python.debugpy
cursor --install-extension ms-python.vscode-pylance
cursor --install-extension ms-toolsai.jupyter
cursor --install-extension ms-azuretools.vscode-docker
cursor --install-extension ms-vscode-remote.remote-containers
cursor --install-extension eamodio.gitlens
cursor --install-extension github.vscode-pull-request-github
cursor --install-extension github.vscode-github-actions
cursor --install-extension donjayamanne.githistory
cursor --install-extension davidanson.vscode-markdownlint
cursor --install-extension hashicorp.terraform
cursor --install-extension ms-kubernetes-tools.vscode-kubernetes-tools
cursor --install-extension vscodevim.vim
```

---

## 2. Clone Repository

```bash
cd ~/Documents
git clone https://github.com/BardGeesaman/Notion_RAG.git RAG
cd RAG
```

---

## 3. Create Conda Environment

### Option A: Exact Replication (Recommended)

Use the exported environment file for exact package versions:

```bash
# Initialize conda for your shell (one-time)
conda init zsh
source ~/.zshrc

# Create environment from export (exact versions)
conda env create -f environment.yml
conda activate myenv

# Install Playwright browsers (for E2E tests)
python -m playwright install
```

### Option B: Fresh Install

Create a new environment and install from requirements:

```bash
# Initialize conda for your shell (one-time)
conda init zsh
source ~/.zshrc

# Create and activate environment
conda create -n myenv python=3.12 -y
conda activate myenv

# Install Python dependencies
python -m pip install -r requirements.txt

# Install Playwright browsers (for E2E tests)
python -m playwright install
```

**Note:** Always use `python -m pip install` instead of `pip install` to ensure packages install in the conda environment.

---

## 4. Configure Environment Variables

```bash
# Copy template
cp .env.example .env

# Edit with your API keys
# Required: OPENAI_API_KEY, PINECONE_API_KEY, ZOTERO_API_KEY
```

Minimum required `.env` contents:
```bash
OPENAI_API_KEY=sk-...
PINECONE_API_KEY=...
ZOTERO_API_KEY=...

POSTGRES_HOST=localhost
POSTGRES_PORT=5432
POSTGRES_DB=amprenta
POSTGRES_USER=postgres
POSTGRES_PASSWORD=postgres
```

---

## 5. Start Docker Services

```bash
# Start Postgres and pgAdmin
docker compose up -d

# Verify containers are running
docker ps
# Expected: amprenta-postgres, amprenta-pgadmin
```

**Access pgAdmin:** http://localhost:5050
- Email: admin@amprenta.local
- Password: admin

---

## 6. Initialize Database

```bash
# Activate conda environment
conda activate myenv

# Run database migrations
alembic upgrade head

# (Optional) Seed test data
python scripts/seed_sar_data.py
python scripts/seed_hts_data.py
```

---

## 7. Verify Setup

```bash
# Validate configuration
python scripts/validate_configuration.py

# Run unit tests
pytest tests/ -v --ignore=tests/e2e

# Start API server (terminal 1)
uvicorn amprenta_rag.api.main:app --reload --port 8000

# Start dashboard (terminal 2)
cd scripts/dashboard && streamlit run app.py
```

---

## 8. (Optional) Start JupyterHub

For Voila dashboards and Jupyter notebooks:

```bash
cd deploy/jupyterhub
docker compose up -d
```

Access at: http://localhost:8000

---

## Quick Reference

| Service | URL | Notes |
|---------|-----|-------|
| FastAPI | http://localhost:8000 | API server |
| API Docs | http://localhost:8000/docs | Swagger UI |
| Streamlit | http://localhost:8501 | Dashboard |
| pgAdmin | http://localhost:5050 | Database admin |
| JupyterHub | http://localhost:8000 | When running |

---

## Troubleshooting

### conda: command not found
```bash
source ~/.zshrc
# or restart terminal
```

### Docker permission denied
- Ensure Docker Desktop is running
- Check Docker Desktop → Settings → General → "Start Docker Desktop when you log in"

### Database connection failed
```bash
# Check if Postgres is running
docker ps | grep postgres

# Restart if needed
docker compose down && docker compose up -d
```

### Missing API keys error
- Verify `.env` file exists in project root
- Check all required keys are set (OPENAI_API_KEY, PINECONE_API_KEY, ZOTERO_API_KEY)

### pip installing to wrong Python
Always use:
```bash
python -m pip install <package>
```
Not:
```bash
pip install <package>  # May use system Python!
```

---

## Daily Workflow

```bash
# Start of day
cd ~/Documents/RAG
conda activate myenv
docker compose up -d

# Start services
uvicorn amprenta_rag.api.main:app --reload --port 8000 &
cd scripts/dashboard && streamlit run app.py &

# End of day
docker compose down
```
