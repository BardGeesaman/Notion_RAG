# Workstation Setup Guide

Quick setup guide for running the Amprenta RAG Platform on any workstation.

---

## Quick Setup (Recommended)

From your cloned repository:

```bash
cd /path/to/RAG
bash scripts/setup_workstation.sh
```

This script:

- Detects your username (`whoami`)
- Creates `.env` from `.env.template` (if it does not exist)
- Replaces `YOUR_USERNAME_HERE` with your actual username
- Installs Python dependencies via `pip install -r requirements.txt`
- Installs Playwright Chromium browser via `playwright install chromium`
- Clears Python caches (`__pycache__`, `*.pyc`)
- Creates the `logs/` directory

> Note: The script **runs real install commands** (`pip`, `playwright`). Review it before running and ensure you are in the correct virtual environment.

---

## Manual Setup

If you prefer manual control or need to debug issues, follow these steps.

### 1. Install Python Dependencies

From the repo root:

```bash
pip install -r requirements.txt
playwright install chromium
```

### 2. Configure Environment

Create your `.env` from the template:

```bash
cp .env.template .env
# Then edit .env and set:
#   POSTGRES_USER=$(whoami)
#   POSTGRES_PASSWORD=your_database_password
#   OPENAI_API_KEY=...
#   PINECONE_API_KEY=...
```

### 3. Initialize Database

Run Postgres migrations:

```bash
alembic upgrade head
```

### 4. Verify Setup

Basic checks:

```bash
# Test core dashboard imports
pytest amprenta_rag/tests/dashboard/test_pages_import.py

# Start dashboard
streamlit run scripts/run_dashboard.py

# Start API
uvicorn amprenta_rag.api.main:app --reload
```

---

## Troubleshooting & Workstation-Specific Issues

### Username Mismatch (Postgres Role Errors)

Symptom:
- Errors like `FATAL: role "bard" does not exist` or similar.

Fix:

- Check your OS username:

```bash
whoami
```

- Update `.env`:

```env
POSTGRES_USER=your_actual_username
```

### Missing Packages

Symptom:
- `ModuleNotFoundError` for libraries like `streamlit`, `playwright`, `psycopg2`, `fastapi`, etc.

Fix:

```bash
pip install streamlit playwright psycopg2-binary scikit-learn fastapi uvicorn
```

If you are using a virtual environment, ensure it is activated before running `pip install`.

### Python Cache Issues

Symptom:
- Tests or imports behave inconsistently after code changes.

Fix:

```bash
find . -type d -name __pycache__ -exec rm -rf {} + 2>/dev/null
find . -name "*.pyc" -delete
```

---

## Switching Between Workstations

Typical flow for moving from Workstation A to Workstation B:

### 1. On Workstation A

Commit and push your latest changes:

```bash
cd /path/to/RAG
git add -A
git commit -m "Your changes"
git push origin main
```

### 2. On Workstation B

Pull and run setup:

```bash
cd /path/to/RAG
git pull origin main
bash scripts/setup_workstation.sh
```

### 3. Verify

```bash
pytest amprenta_rag/tests/dashboard/test_pages_import.py
streamlit run scripts/run_dashboard.py
```

---

## Environment Variables Checklist

For full functionality, ensure these are set (typically via `.env`):

- **POSTGRES_HOST**: e.g., `localhost`
- **POSTGRES_PORT**: e.g., `5432`
- **POSTGRES_DB**: `amprenta` (canonical local database name)
- **POSTGRES_USER**: your OS username (or DB role)
- **POSTGRES_PASSWORD**: password for the Postgres user

- **OPENAI_API_KEY**: OpenAI API key
- **PINECONE_API_KEY**: Pinecone API key
- **PINECONE_INDEX_NAME**: e.g., `amprenta-rag`
- **PINECONE_ENVIRONMENT**: Pinecone environment/region

Optional (have defaults in code, but can be overridden):

- **FEATURE_CACHE_ENABLED**
- **FEATURE_CACHE_TTL**
- **FEATURE_CACHE_MAX_SIZE**
- **AUTO_LINK_ENABLED**
- **AUTO_LINK_MIN_CONFIDENCE**

Keep `.env.template` under version control and `.env` out of version control so each workstation can safely customize credentials without leaking secrets.


