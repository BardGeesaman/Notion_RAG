# Local Development Setup - Amprenta Multi-Omics Platform

**Status**: Production-ready local development environment

**Purpose**: Guide for setting up the local Postgres-first development environment with Streamlit UI and multi-agent workflows.

---

## Overview

This local workspace uses:
- **Postgres** as system of record for all multi-omics data
- **Docker Compose** for local Postgres instance
- **Alembic** for database migrations
- **Streamlit** for interactive UI
- **FastAPI** for REST API (optional)
- **SQLite** for chemistry/HTS staging only

**Notion is NO LONGER part of the RAG loop** - it's an optional documentation layer only.

---

## Prerequisites

### Required Software

1. **Python 3.9+**
   ```bash
   python --version  # Should be 3.9 or higher
   ```

2. **Docker & Docker Compose**
   - Install Docker Desktop: https://www.docker.com/products/docker-desktop
   - Verify installation:
   ```bash
   docker --version
   docker-compose --version
   ```

3. **Git** (already installed if you cloned the repo)

### Optional Tools

- **pgAdmin** (included in Docker Compose, accessible at http://localhost:5050)
- **psql** command-line client (for advanced DB operations)

---

## Step 1: Clone Repository

If you haven't already:

```bash
git clone <your-repo-url> amprenta_rag_local
cd amprenta_rag_local
```

---

## Step 2: Start Postgres with Docker Compose

### Start Services

```bash
# Start Postgres + pgAdmin in detached mode
docker-compose up -d

# Verify services are running
docker-compose ps
```

Expected output:
```
NAME                 IMAGE                    STATUS
amprenta-postgres    postgres:15-alpine       Up (healthy)
amprenta-pgadmin     dpage/pgadmin4:latest    Up
```

### Access Postgres

- **Host**: `localhost`
- **Port**: `5432`
- **Database**: `amprenta_rag`
- **User**: `postgres`
- **Password**: `postgres`

### Access pgAdmin (Optional)

1. Open browser: http://localhost:5050
2. Login:
   - Email: `admin@amprenta.local`
   - Password: `admin`
3. Add server connection:
   - Host: `postgres` (Docker internal network)
   - Port: `5432`
   - Database: `amprenta_rag`
   - Username: `postgres`
   - Password: `postgres`

### Stop Services

```bash
# Stop services (preserves data)
docker-compose stop

# Stop and remove containers (preserves data in volumes)
docker-compose down

# Stop and remove everything including data (DESTRUCTIVE)
docker-compose down -v
```

---

## Step 3: Python Environment Setup

### Create Virtual Environment

```bash
# Create virtual environment
python -m venv venv

# Activate virtual environment
# On macOS/Linux:
source venv/bin/activate

# On Windows:
venv\Scripts\activate
```

### Install Dependencies

```bash
# Install core dependencies
pip install -r requirements.txt

# Install Streamlit dependencies (if building UI)
pip install streamlit plotly pandas
```

### Optional: Structural Biology Dependencies

**Note**: Structural biology features (protein structure preparation, molecular docking) are optional and require additional dependencies that may be challenging to install via pip on some platforms.

#### Recommended Installation (conda-forge)

```bash
# Install via conda-forge (most reliable)
conda install -c conda-forge pdbfixer openmm
```

#### Alternative Installation (pip)

```bash
# Install via pip (may fail on some platforms)
pip install -r requirements-structural.txt
```

#### What requires structural dependencies:

- `amprenta_rag.structural.prep` - Protein structure preparation (PDBFixer)
- Molecular docking workflows
- Structure-based drug design features

#### Troubleshooting

If structural biology dependencies fail to install:
1. The core platform will work without them
2. Structural features will show helpful error messages
3. Consider using conda-forge instead of pip
4. On Apple Silicon (M1/M2), conda-forge is strongly recommended

### Optional: OCR Support

For extracting text from scanned PDFs and images:

**Linux:**
```bash
apt-get install tesseract-ocr poppler-utils
```

**macOS:**
```bash
brew install tesseract poppler
```

**Windows:**
- Tesseract: Download from https://github.com/UB-Mannheim/tesseract/wiki
- Poppler: Download from https://github.com/oschwartz10612/poppler-windows/releases

**Note**: OCR functionality gracefully degrades if tesseract is not installed - the platform will continue to work but scanned PDFs will not be processed.

---

## Step 4: Environment Configuration

### Copy .env.example to .env

```bash
cp .env.example .env
```

### Edit .env File

Open `.env` in your editor and configure:

#### Required: Postgres Configuration
```env
# Postgres Database (Local Docker)
POSTGRES_HOST=localhost
POSTGRES_PORT=5432
POSTGRES_DB=amprenta
POSTGRES_USER=postgres
POSTGRES_PASSWORD=postgres
POSTGRES_ECHO=false

# Use Postgres as System of Record (REQUIRED)
USE_POSTGRES_AS_SOT=true
```

#### Required: API Keys
```env
# OpenAI (for embeddings and reasoning)
OPENAI_API_KEY=sk-...

# Pinecone (for vector storage)
PINECONE_API_KEY=...
```

#### Optional: Notion (for legacy migration only)
```env
# Notion Integration (OPTIONAL - only needed for migration)
ENABLE_NOTION_SYNC=false
NOTION_API_KEY=secret_...

# Notion Database IDs (only if you're migrating data)
NOTION_EXP_DATA_DB_ID=...
NOTION_SIGNATURE_DB_ID=...
# ... other database IDs
```

#### Optional: Additional Services
```env
# Zotero (if using literature ingestion)
ZOTERO_API_KEY=...

# Feature Caching
FEATURE_CACHE_TTL_SECONDS=3600
FEATURE_CACHE_ENABLE_PERSISTENCE=true
```

---

## Step 5: Database Migrations

### Apply Alembic Migrations

```bash
# Run all pending migrations
alembic upgrade head
```

Expected output:
```
INFO  [alembic.runtime.migration] Running upgrade  -> 80008f433a68, Initial schema: programs, experiments, datasets, features, signatures
INFO  [alembic.runtime.migration] Running upgrade 80008f433a68 -> 0c9c72e35979, add missing metadata fields to dataset
INFO  [alembic.runtime.migration] Running upgrade 0c9c72e35979 -> d2d66c54e753, add signature metadata fields
```

### Verify Database Schema

```bash
# Validate Postgres connection and schema
python scripts/validate_postgres.py
```

Expected output:
```
✓ Postgres connection successful
✓ Database 'amprenta_rag' exists
✓ All tables created successfully:
  - programs
  - experiments
  - datasets
  - features
  - signatures
  - signature_components
  - (junction tables)
✓ Schema validation complete
```

---

## Step 6: Verify Installation

### Run Health Check

```bash
python scripts/validate_postgres.py
```

### Test Ingestion (Optional)

```bash
# Test ingesting a sample dataset (if you have test data)
python scripts/ingest_lipidomics.py --file data/test_lipidomics.csv --create-page
```

### Start Streamlit UI (Once Phase 4 is complete)

```bash
streamlit run streamlit_app/app.py
```

Access at: http://localhost:8501

---

## Common Operations

### Database Management

#### View Alembic Migration History
```bash
alembic history
alembic current
```

#### Rollback Migration
```bash
alembic downgrade -1  # Rollback one migration
alembic downgrade <revision>  # Rollback to specific revision
```

#### Create New Migration (Advanced)
```bash
alembic revision --autogenerate -m "description of changes"
alembic upgrade head
```

#### Reset Database (DESTRUCTIVE)
```bash
# Drop all tables and re-create
alembic downgrade base
alembic upgrade head
```

### Docker Operations

#### View Logs
```bash
# All services
docker-compose logs -f

# Postgres only
docker-compose logs -f postgres
```

#### Restart Services
```bash
docker-compose restart
```

#### Rebuild Containers
```bash
docker-compose down
docker-compose up -d --build
```

#### Backup Database
```bash
docker exec amprenta-postgres pg_dump -U postgres amprenta_rag > backup.sql
```

#### Restore Database
```bash
cat backup.sql | docker exec -i amprenta-postgres psql -U postgres amprenta_rag
```

---

## Troubleshooting

### Issue: Docker Compose Fails to Start

**Symptom**: `Cannot start service postgres: port is already allocated`

**Solution**:
```bash
# Check if another Postgres instance is running
lsof -i :5432

# Stop the conflicting service or change port in docker-compose.yml
```

### Issue: Alembic Migrations Fail

**Symptom**: `Target database is not up to date`

**Solution**:
```bash
# Check current migration
alembic current

# View migration history
alembic history

# Force upgrade to latest
alembic upgrade head
```

### Issue: Cannot Connect to Postgres

**Symptom**: `could not connect to server: Connection refused`

**Solution**:
1. Verify Docker containers are running: `docker-compose ps`
2. Check Postgres logs: `docker-compose logs postgres`
3. Verify `.env` has correct connection settings
4. Try restarting services: `docker-compose restart`

### Issue: Streamlit Not Finding Database

**Symptom**: Streamlit shows database connection errors

**Solution**:
1. Ensure virtual environment is activated
2. Check `.env` file is in repo root
3. Verify `USE_POSTGRES_AS_SOT=true` in `.env`
4. Run `python scripts/validate_postgres.py` to verify connection

---

## Architecture Notes

### Postgres as System of Record

- **All multi-omics data** (programs, experiments, datasets, features, signatures) lives in Postgres
- **Ingestion pipelines** write directly to Postgres (not Notion)
- **RAG text generation** queries Postgres (10-100x faster than Notion API)
- **Streamlit UI** connects to Postgres for all data

### Notion (Optional)

- **Notion is NO LONGER required** for RAG or ingestion
- If you have existing Notion data, use migration scripts in Phase 5
- Notion can be used as a documentation layer (read-only) if desired
- Set `ENABLE_NOTION_SYNC=false` to disable Notion entirely

### SQLite (Chemistry Only)

- SQLite is used ONLY for chemistry/HTS prototyping (`data/chemistry/chemistry.db`)
- Multi-omics data does NOT use SQLite
- Chemistry data can be promoted to Postgres later via scripts

### Multi-Agent Workflows

- Use agent prompts in `agents/` folder for efficient development:
  - `@architect.md` - Design phase
  - `@implementor.md` - Implementation phase
  - `@reviewer.md` - Code review phase
  - `@testing.md` - Test creation phase
- See `docs/AGENT_WORKFLOW.md` for detailed usage (Phase 6)

---

## Next Steps

### For New Users

1. ✅ Complete Steps 1-5 above
2. ✅ Run `python scripts/validate_postgres.py` to verify setup
3. 📖 Read `context/MASTER_CONTEXT_FOR_NEW_CHAT.md` for system overview
4. 🚀 Start ingesting data or exploring with Streamlit

### For Notion Migration

1. Complete local setup (Steps 1-5)
2. Run `python scripts/export_notion_to_json.py` to export Notion data
3. Run `python scripts/import_notion_json_to_postgres.py` to import to Postgres
4. Validate with `python scripts/validate_migration.py`
5. See Phase 5 documentation for details

### For Development

1. Review multi-agent workflow: `docs/AGENT_WORKFLOW.md` (once created)
2. Follow coding patterns in existing modules
3. Use Alembic for schema changes
4. Add tests for new features
5. Update context documents when making significant changes

---

## Quick Reference

### Start Development Session

```bash
# 1. Start Postgres
docker-compose up -d

# 2. Activate virtual environment
source venv/bin/activate  # macOS/Linux
# venv\Scripts\activate  # Windows

# 3. Verify database
python scripts/validate_postgres.py

# 4. Start Streamlit (optional)
streamlit run streamlit_app/app.py
```

### Stop Development Session

```bash
# Stop Postgres (preserves data)
docker-compose stop

# Deactivate virtual environment
deactivate
```

---

**Last Updated**: 2025-12-06

**For Questions**: See `context/MASTER_CONTEXT_FOR_NEW_CHAT.md` or create an issue in the repository.

