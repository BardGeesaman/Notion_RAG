# Presentation Layer Guide - Postgres Migration

**Status**: FastAPI API Ready, Frontend Pending  
**Last Updated**: 2025-01-XX

## Current State

With Postgres as primary and Notion dual-write disabled:

- ✅ **Postgres**: Primary database (source of truth)
- ✅ **FastAPI**: REST API layer (built and ready)
- ⏳ **Frontend UI**: Not yet built (pending)
- ❌ **Notion**: Disabled (was presentation layer)

## Available Presentation Layers

### 1. FastAPI REST API (Current - Built)

**Status**: ✅ Fully operational

**Access**:
- API Server: `http://localhost:8000`
- Interactive Docs: `http://localhost:8000/docs` (Swagger UI)
- Alternative Docs: `http://localhost:8000/redoc` (ReDoc)

**Start Server**:
```bash
python scripts/run_api_server.py
```

**Available Endpoints**:
- `GET /api/v1/programs` - List all programs
- `POST /api/v1/programs` - Create program
- `GET /api/v1/programs/{id}` - Get program by ID
- `PATCH /api/v1/programs/{id}` - Update program
- `DELETE /api/v1/programs/{id}` - Delete program

- `GET /api/v1/experiments` - List experiments
- `POST /api/v1/experiments` - Create experiment
- `GET /api/v1/experiments/{id}` - Get experiment

- `GET /api/v1/datasets` - List datasets
- `POST /api/v1/datasets` - Create dataset
- `GET /api/v1/datasets/{id}` - Get dataset

- `GET /api/v1/features` - List features
- `GET /api/v1/signatures` - List signatures

**Example Usage**:
```bash
# List all datasets
curl http://localhost:8000/api/v1/datasets

# Create a program
curl -X POST http://localhost:8000/api/v1/programs \
  -H "Content-Type: application/json" \
  -d '{"name": "ALS Program", "description": "ALS research"}'

# Get dataset by ID
curl http://localhost:8000/api/v1/datasets/{uuid}
```

### 2. Interactive API Documentation (Swagger UI)

**Status**: ✅ Available

**Access**: `http://localhost:8000/docs`

**Features**:
- Browse all endpoints
- Test API calls directly in browser
- See request/response schemas
- Try out operations without writing code

### 3. Direct Database Access (Postgres)

**Status**: ✅ Available

**Access**: Direct SQL queries

```bash
# Connect to database
psql amprenta_rag

# Query datasets
SELECT name, omics_type, created_at FROM datasets;

# Query with relationships
SELECT d.name, p.name as program_name
FROM datasets d
JOIN dataset_program dp ON d.id = dp.dataset_id
JOIN programs p ON dp.program_id = p.id;
```

### 4. Python API Client (Recommended)

**Status**: ✅ Available

**Usage**:
```python
import requests

# List datasets
response = requests.get("http://localhost:8000/api/v1/datasets")
datasets = response.json()

# Create dataset
new_dataset = {
    "name": "Test Dataset",
    "omics_type": "lipidomics",
    "description": "Test description"
}
response = requests.post(
    "http://localhost:8000/api/v1/datasets",
    json=new_dataset
)
```

## What's Missing: Frontend UI

### Current Gap

- ❌ **Web UI**: No browser-based interface yet
- ❌ **Dashboard**: No visual dashboard
- ❌ **Data Browser**: No GUI for exploring data
- ❌ **Visualizations**: No charts/graphs

### Planned Frontend (From Roadmap)

**Option 1: Next.js/React Frontend** (Recommended)
- Modern, fast, scalable
- Can consume FastAPI REST API
- Rich interactive UI
- Status: ⏳ Pending (Tier 3.1 Phase 6)

**Option 2: Streamlit Dashboard** (Quick Option)
- Python-based, easy to build
- Good for internal tools
- Can connect to Postgres directly
- Status: ⏳ Not started

**Option 3: Jupyter Notebooks** (Analysis)
- For data exploration
- Can use FastAPI or direct Postgres
- Status: ✅ Can be used now

## Recommended Workflow (Current)

### For Data Ingestion
```bash
# Use CLI scripts (creates in Postgres)
python scripts/ingest_lipidomics.py --file data.csv --create-page
```

### For Data Access
```bash
# Option 1: FastAPI (recommended)
python scripts/run_api_server.py
# Then visit http://localhost:8000/docs

# Option 2: Direct Postgres
psql amprenta_rag -c "SELECT * FROM datasets;"

# Option 3: Python script
python -c "
from amprenta_rag.database.base import get_db
from amprenta_rag.database.models import Dataset
db = next(get_db())
datasets = db.query(Dataset).all()
for d in datasets:
    print(f'{d.name}: {d.omics_type}')
"
```

### For Data Exploration
```python
# Jupyter notebook
import pandas as pd
from sqlalchemy import create_engine

engine = create_engine('postgresql://user:pass@localhost/amprenta_rag')
df = pd.read_sql('SELECT * FROM datasets', engine)
df.head()
```

## Next Steps

### Immediate (Use What's Built)
1. **Use FastAPI**: Start server, use Swagger UI
2. **Use Python**: Write scripts to query Postgres
3. **Use CLI**: Continue using ingestion scripts

### Short-term (Build Frontend)
1. **Option A**: Build Streamlit dashboard (quick)
2. **Option B**: Build Next.js frontend (production-ready)
3. **Option C**: Use existing tools (Jupyter, Postgres clients)

### Long-term (Full Stack)
1. Complete Next.js/React frontend
2. Add authentication
3. Add visualizations
4. Add real-time updates

## Migration Path

**Before (Notion Primary)**:
- Notion = Database + Presentation Layer
- All data visible in Notion
- Human-friendly interface

**Now (Postgres Primary)**:
- Postgres = Database
- FastAPI = API Layer
- ⏳ Frontend = Presentation Layer (to be built)

**Future (Full Stack)**:
- Postgres = Database
- FastAPI = API Layer
- Next.js = Frontend
- Notion = Optional documentation layer

## Summary

**Current Presentation Options**:
1. ✅ **FastAPI Swagger UI** - Best for testing/exploration
2. ✅ **Postgres CLI** - Best for direct queries
3. ✅ **Python Scripts** - Best for automation
4. ⏳ **Web Frontend** - To be built

**Recommendation**: Use FastAPI Swagger UI (`/docs`) for now, plan to build a proper frontend next.

