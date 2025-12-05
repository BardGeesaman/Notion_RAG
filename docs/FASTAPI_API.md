# FastAPI Service Layer

**Status**: Foundation Complete, Core Endpoints Implemented  
**Roadmap**: `context/UNIFIED_STRATEGIC_ROADMAP.md` Section 3.1, Phase 3

## Overview

The FastAPI service layer provides a REST API for accessing the multi-omics platform data through Postgres. This is part of TIER 3 architecture evolution, moving from Notion-as-database to Postgres-as-database.

## Current Status

### ✅ Implemented

- **FastAPI Application Structure**
  - Main application with CORS middleware
  - Health check endpoints
  - Router organization

- **Programs API** (Complete)
  - `POST /api/v1/programs` - Create program
  - `GET /api/v1/programs` - List programs
  - `GET /api/v1/programs/{id}` - Get program
  - `PATCH /api/v1/programs/{id}` - Update program
  - `DELETE /api/v1/programs/{id}` - Delete program

- **Experiments API** (Complete)
  - `POST /api/v1/experiments` - Create experiment
  - `GET /api/v1/experiments` - List experiments
  - `GET /api/v1/experiments/{id}` - Get experiment
  - `PATCH /api/v1/experiments/{id}` - Update experiment
  - `DELETE /api/v1/experiments/{id}` - Delete experiment

### 🚧 Partially Implemented (Stubs)

- **Datasets API** - Stub endpoints created, services TODO
- **Features API** - Stub endpoints created, services TODO
- **Signatures API** - Stub endpoints created, services TODO

## Running the API Server

### Development Server

```bash
python scripts/run_api_server.py
```

### With Custom Options

```bash
python scripts/run_api_server.py --host 0.0.0.0 --port 8000 --reload
```

### Using Uvicorn Directly

```bash
uvicorn amprenta_rag.api.main:app --host 127.0.0.1 --port 8000 --reload
```

## API Documentation

Once the server is running:

- **Swagger UI**: http://127.0.0.1:8000/docs
- **ReDoc**: http://127.0.0.1:8000/redoc

## Architecture

```
amprenta_rag/api/
├── main.py              # FastAPI application
├── schemas.py           # Pydantic request/response schemas
├── dependencies.py      # Database session dependencies
├── routers/            # API route handlers
│   ├── programs.py     # ✅ Complete
│   ├── experiments.py  # ✅ Complete
│   ├── datasets.py     # 🚧 Stub
│   ├── features.py     # 🚧 Stub
│   └── signatures.py   # 🚧 Stub
└── services/           # CRUD business logic
    ├── programs.py     # ✅ Complete
    └── experiments.py  # ✅ Complete
```

## Request/Response Examples

### Create a Program

```bash
curl -X POST "http://127.0.0.1:8000/api/v1/programs" \
  -H "Content-Type: application/json" \
  -d '{
    "name": "ALS Research Program",
    "description": "Research into ALS biomarkers",
    "disease": ["ALS"]
  }'
```

### List Programs

```bash
curl "http://127.0.0.1:8000/api/v1/programs?skip=0&limit=10"
```

### Create an Experiment

```bash
curl -X POST "http://127.0.0.1:8000/api/v1/experiments" \
  -H "Content-Type: application/json" \
  -d '{
    "name": "CSF Lipidomics Study",
    "type": "patient",
    "description": "Analysis of CSF lipid profiles",
    "matrix": ["CSF"],
    "disease": ["ALS"],
    "program_ids": ["<program-uuid>"]
  }'
```

## Next Steps

1. **Complete CRUD Services**
   - Implement datasets service
   - Implement features service
   - Implement signatures service

2. **Complete API Routers**
   - Complete datasets router
   - Complete features router
   - Complete signatures router

3. **Add Authentication**
   - Basic authentication
   - OIDC/OAuth2 (future)

4. **Add Advanced Features**
   - Pagination helpers
   - Filtering/searching
   - Bulk operations
   - Validation improvements

## Configuration

The API uses the same Postgres configuration as the database layer. See `amprenta_rag/config.py` and `.env` file for database connection settings.

## Error Handling

The API uses standard HTTP status codes:
- `200` - Success
- `201` - Created
- `404` - Not Found
- `500` - Internal Server Error

Error responses follow standard FastAPI/Pydantic error format.

## Testing

To test the API, you can:

1. Use the interactive Swagger UI at `/docs`
2. Use curl commands (examples above)
3. Write integration tests (TODO)

## Dependencies

- `fastapi>=0.109.0` - Web framework
- `uvicorn[standard]>=0.27.0` - ASGI server
- `pydantic>=2.5.0` - Data validation
- `sqlalchemy>=2.0.0` - Database ORM

