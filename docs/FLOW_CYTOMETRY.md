# Flow Cytometry / FACS Data Analysis

Comprehensive guide to flow cytometry data ingestion, analysis, and visualization in the Amprenta RAG platform.

## Overview

The flow cytometry module provides end-to-end support for Flow Cytometry Standard (FCS) file analysis, including:

- **FCS File Parsing**: Support for FCS 2.0, 3.0, and 3.1 formats
- **Data Transformation**: Logicle, arcsinh, and compensation transforms
- **Interactive Gating**: Polygon, rectangle, quadrant, and boolean gate creation
- **Population Analysis**: Automated statistics calculation and export
- **Dashboard Interface**: Web-based analysis and visualization
- **API Integration**: RESTful endpoints for programmatic access

## Features

### ✅ Implemented Features

| Component | Feature | Description |
|-----------|---------|-------------|
| **Parser** | FCS File Support | FCS 2.0, 3.0, 3.1 format parsing |
| **Parser** | Metadata Extraction | Cytometer settings, acquisition parameters |
| **Parser** | Event Data Loading | Efficient columnar data access |
| **Transforms** | Logicle Transform | Bi-exponential scaling for flow data |
| **Transforms** | Arcsinh Transform | Inverse hyperbolic sine scaling |
| **Transforms** | Compensation | Spectral overlap correction |
| **Transforms** | Subsampling | Performance optimization for large datasets |
| **Gating** | Polygon Gates | Multi-vertex region selection |
| **Gating** | Rectangle Gates | Rectangular region selection |
| **Gating** | Quadrant Gates | Four-region division |
| **Gating** | Boolean Logic | AND, OR, NOT gate combinations |
| **Analysis** | Population Stats | Event counts, percentages, medians, CVs |
| **Analysis** | Hierarchical Gating | Parent-child gate relationships |
| **Storage** | Parquet Format | Efficient columnar event storage |
| **Storage** | PostgreSQL Schema | Metadata and gate definitions |
| **API** | RESTful Endpoints | 9 endpoints for complete workflow |
| **Dashboard** | 4-Tab Interface | Upload, visualization, gating, statistics |

## Quick Start

### 1. Upload FCS File

```python
# Via API
import requests

with open("sample.fcs", "rb") as f:
    response = requests.post(
        "http://localhost:8000/api/v1/flow-cytometry/upload",
        files={"file": f}
    )

dataset = response.json()
print(f"Dataset ID: {dataset['flow_dataset_id']}")
```

### 2. Create Gates

```python
# Polygon gate for lymphocytes
gate_data = {
    "gate_name": "Lymphocytes",
    "gate_type": "polygon",
    "x_parameter": "FSC-A",
    "y_parameter": "SSC-A",
    "gate_definition": {
        "vertices": [
            [20000, 15000], [80000, 15000], 
            [120000, 60000], [40000, 100000]
        ]
    }
}

response = requests.post(
    f"http://localhost:8000/api/v1/flow-cytometry/datasets/{dataset_id}/gates",
    json=gate_data
)
```

### 3. Analyze Populations

```python
# Get population statistics
response = requests.get(
    f"http://localhost:8000/api/v1/flow-cytometry/datasets/{dataset_id}/populations"
)

populations = response.json()
for pop in populations:
    print(f"{pop['population_name']}: {pop['n_events']:,} events ({pop['pct_of_total']:.1f}%)")
```

## FCS File Requirements

### Supported Formats

| Version | Support | Notes |
|---------|---------|-------|
| FCS 2.0 | ✅ Full | Legacy format support |
| FCS 3.0 | ✅ Full | Most common format |
| FCS 3.1 | ✅ Full | Latest standard |

### File Size Limits

- **Maximum file size**: 500MB
- **Recommended size**: <100MB for optimal performance
- **Event count**: Up to 10M events supported
- **Parameter count**: Up to 50 parameters

### Required Metadata

FCS files must contain these standard keywords:

```
$DATATYPE: Data format (I, F, D)
$MODE: List mode (L) or histogram (H)
$NEXTDATA: Offset to next dataset
$TOT: Total number of events
$PAR: Number of parameters
$PnN: Parameter n name
$PnR: Parameter n range
```

## Data Transformations

### Logicle Transform

Bi-exponential transformation optimized for flow cytometry data:

```python
from amprenta_rag.flow_cytometry.transforms import logicle_transform, auto_logicle_params

# Auto-detect parameters
T, W, M, A = auto_logicle_params(data)

# Apply transform
transformed_data = logicle_transform(data, T=T, W=W, M=M, A=A)
```

**Parameters:**
- `T`: Top of scale value (typically 262144)
- `W`: Width of linear region (typically 0.5-2.0)
- `M`: Number of decades (typically 4.5)
- `A`: Additional decades of negative data (typically 0)

### Arcsinh Transform

Inverse hyperbolic sine transformation:

```python
from amprenta_rag.flow_cytometry.transforms import arcsinh_transform

# Apply arcsinh with cofactor
transformed_data = arcsinh_transform(data, cofactor=150)
```

**Use cases:**
- CyTOF (mass cytometry) data
- When logicle parameters are unavailable
- Simple alternative to logicle

### Compensation

Spectral overlap correction:

```python
from amprenta_rag.flow_cytometry.transforms import apply_compensation

# Compensation matrix (parameter x parameter)
comp_matrix = np.array([
    [1.0, 0.1, 0.05],  # FITC spillover
    [0.02, 1.0, 0.15], # PE spillover  
    [0.0, 0.08, 1.0]   # PerCP spillover
])

compensated_data = apply_compensation(data, comp_matrix, ["FITC-A", "PE-A", "PerCP-A"])
```

## Gating Workflows

### Polygon Gates

Multi-vertex regions for complex cell populations:

```python
from amprenta_rag.flow_cytometry.gating import apply_polygon_gate

# Define lymphocyte gate
vertices = [
    [20000, 15000], [80000, 15000], [120000, 60000], 
    [100000, 120000], [40000, 100000], [15000, 40000]
]

mask = apply_polygon_gate(events, x_idx=0, y_idx=1, vertices=vertices)
lymphocytes = events[mask]
```

### Rectangle Gates

Simple rectangular regions:

```python
from amprenta_rag.flow_cytometry.gating import apply_rectangle_gate

# CD3+ T cell gate
bounds = {"x_min": 1000, "x_max": 100000, "y_min": 10000, "y_max": 150000}
mask = apply_rectangle_gate(events, x_idx=2, y_idx=1, bounds=bounds)
```

### Quadrant Gates

Four-region division:

```python
from amprenta_rag.flow_cytometry.gating import apply_quadrant_gate

# CD4/CD8 analysis
quadrants = apply_quadrant_gate(events, x_idx=3, y_idx=4, x_thresh=1000, y_thresh=1000)

cd4_single = events[quadrants["Q1"]]  # CD4+CD8-
cd8_single = events[quadrants["Q4"]]  # CD4-CD8+
double_pos = events[quadrants["Q2"]]  # CD4+CD8+
double_neg = events[quadrants["Q3"]]  # CD4-CD8-
```

### Boolean Gates

Combine multiple gates with logical operations:

```python
from amprenta_rag.flow_cytometry.gating import apply_boolean_gate

# T helper cells = CD3+ AND CD4+ AND CD8-
masks = {
    "cd3_gate": cd3_mask,
    "cd4_gate": cd4_mask, 
    "cd8_gate": cd8_mask
}

# CD3+ AND CD4+ (intermediate)
cd3_cd4_mask = apply_boolean_gate(
    {"cd3": masks["cd3_gate"], "cd4": masks["cd4_gate"]}, 
    operator="AND", 
    operand_ids=["cd3", "cd4"]
)

# Final: CD3+CD4+ AND NOT CD8+
th_mask = apply_boolean_gate(
    {"cd3_cd4": cd3_cd4_mask, "cd8": masks["cd8_gate"]},
    operator="AND_NOT",
    operand_ids=["cd3_cd4", "cd8"]
)
```

## Population Statistics

### Basic Statistics

```python
from amprenta_rag.flow_cytometry.gating import compute_population_stats

stats = compute_population_stats(
    events=events,
    mask=gate_mask,
    param_names=["FSC-A", "SSC-A", "FITC-A", "PE-A"],
    parent_event_count=len(parent_events),
    total_events=len(all_events)
)

print(f"Events: {stats.n_events:,}")
print(f"% of Parent: {stats.pct_of_parent:.1f}%")
print(f"% of Total: {stats.pct_of_total:.1f}%")
print(f"Medians: {dict(zip(param_names, stats.median_values))}")
```

### Available Metrics

| Metric | Description | Use Case |
|--------|-------------|----------|
| `n_events` | Event count in gate | Population size |
| `pct_of_parent` | Percentage of parent population | Gating efficiency |
| `pct_of_total` | Percentage of total events | Population frequency |
| `median_values` | Median fluorescence intensity (MFI) | Expression levels |
| `mean_values` | Mean fluorescence intensity | Expression levels |
| `cv_values` | Coefficient of variation (%) | Population homogeneity |

## API Reference

### Authentication

All API endpoints require authentication via JWT token:

```bash
curl -H "Authorization: Bearer <token>" \
  http://localhost:8000/api/v1/flow-cytometry/datasets
```

### Endpoints

#### 1. Upload FCS File

```http
POST /api/v1/flow-cytometry/upload
Content-Type: multipart/form-data

file: <FCS file>
```

**Response:**
```json
{
  "flow_dataset_id": "uuid",
  "dataset_id": "uuid", 
  "filename": "sample.fcs",
  "file_size_bytes": 1048576,
  "processing_status": "completed",
  "message": "FCS file uploaded successfully"
}
```

#### 2. List Datasets

```http
GET /api/v1/flow-cytometry/datasets?offset=0&limit=50&processing_status=completed
```

**Response:**
```json
{
  "items": [
    {
      "id": "uuid",
      "processing_status": "completed",
      "n_events": 100000,
      "n_parameters": 12,
      "cytometer_model": "BD FACSCanto II",
      "sample_id": "SAMPLE_001"
    }
  ],
  "total": 1,
  "offset": 0,
  "limit": 50
}
```

#### 3. Get Dataset Details

```http
GET /api/v1/flow-cytometry/datasets/{dataset_id}
```

#### 4. Get Parameters

```http
GET /api/v1/flow-cytometry/datasets/{dataset_id}/parameters
```

**Response:**
```json
[
  {
    "id": "uuid",
    "parameter_name": "FSC-A",
    "parameter_index": 0,
    "range_min": 0.0,
    "range_max": 262144.0,
    "display_name": "Forward Scatter Area"
  }
]
```

#### 5. Get Events (Paginated)

```http
GET /api/v1/flow-cytometry/datasets/{dataset_id}/events?offset=0&limit=10000&subsample=true
```

**Response:**
```json
{
  "events": [
    {
      "FSC-A": 45231.2,
      "SSC-A": 23451.8,
      "FITC-A": 1234.5
    }
  ],
  "total_events": 100000,
  "returned_events": 10000
}
```

#### 6. Create Gate

```http
POST /api/v1/flow-cytometry/datasets/{dataset_id}/gates
Content-Type: application/json

{
  "gate_name": "Lymphocytes",
  "gate_type": "polygon",
  "x_parameter": "FSC-A",
  "y_parameter": "SSC-A", 
  "gate_definition": {
    "vertices": [[20000, 15000], [80000, 15000]]
  }
}
```

#### 7. List Gates

```http
GET /api/v1/flow-cytometry/datasets/{dataset_id}/gates
```

#### 8. Update Gate

```http
PUT /api/v1/flow-cytometry/gates/{gate_id}
Content-Type: application/json

{
  "gate_name": "Updated Name",
  "is_active": false
}
```

#### 9. Delete Gate

```http
DELETE /api/v1/flow-cytometry/gates/{gate_id}
```

#### 10. Get Population Statistics

```http
GET /api/v1/flow-cytometry/datasets/{dataset_id}/populations
```

**Response:**
```json
[
  {
    "id": "uuid",
    "population_name": "Lymphocytes Population",
    "n_events": 45000,
    "pct_of_parent": 85.2,
    "pct_of_total": 45.0,
    "median_values": {"FSC-A": 52341, "SSC-A": 23451},
    "mean_values": {"FSC-A": 54123, "SSC-A": 25678},
    "cv_values": {"FSC-A": 15.2, "SSC-A": 18.5}
  }
]
```

## Dashboard User Guide

### Accessing the Dashboard

1. Start the Streamlit dashboard:
   ```bash
   streamlit run scripts/dashboard/app.py
   ```

2. Navigate to: **Analysis → Flow Cytometry**

### Tab 1: Upload 📤

**Upload FCS Files:**
1. Click "Choose an FCS file"
2. Select your .fcs file
3. Click "Upload and Process"
4. Monitor processing status

**Existing Datasets:**
- View all uploaded datasets
- Select datasets for analysis
- Check processing status and metadata

### Tab 2: Scatter Plots 📊

**Create Visualizations:**
1. Select X and Y parameters
2. Choose transformation (linear, log, logicle, arcsinh)
3. Toggle density coloring
4. Click "Generate Plot"

**Features:**
- Interactive Plotly charts
- Zoom and pan capabilities
- Automatic subsampling for performance
- Real-time parameter updates

### Tab 3: Gating 🎯

**Create Gates:**
1. Select gate type (polygon, rectangle, quadrant)
2. Enter gate name
3. Define gate parameters:
   - **Polygon**: JSON vertex coordinates
   - **Rectangle**: Min/max bounds
   - **Quadrant**: X/Y thresholds
4. Click "Create Gate"

**Manage Gates:**
- View existing gates table
- Deactivate or delete gates
- Edit gate properties

### Tab 4: Statistics 📈

**View Population Data:**
- Population statistics table
- Summary metrics
- Distribution pie charts
- Channel statistics (expandable)

**Export Data:**
- Download CSV files
- Population summaries
- Statistical reports

## Troubleshooting

### Common Issues

#### FCS File Upload Errors

**Error: "Invalid file format"**
- **Cause**: File extension is not .fcs
- **Solution**: Rename file with .fcs extension

**Error: "Failed to parse FCS file"**
- **Cause**: Corrupted or non-standard FCS file
- **Solutions**:
  - Verify file integrity
  - Check FCS version compatibility
  - Try re-exporting from original software

**Error: "File size too large"**
- **Cause**: File exceeds 500MB limit
- **Solution**: Subsample data in acquisition software

#### Processing Issues

**Status stuck at "pending"**
- **Cause**: Background processing error
- **Solutions**:
  - Check server logs
  - Restart processing service
  - Verify disk space

**Error: "Failed to load events"**
- **Cause**: Parquet file corruption or missing
- **Solution**: Re-upload and re-process FCS file

#### Gating Problems

**Gate creation fails**
- **Cause**: Invalid gate definition
- **Solutions**:
  - Check parameter names match dataset
  - Verify coordinate ranges
  - Use valid JSON format for polygons

**Population statistics missing**
- **Cause**: Gate processing not completed
- **Solution**: Wait for background processing to complete

### Performance Optimization

#### Large Datasets (>1M events)

1. **Enable subsampling**:
   ```python
   # API request with subsampling
   response = requests.get(
       f"/api/v1/flow-cytometry/datasets/{dataset_id}/events?subsample=true&limit=50000"
   )
   ```

2. **Use pagination**:
   ```python
   # Process in chunks
   offset = 0
   limit = 10000
   while True:
       response = requests.get(
           f"/api/v1/flow-cytometry/datasets/{dataset_id}/events?offset={offset}&limit={limit}"
       )
       events = response.json()["events"]
       if not events:
           break
       # Process chunk
       offset += limit
   ```

3. **Optimize transformations**:
   - Apply transforms during upload, not visualization
   - Cache transformed data
   - Use appropriate data types (float32 vs float64)

#### Memory Management

- Monitor memory usage during large file processing
- Use streaming for event data access
- Implement garbage collection for temporary objects

### Error Codes

| Code | Error | Description | Solution |
|------|-------|-------------|----------|
| 400 | Bad Request | Invalid parameters or file format | Check request format |
| 404 | Not Found | Dataset or gate not found | Verify IDs exist |
| 413 | File Too Large | File exceeds size limit | Reduce file size |
| 422 | Validation Error | Invalid gate definition or parameters | Fix input validation |
| 500 | Internal Error | Server processing error | Check logs, contact admin |

### Getting Help

1. **Check logs**: Application logs contain detailed error information
2. **API documentation**: Swagger UI at `/docs` endpoint
3. **Integration tests**: Run test suite to verify functionality
4. **Seed data**: Use demo data to test workflows

## Development

### Running Tests

```bash
# Unit tests
pytest amprenta_rag/tests/test_flow_cytometry_*.py -v

# Integration tests  
pytest amprenta_rag/tests/integration/test_flow_cytometry_integration.py -v

# API tests
pytest amprenta_rag/tests/api/test_flow_cytometry_api.py -v

# E2E dashboard tests
pytest amprenta_rag/tests/e2e/test_flow_cytometry_page.py -v
```

### Seed Demo Data

```bash
# Create synthetic demo data
python scripts/seed_flow_cytometry_data.py

# Reset and recreate
python scripts/seed_flow_cytometry_data.py --reset
```

### Code Structure

```
amprenta_rag/flow_cytometry/
├── fcs_parser.py          # FCS file parsing
├── transforms.py          # Data transformations  
├── gating.py             # Gating algorithms
└── ingest_service.py     # Orchestration service

amprenta_rag/api/routers/
└── flow_cytometry.py     # API endpoints

scripts/dashboard/pages/
└── flow_cytometry.py     # Dashboard interface

amprenta_rag/database/
└── models_flow_cytometry.py  # Database schema
```

### Contributing

1. Follow existing code patterns and style
2. Add comprehensive tests for new features
3. Update documentation for API changes
4. Run full test suite before submitting
5. Follow semantic commit message format

## Changelog

### Version 1.0.0 (2025-12-30)

**Initial Release - Complete Flow Cytometry Support**

✅ **Core Features:**
- FCS 2.0/3.0/3.1 file parsing
- Logicle and arcsinh transformations
- Polygon, rectangle, quadrant, and boolean gating
- Population statistics calculation
- Hierarchical gating support

✅ **API Endpoints:**
- 9 RESTful endpoints for complete workflow
- File upload with validation
- Paginated event data access
- Full gate CRUD operations
- Population statistics retrieval

✅ **Dashboard Interface:**
- 4-tab Streamlit interface
- Interactive plotting with Plotly
- Real-time gate creation and management
- CSV export functionality

✅ **Infrastructure:**
- PostgreSQL database schema
- Parquet-based event storage
- Background processing with threading
- Comprehensive test coverage (71+ tests)

**Statistics:**
- **6 Batches** completed
- **71+ Tests** passing
- **15 E2E Tests** for dashboard
- **6 Integration Tests** for complex workflows
- **Complete Documentation** with troubleshooting guide
