# Streamlit Dashboard Guide

**Status**: ✅ Built and Ready  
**Last Updated**: 2025-01-XX

## Overview

The Streamlit dashboard provides a visual, browser-based interface for exploring and browsing data stored in Postgres. It's a quick, user-friendly alternative to the FastAPI Swagger UI.

## Features

- 📊 **Overview Dashboard**: Statistics and recent datasets
- 📁 **Datasets Browser**: Browse, filter, and search datasets
- 🏥 **Programs Browser**: View and search programs
- 🔬 **Experiments Browser**: Explore experiments
- 🧪 **Features Browser**: Browse biological features
- 📋 **Signatures Browser**: View multi-omics signatures

## Installation

### Install Dependencies

```bash
pip install streamlit pandas
# Or install all requirements
pip install -r requirements.txt
```

## Usage

### Start the Dashboard

```bash
# Option 1: Use the script
python scripts/run_dashboard.py

# Option 2: Use streamlit directly
streamlit run scripts/run_dashboard.py
```

The dashboard will open in your browser at `http://localhost:8501`

### Navigation

Use the sidebar to navigate between pages:
- **Overview**: Statistics and summary
- **Datasets**: Browse datasets with filters
- **Programs**: View programs
- **Experiments**: Explore experiments
- **Features**: Browse features by type
- **Signatures**: View signatures

## Features by Page

### Overview Page

- Total counts for all entity types
- Bar chart of datasets by omics type
- Table of recent datasets

### Datasets Page

- Filter by omics type
- Search by name
- View dataset details:
  - ID, name, omics type
  - Description
  - File paths
  - Disease associations
  - Creation/update timestamps

### Programs Page

- Search by name
- View program details:
  - ID, name, description
  - Disease associations
  - Related datasets count

### Experiments Page

- Search by name
- View experiment details:
  - ID, name, type
  - Description
  - Disease, matrix, model systems
  - Related datasets count

### Features Page

- Filter by feature type (gene, protein, metabolite, lipid)
- Search by name
- View feature details:
  - Name, type, normalized name

### Signatures Page

- Search by name
- View signature details:
  - ID, name, description
  - Modalities
  - Component count

## Configuration

The dashboard connects directly to Postgres using the configuration from `amprenta_rag/config.py`. Make sure your `.env` file has:

```bash
POSTGRES_HOST=localhost
POSTGRES_PORT=5432
POSTGRES_DB=amprenta_rag
POSTGRES_USER=your_user
POSTGRES_PASSWORD=your_password
```

## Customization

### Add New Pages

Edit `scripts/run_dashboard.py` and add new pages in the sidebar:

```python
page = st.sidebar.radio(
    "Select Page",
    ["Overview", "Datasets", "Your New Page"],
)

if page == "Your New Page":
    st.header("Your New Page")
    # Your code here
```

### Add Visualizations

Use Streamlit's charting capabilities:

```python
import pandas as pd
import streamlit as st

# Bar chart
st.bar_chart(df)

# Line chart
st.line_chart(df)

# Area chart
st.area_chart(df)
```

## Performance

- **Caching**: Database sessions are cached for performance
- **Pagination**: Features page limits to 100 results (use filters)
- **Lazy Loading**: Data is loaded on-demand per page

## Troubleshooting

### Dashboard Won't Start

1. Check Postgres connection:
   ```bash
   python scripts/validate_postgres_setup.py
   ```

2. Check dependencies:
   ```bash
   pip install streamlit pandas
   ```

3. Check port availability:
   - Default port is 8501
   - Change with: `streamlit run scripts/run_dashboard.py --server.port 8502`

### No Data Showing

1. Verify Postgres has data:
   ```bash
   psql amprenta_rag -c "SELECT COUNT(*) FROM datasets;"
   ```

2. Check database connection in dashboard logs

### Slow Performance

- Use filters to narrow down results
- The dashboard queries Postgres directly (should be fast)
- If slow, check Postgres indexes

## Next Steps

### Enhancements

1. **Add Charts**: More visualizations (pie charts, histograms)
2. **Add Filters**: More advanced filtering options
3. **Add Export**: Export data to CSV/Excel
4. **Add Editing**: Edit records directly from dashboard
5. **Add Search**: Global search across all entities

### Integration

- Connect to FastAPI for updates
- Add authentication
- Add user preferences
- Add data export functionality

## Comparison: Streamlit vs FastAPI Swagger UI

| Feature | Streamlit | FastAPI Swagger |
|---------|-----------|-----------------|
| **Visual Interface** | ✅ Rich UI | ⚠️ JSON only |
| **Ease of Use** | ✅ Very easy | ⚠️ Technical |
| **Data Exploration** | ✅ Excellent | ⚠️ Limited |
| **API Testing** | ❌ No | ✅ Yes |
| **Programmatic Access** | ❌ No | ✅ Yes |
| **Customization** | ✅ Easy | ⚠️ Limited |

**Recommendation**: Use Streamlit for browsing/exploring, FastAPI for programmatic access.

