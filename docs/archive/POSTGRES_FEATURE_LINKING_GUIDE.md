# Postgres Feature Linking Guide

**Last Updated**: 2025-12-04  
**Status**: ✅ Complete - All omics pipelines integrated

## Overview

Feature linking in Postgres creates Feature records and links them to datasets via the association table. This provides fast, scalable feature tracking without relying on Notion API calls.

## Architecture

### Postgres Schema

```
features
├── id (UUID, primary key)
├── name (String, indexed)
├── feature_type (String, indexed) - gene, protein, metabolite, lipid
├── normalized_name (String, indexed)
├── aliases (Array[String])
├── external_ids (JSON)
└── notion_page_id (String, optional) - for migration support

dataset_feature (association table)
├── dataset_id (UUID, foreign key)
└── feature_id (UUID, foreign key)
```

### Feature Types

- **GENE**: Transcriptomics features
- **PROTEIN**: Proteomics features
- **METABOLITE**: Metabolomics features
- **LIPID**: Lipidomics features

## Usage

### Automatic Linking (Recommended)

Features are automatically linked during dataset ingestion:

```python
# Metabolomics ingestion automatically links features
from amprenta_rag.ingestion.metabolomics.ingestion import ingest_metabolomics_file

dataset_id = ingest_metabolomics_file(
    file_path="data/metabolomics_sample.csv",
    create_page=False,  # Postgres-only mode
)
# Features are automatically created and linked!
```

### Manual Linking

You can also link features manually:

```python
from amprenta_rag.ingestion.features.postgres_linking import (
    find_or_create_feature_in_postgres,
    link_feature_to_dataset_in_postgres,
    batch_link_features_to_dataset_in_postgres,
)
from amprenta_rag.models.domain import FeatureType
from uuid import UUID

# Single feature
feature = find_or_create_feature_in_postgres(
    name="Glucose",
    feature_type=FeatureType.METABOLITE,
)

# Link to dataset
link_feature_to_dataset_in_postgres(
    feature_id=feature.id,
    dataset_id=UUID("your-dataset-id"),
)

# Batch linking
features = [
    ("Glucose", FeatureType.METABOLITE),
    ("Lactate", FeatureType.METABOLITE),
]
results = batch_link_features_to_dataset_in_postgres(
    features=features,
    dataset_id=dataset_id,
    max_workers=10,
)
```

### Querying Features

```python
from amprenta_rag.ingestion.features.postgres_linking import (
    get_dataset_features_from_postgres,
)
from amprenta_rag.database.base import get_db
from amprenta_rag.database.models import Dataset, Feature

# Get features for a dataset
db = next(get_db())
dataset = db.query(Dataset).filter(Dataset.id == dataset_id).first()

# Via relationship (easiest)
features = dataset.features  # List of Feature objects

# Via query function
from amprenta_rag.models.domain import FeatureType
features = get_dataset_features_from_postgres(
    dataset_id=dataset_id,
    feature_type=FeatureType.METABOLITE,  # Optional filter
)
```

## Configuration

Feature linking is controlled by configuration:

```bash
# Enable feature linking (default: true)
ENABLE_FEATURE_LINKING=true

# Maximum parallel workers for batch linking (default: 10)
FEATURE_LINKING_MAX_WORKERS=10

# Use Postgres as source of truth (required for Postgres feature linking)
USE_POSTGRES_AS_SOT=true
```

## Integration

### Ingestion Pipelines

All omics pipelines now support Postgres feature linking:

1. **Metabolomics**: Links metabolites automatically
2. **Proteomics**: Links proteins automatically
3. **Transcriptomics**: Links genes automatically
4. **Lipidomics**: Links lipid species automatically

### Dashboard

The Streamlit dashboard shows:
- Feature count per dataset
- Feature type breakdown
- Sample feature names

### API

Features are accessible via FastAPI:

```bash
# Get dataset features
GET /api/v1/datasets/{dataset_id}/features

# Get all features
GET /api/v1/features?feature_type=metabolite

# Get feature details
GET /api/v1/features/{feature_id}
```

## Normalization

Features are automatically normalized based on type:

- **Genes**: Uppercase, remove species suffixes
- **Proteins**: Uppercase, remove isoform suffixes
- **Metabolites**: Lowercase, strip adducts
- **Lipids**: Canonical species format

## Performance

### Batch Linking

- **10-100x faster** than Notion API calls
- Parallel processing (configurable workers)
- Automatic deduplication
- Transaction-safe

### Querying

- Direct SQL queries (no API calls)
- Indexed lookups (feature_type, normalized_name)
- Efficient joins via association table

## Migration from Notion

Features can exist in both Postgres and Notion:

- **Postgres**: Primary database (fast, scalable)
- **Notion**: Documentation layer (optional)

If `ENABLE_NOTION_SYNC=true`, features are linked to both:
1. Postgres (primary)
2. Notion (documentation)

## Troubleshooting

### Features Not Linking

**Symptom**: Features appear in file but not linked to dataset

**Solutions**:
1. Check `ENABLE_FEATURE_LINKING=true`
2. Check `USE_POSTGRES_AS_SOT=true`
3. Check logs for errors
4. Verify Postgres connection

### Duplicate Features

**Symptom**: Same feature appears multiple times

**Solution**:
- Features are deduplicated by `normalized_name` + `feature_type`
- Check normalization functions are consistent

### Performance Issues

**Symptom**: Slow batch linking

**Solutions**:
1. Increase `FEATURE_LINKING_MAX_WORKERS`
2. Use batch linking instead of individual calls
3. Check Postgres connection pool size

## Examples

### Complete Ingestion with Feature Linking

```python
from amprenta_rag.ingestion.metabolomics.ingestion import ingest_metabolomics_file

# Ingest with automatic feature linking
dataset_id = ingest_metabolomics_file(
    file_path="data/metabolomics.csv",
    create_page=False,  # Postgres-only
)

# Features are automatically:
# 1. Normalized
# 2. Created in Postgres (if new)
# 3. Linked to dataset
```

### Query Features for Dataset

```python
from amprenta_rag.database.base import get_db
from amprenta_rag.database.models import Dataset

db = next(get_db())
dataset = db.query(Dataset).filter(Dataset.name == "My Dataset").first()

# Get all features
features = dataset.features

# Filter by type
metabolites = [f for f in features if f.feature_type == "metabolite"]

print(f"Dataset has {len(features)} features:")
for feature in features[:10]:
    print(f"  - {feature.name} ({feature.feature_type})")
```

## Next Steps

1. **Feature Browser in Dashboard**: Browse all features with dataset counts
2. **Feature Statistics**: Count features per type, per dataset
3. **Feature Search**: Search features across all datasets
4. **Feature Relationships**: Track which datasets share features

## Resources

- **Module**: `amprenta_rag/ingestion/features/postgres_linking.py`
- **Models**: `amprenta_rag/database/models.py` (Feature, dataset_feature_assoc)
- **Domain Models**: `amprenta_rag/models/domain.py` (FeatureType enum)

