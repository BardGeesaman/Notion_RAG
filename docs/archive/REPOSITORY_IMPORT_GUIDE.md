# Repository Import Guide - All Omics Types

## Overview

This guide explains how to discover and import all types of omics data from public repositories:
- **GEO** - Transcriptomics (gene expression)
- **PRIDE** - Proteomics
- **MetaboLights** - Metabolomics
- **Metabolomics Workbench (MW)** - Metabolomics and Lipidomics

---

## Quick Start

### 1. Discover Studies

```bash
# Discover all omics types for a disease
python scripts/import_all_omics_repositories.py \
    --disease "ALS" \
    --discover-only \
    --output discovered_studies.json

# Discover specific omics type
python scripts/import_all_omics_repositories.py \
    --omics-type metabolomics \
    --keywords "Alzheimer" "amyloid" \
    --discover-only \
    --output metabolomics_studies.json

# Discover from specific repository
python scripts/import_all_omics_repositories.py \
    --repository MW_LIPIDOMICS \
    --keywords "ceramide" \
    --discover-only \
    --output lipidomics_studies.json
```

### 2. Import Studies

```bash
# Import discovered studies (Postgres-only, no Notion)
python scripts/import_all_omics_repositories.py \
    --input discovered_studies.json \
    --import \
    --ingest

# Import with Notion sync (optional)
python scripts/import_all_omics_repositories.py \
    --input discovered_studies.json \
    --import \
    --ingest \
    --create-notion
```

### 3. Discover and Import in One Step

```bash
# Discover and import transcriptomics studies
python scripts/import_all_omics_repositories.py \
    --omics-type transcriptomics \
    --disease "Parkinson" \
    --import \
    --ingest

# Discover and import all omics types for a disease
python scripts/import_all_omics_repositories.py \
    --disease "Alzheimer" \
    --max-results 20 \
    --import \
    --ingest
```

---

## Available Repositories

List all available repositories:

```bash
python scripts/import_all_omics_repositories.py --list-repositories
```

Output:
```
Available repositories:
  - MW
  - MW_LIPIDOMICS
  - MW_METABOLOMICS
  - GEO
  - PRIDE
  - MetaboLights
```

---

## Discovery Options

### By Disease

```bash
python scripts/import_all_omics_repositories.py \
    --disease "ALS" \
    --discover-only \
    --output als_studies.json
```

### By Keywords

```bash
python scripts/import_all_omics_repositories.py \
    --keywords "Alzheimer" "amyloid" "tau" \
    --discover-only \
    --output alzheimer_studies.json
```

### By Omics Type

```bash
# Transcriptomics (GEO)
python scripts/import_all_omics_repositories.py \
    --omics-type transcriptomics \
    --keywords "ALS" \
    --discover-only

# Proteomics (PRIDE)
python scripts/import_all_omics_repositories.py \
    --omics-type proteomics \
    --keywords "Alzheimer" \
    --discover-only

# Metabolomics (MetaboLights, MW_METABOLOMICS)
python scripts/import_all_omics_repositories.py \
    --omics-type metabolomics \
    --keywords "Parkinson" \
    --discover-only

# Lipidomics (MW_LIPIDOMICS)
python scripts/import_all_omics_repositories.py \
    --omics-type lipidomics \
    --keywords "ceramide" \
    --discover-only
```

### Combined Filters

```bash
python scripts/import_all_omics_repositories.py \
    --disease "ALS" \
    --organism "Homo sapiens" \
    --sample-type "CSF" \
    --omics-type metabolomics \
    --max-results 30 \
    --discover-only \
    --output als_csf_metabolomics.json
```

---

## Import Options

### Basic Import (Postgres-only)

```bash
python scripts/import_all_omics_repositories.py \
    --input discovered_studies.json \
    --import
```

This creates datasets in Postgres only (no Notion).

### Import with Ingestion

```bash
python scripts/import_all_omics_repositories.py \
    --input discovered_studies.json \
    --import \
    --ingest
```

This:
1. Creates datasets in Postgres
2. Extracts features
3. Embeds to Pinecone
4. Enables RAG queries

### Import with Notion Sync (Optional)

```bash
python scripts/import_all_omics_repositories.py \
    --input discovered_studies.json \
    --import \
    --ingest \
    --create-notion
```

This also creates/updates Notion pages for documentation.

### Stop on Error

```bash
python scripts/import_all_omics_repositories.py \
    --input discovered_studies.json \
    --import \
    --stop-on-error
```

Stops batch import on first error (useful for debugging).

---

## Complete Workflow Examples

### Example 1: Import All Omics for ALS

```bash
# Step 1: Discover studies
python scripts/import_all_omics_repositories.py \
    --disease "ALS" \
    --max-results 50 \
    --discover-only \
    --output als_all_omics.json

# Step 2: Review discovered studies
cat als_all_omics.json | jq '.summary'

# Step 3: Import all discovered studies
python scripts/import_all_omics_repositories.py \
    --input als_all_omics.json \
    --import \
    --ingest
```

### Example 2: Import Specific Omics Types

```bash
# Transcriptomics
python scripts/import_all_omics_repositories.py \
    --omics-type transcriptomics \
    --disease "Alzheimer" \
    --import \
    --ingest

# Proteomics
python scripts/import_all_omics_repositories.py \
    --omics-type proteomics \
    --disease "Alzheimer" \
    --import \
    --ingest

# Metabolomics
python scripts/import_all_omics_repositories.py \
    --omics-type metabolomics \
    --disease "Alzheimer" \
    --import \
    --ingest

# Lipidomics
python scripts/import_all_omics_repositories.py \
    --omics-type lipidomics \
    --disease "Alzheimer" \
    --import \
    --ingest
```

### Example 3: Targeted Discovery and Import

```bash
# Discover lipidomics studies about ceramides in CSF
python scripts/import_all_omics_repositories.py \
    --repository MW_LIPIDOMICS \
    --keywords "ceramide" \
    --sample-type "CSF" \
    --import \
    --ingest
```

---

## Output Files

### Discovery Results JSON Format

```json
{
  "results": {
    "GEO": ["GSE12345", "GSE67890"],
    "PRIDE": ["PXD001234", "PXD005678"],
    "MetaboLights": ["MTBLS123", "MTBLS456"],
    "MW_LIPIDOMICS": ["ST004168", "ST004169"],
    "MW_METABOLOMICS": ["ST001234", "ST001567"]
  },
  "summary": {
    "total_repositories": 5,
    "total_studies": 10,
    "by_repository": {
      "GEO": 2,
      "PRIDE": 2,
      "MetaboLights": 2,
      "MW_LIPIDOMICS": 2,
      "MW_METABOLOMICS": 2
    }
  }
}
```

---

## Study Grouping

### What is Study Grouping?

**Study grouping** automatically creates **Experiment** entities to organize related datasets from the same repository study.

**Why It Matters**:
- Multi-dataset studies stay organized
- Relationships preserved from source repository
- Easy to find all datasets from one study
- Maintains experimental context (case vs control, replicates, cohorts)

### How It Works

When you import a repository study, the system:

1. **Creates an Experiment** - One per repository study
2. **Extracts metadata** - From repository study record (title, description, methods, etc.)
3. **Links all datasets** - All datasets from that study automatically linked to the Experiment

### Example: Multi-Dataset Study

**Repository**: Metabolomics Workbench, Study ID `ST004396`

**Study Title**: "Lipidomics analysis of CSF from ALS patients"

**Datasets in Study**:
- ALS patient samples (n=30)
- Healthy control samples (n=30)  
- Independent validation cohort (n=15)

**What Gets Created**:

```
Experiment: ST004396 - ALS CSF Lipidomics
├── Description: Multi-cohort lipidomics study...
├── Disease: [ALS]
├── Sample Type: [CSF]
├── Organism: [Human]
└── Linked Datasets:
    ├── Dataset 1: ST004396_ALS_patients (n=30)
    ├── Dataset 2: ST004396_controls (n=30)
    └── Dataset 3: ST004396_validation (n=15)
```

**Result**: All three datasets grouped under one Experiment, making it easy to:
- Find all datasets from this study
- Understand case/control relationships
- Compare results across cohorts
- Maintain study context for interpretation

### How to Enable Study Grouping

Study grouping is enabled by using the `--create-experiment` flag:

```bash
# Single study import with experiment creation
python scripts/ingest_from_repository.py \
    --repository MW \
    --study-id ST004396 \
    --create-experiment \
    --ingest

# Batch import with experiment creation
python scripts/import_all_omics_repositories.py \
    --disease "ALS" \
    --import \
    --ingest \
    --create-experiment
```

**Note**: If `--create-experiment` is NOT used, datasets are still imported but not grouped under an Experiment entity.

### Viewing Grouped Studies

#### In Dashboard

1. Navigate to **Experiments** page
2. Select experiment by name
3. View linked datasets in "Datasets" section

#### In Postgres

```sql
-- Find all experiments with their linked datasets
SELECT 
    e.name AS experiment_name,
    e.description,
    COUNT(d.id) AS dataset_count
FROM experiments e
LEFT JOIN dataset_experiment de ON e.id = de.experiment_id
LEFT JOIN datasets d ON de.dataset_id = d.id
GROUP BY e.id, e.name, e.description
ORDER BY dataset_count DESC;

-- Get all datasets for a specific experiment
SELECT 
    d.name AS dataset_name,
    d.omics_type,
    d.data_origin,
    d.created_at
FROM datasets d
JOIN dataset_experiment de ON d.id = de.dataset_id
JOIN experiments e ON de.experiment_id = e.id
WHERE e.name LIKE '%ST004396%'
ORDER BY d.created_at;
```

### Study Grouping by Repository

Different repositories structure studies differently:

#### Metabolomics Workbench (MW)

- **Study ID**: `ST001234`
- **Grouping**: All datasets with same Study ID grouped together
- **Typical Structure**: Case datasets + control datasets + metadata

**Example**:
```bash
python scripts/ingest_from_repository.py \
    --repository MW \
    --study-id ST004396 \
    --create-experiment \
    --ingest
```

Result: One Experiment linking all MW datasets from ST004396

#### GEO (Transcriptomics)

- **Study ID**: `GSE12345`
- **Grouping**: All samples/datasets from same GEO Series
- **Typical Structure**: Expression matrices for different conditions/timepoints

**Example**:
```bash
python scripts/ingest_from_repository.py \
    --repository GEO \
    --study-id GSE12345 \
    --create-experiment \
    --ingest
```

Result: One Experiment linking all datasets from GSE12345

#### PRIDE (Proteomics)

- **Study ID**: `PXD012345`
- **Grouping**: All datasets from same PRIDE project
- **Typical Structure**: Protein abundance matrices for different samples

**Example**:
```bash
python scripts/ingest_from_repository.py \
    --repository PRIDE \
    --study-id PXD012345 \
    --create-experiment \
    --ingest
```

Result: One Experiment linking all datasets from PXD012345

#### MetaboLights

- **Study ID**: `MTBLS123`
- **Grouping**: All datasets from same MetaboLights study
- **Typical Structure**: Metabolite abundance files + sample metadata

**Example**:
```bash
python scripts/ingest_from_repository.py \
    --repository MetaboLights \
    --study-id MTBLS123 \
    --create-experiment \
    --ingest
```

Result: One Experiment linking all datasets from MTBLS123

### Multi-Dataset Study Scenarios

#### Scenario 1: Case vs Control

**Study**: Alzheimer's disease metabolomics

**Datasets**:
- Dataset 1: Alzheimer's patient plasma (n=50)
- Dataset 2: Healthy control plasma (n=50)

**Experiment**:
```
Name: MTBLS789 - Alzheimer's Plasma Metabolomics
Linked Datasets: 2
  - MTBLS789_AD_patients (n=50)
  - MTBLS789_controls (n=50)
```

**Query Example**:
```python
# Find all datasets from this study
experiment = db.query(Experiment).filter(
    Experiment.name.like('%MTBLS789%')
).first()

datasets = experiment.datasets
print(f"Found {len(datasets)} datasets from study")
```

#### Scenario 2: Multi-Cohort Study

**Study**: ALS CSF lipidomics with validation cohort

**Datasets**:
- Dataset 1: Discovery cohort - ALS patients (n=30)
- Dataset 2: Discovery cohort - Controls (n=30)
- Dataset 3: Validation cohort - ALS patients (n=15)
- Dataset 4: Validation cohort - Controls (n=15)

**Experiment**:
```
Name: ST004396 - ALS CSF Lipidomics
Linked Datasets: 4
  - ST004396_discovery_ALS (n=30)
  - ST004396_discovery_controls (n=30)
  - ST004396_validation_ALS (n=15)
  - ST004396_validation_controls (n=15)
```

**Analysis Example**: Compare discovery vs validation results by querying all 4 datasets together

#### Scenario 3: Time Course Study

**Study**: Drug response transcriptomics

**Datasets**:
- Dataset 1: Baseline (t=0h)
- Dataset 2: Early response (t=4h)
- Dataset 3: Mid response (t=8h)
- Dataset 4: Late response (t=24h)

**Experiment**:
```
Name: GSE98765 - Drug Response Time Course
Linked Datasets: 4
  - GSE98765_t0h_baseline
  - GSE98765_t4h_early
  - GSE98765_t8h_mid
  - GSE98765_t24h_late
```

**Analysis Example**: Track gene expression changes over time by analyzing all 4 timepoint datasets

### Editing Experiment Metadata

After import, you can edit Experiment metadata via the dashboard:

1. Navigate to **Data Management** → **Edit Metadata** tab
2. Select entity type: **Experiment**
3. Choose experiment from dropdown
4. Edit fields:
   - Name
   - Description
   - Disease
   - Sample Type
   - Organism
   - Model System
   - Methods
   - Summary
   - Results
   - Date

See [Metadata Editing Guide](METADATA_EDITING.md) for detailed instructions.

### Linking Datasets to Experiments Manually

If you imported datasets without `--create-experiment`, you can link them manually:

#### Via Dashboard

1. Navigate to **Data Management** → **Link Entities** tab
2. Select link type: **Dataset → Experiment**
3. Choose dataset from dropdown
4. Select one or more experiments (multiselect)
5. Click "🔗 Link Dataset to Experiments"

#### Via Python API

```python
from amprenta_rag.ingestion.postgres_program_experiment_linking import (
    link_dataset_to_programs_and_experiments_in_postgres
)
from amprenta_rag.database.base import get_db

db = next(get_db())

# Link dataset to experiment
result = link_dataset_to_programs_and_experiments_in_postgres(
    dataset_id="your-dataset-uuid",
    experiment_ids=["experiment-uuid-1", "experiment-uuid-2"],
    db=db
)

db.commit()
print(f"Linked to {result['experiments_linked']} experiment(s)")
```

---

## What Happens During Import

1. **Repository Harvest**
   - Fetches study metadata from repository API
   - Extracts omics type, disease, organism, sample type
   - Retrieves study files/metadata

2. **Postgres Dataset Creation**
   - Creates `Dataset` record in Postgres
   - Stores metadata (name, description, omics_type, etc.)
   - Stores external repository IDs

3. **Experiment Creation** (if `--create-experiment`)
   - Creates `Experiment` record with study metadata
   - Links all datasets from study to experiment
   - Preserves study context and relationships

4. **Feature Extraction** (if `--ingest`)
   - Extracts features from study files
   - Links features to dataset in Postgres
   - Normalizes feature names

5. **Pinecone Embedding** (if `--ingest`)
   - Generates text content from dataset metadata
   - Chunks and embeds text
   - Upserts vectors to Pinecone

6. **Notion Sync** (if `--create-notion`)
   - Creates/updates Notion Dataset pages
   - Syncs metadata to Notion

---

## Tips and Best Practices

### 1. Start with Discovery

Always run discovery first to see what's available:

```bash
python scripts/import_all_omics_repositories.py \
    --disease "ALS" \
    --discover-only \
    --output preview.json
```

### 2. Filter Before Importing

Use filters to narrow down results:

```bash
# Only human CSF studies
python scripts/import_all_omics_repositories.py \
    --disease "ALS" \
    --organism "Homo sapiens" \
    --sample-type "CSF" \
    --import \
    --ingest
```

### 3. Use Study Grouping for Multi-Dataset Studies

Always use `--create-experiment` to preserve study structure:

```bash
# Repository studies with multiple datasets
python scripts/import_all_omics_repositories.py \
    --disease "ALS" \
    --import \
    --ingest \
    --create-experiment
```

**Benefits**:
- Maintains case/control relationships
- Preserves cohort structure
- Easier to find related datasets
- Better experimental context

### 4. Import in Batches

For large discovery results, import in smaller batches:

```bash
# Import first 10 studies
python scripts/batch_import_repository_studies.py \
    --studies $(jq -r '.results.GEO[:10] | .[] | "GEO:\(.)"' preview.json) \
    --create-postgres \
    --ingest
```

### 5. Monitor Progress

The script provides detailed logging. Watch for:
- `✓ Successfully imported` - Study imported successfully
- `✗ Failed to import` - Study import failed (check logs)

---

## Troubleshooting

### No Studies Found

- Try broader keywords
- Remove filters
- Check repository availability

### Import Failures

- Check logs for specific errors
- Verify repository study IDs are valid
- Ensure Postgres is running and accessible

### Slow Performance

- Reduce `--max-results` for discovery
- Import studies in smaller batches
- Use `--stop-on-error` to catch issues early

---

## Related Scripts

- `scripts/discover_omics_studies.py` - Standalone discovery script
- `scripts/harvest_repository_study.py` - Import single study
- `scripts/batch_import_repository_studies.py` - Batch import from file

---

## Next Steps

After importing studies:

1. **Browse in Dashboard**
   ```bash
   python scripts/run_dashboard.py
   ```

2. **Query with RAG**
   ```bash
   python scripts/rag_query.py "What datasets are available for ALS?"
   ```

3. **View in Postgres**
   ```sql
   SELECT name, omics_type, disease, organism 
   FROM datasets 
   WHERE omics_type = 'metabolomics';
   ```

---

## Support

For issues or questions:
- Check logs for detailed error messages
- Review repository-specific documentation
- Verify configuration in `.env` file

