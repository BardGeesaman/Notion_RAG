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

## What Happens During Import

1. **Repository Harvest**
   - Fetches study metadata from repository API
   - Extracts omics type, disease, organism, sample type
   - Retrieves study files/metadata

2. **Postgres Dataset Creation**
   - Creates `Dataset` record in Postgres
   - Stores metadata (name, description, omics_type, etc.)
   - Stores external repository IDs

3. **Feature Extraction** (if `--ingest`)
   - Extracts features from study files
   - Links features to dataset in Postgres
   - Normalizes feature names

4. **Pinecone Embedding** (if `--ingest`)
   - Generates text content from dataset metadata
   - Chunks and embeds text
   - Upserts vectors to Pinecone

5. **Notion Sync** (if `--create-notion`)
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

### 3. Import in Batches

For large discovery results, import in smaller batches:

```bash
# Import first 10 studies
python scripts/batch_import_repository_studies.py \
    --studies $(jq -r '.results.GEO[:10] | .[] | "GEO:\(.)"' preview.json) \
    --create-postgres \
    --ingest
```

### 4. Monitor Progress

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

