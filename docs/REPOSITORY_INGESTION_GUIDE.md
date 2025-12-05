# Repository Ingestion Guide

**Last Updated**: 2025-12-04  
**Status**: ✅ Ready to Use

## Overview

The platform supports importing data from multiple public omics repositories. This allows you to discover and harvest studies directly into your system.

## Available Repositories

### 1. GEO (Gene Expression Omnibus)
- **Omics Type**: Transcriptomics
- **Study ID Format**: `GSE12345`
- **Source**: NCBI GEO database
- **Use Case**: RNA-seq, microarray gene expression data

### 2. PRIDE (PRIDE Archive)
- **Omics Type**: Proteomics
- **Project ID Format**: `PXD012345`
- **Source**: EBI PRIDE Archive
- **Use Case**: Mass spectrometry proteomics data

### 3. MetaboLights
- **Omics Type**: Metabolomics
- **Study ID Format**: `MTBLS123`
- **Source**: EBI MetaboLights
- **Use Case**: Metabolomics studies

### 4. MW (Metabolomics Workbench)
- **Omics Types**: Metabolomics & Lipidomics
- **Study ID Format**: `ST001111`
- **Source**: Metabolomics Workbench
- **Use Case**: Metabolomics and lipidomics studies

## Workflow

### Step 1: Discover Studies

Search for studies matching your criteria:

```bash
# Search GEO for Alzheimer studies
python scripts/discover_omics_studies.py \
  --repository GEO \
  --keywords "Alzheimer" \
  --max-results 10

# Search Metabolomics Workbench for ALS studies
python scripts/discover_omics_studies.py \
  --repository MW \
  --keywords "ALS" \
  --max-results 10

# Search across all repositories
python scripts/discover_omics_studies.py \
  --keywords "dementia" \
  --omics-type metabolomics \
  --max-results 20
```

**Options**:
- `--repository`: Specific repository (GEO, PRIDE, MetaboLights, MW, MW_LIPIDOMICS, MW_METABOLOMICS)
- `--keywords`: Search terms (space-separated)
- `--omics-type`: Filter by type (transcriptomics, proteomics, metabolomics, lipidomics)
- `--disease`: Filter by disease
- `--organism`: Filter by organism (e.g., "Homo sapiens")
- `--sample-type`: Filter by sample type (e.g., "CSF", "plasma")
- `--max-results`: Maximum number of results (default: 50)
- `--output`: Save results to JSON file

### Step 2: Harvest a Study

Harvest metadata and optionally create dataset pages:

```bash
# Dry run (just show what would be done)
python scripts/harvest_repository_study.py \
  --study-id GSE12345 \
  --repository GEO \
  --dry-run

# Create Notion page
python scripts/harvest_repository_study.py \
  --study-id GSE12345 \
  --repository GEO \
  --create-notion

# Create Notion page and trigger ingestion
python scripts/harvest_repository_study.py \
  --study-id GSE12345 \
  --repository GEO \
  --create-notion \
  --ingest
```

**Options**:
- `--study-id`: Repository-specific study identifier (required)
- `--repository`: Repository name (required)
- `--create-notion`: Create/update Notion Dataset page
- `--ingest`: Trigger dataset ingestion after creating page
- `--dry-run`: Only print what would be done (no changes)

## Examples

### Example 1: Discover and Harvest GEO Study

```bash
# 1. Search for Alzheimer transcriptomics studies
python scripts/discover_omics_studies.py \
  --repository GEO \
  --keywords "Alzheimer transcriptomics" \
  --max-results 5

# Output shows study IDs like:
# GEO: 5 studies
#   1. GSE12345
#   2. GSE67890
#   ...

# 2. Harvest a specific study
python scripts/harvest_repository_study.py \
  --study-id GSE12345 \
  --repository GEO \
  --create-notion \
  --ingest
```

### Example 2: Harvest Metabolomics Workbench Study

```bash
# Harvest MW study (automatically detects omics type)
python scripts/harvest_repository_study.py \
  --study-id ST001111 \
  --repository MW \
  --create-notion \
  --ingest
```

### Example 3: Batch Discovery and Export

```bash
# Discover and save results
python scripts/discover_omics_studies.py \
  --repository MW \
  --keywords "ALS" \
  --max-results 50 \
  --output als_studies.json

# Results saved to als_studies.json for later processing
```

## What Happens During Harvest

1. **Fetch Metadata**: Retrieves study metadata from repository API
2. **Create/Update Notion Page**: Creates Dataset page in Notion (if `--create-notion`)
3. **Add mwTab Data**: For MW studies, adds mwTab data to page
4. **Trigger Ingestion**: If `--ingest`, triggers full dataset ingestion:
   - Downloads data files
   - Extracts features
   - Creates Postgres dataset
   - Links features
   - Embeds in Pinecone

## Integration with Postgres

**Note**: The harvest script currently creates Notion pages. When `--ingest` is used, the ingestion pipeline will:
- Create Postgres dataset (if `USE_POSTGRES_AS_SOT=true`)
- Link features to Postgres
- Optionally sync to Notion

## List Available Repositories

```bash
python scripts/discover_omics_studies.py --list-repositories
```

## Error Handling

The scripts include:
- Rate limiting (respects repository API limits)
- Retry logic for failed requests
- Graceful error handling
- Detailed logging

## Configuration

No additional configuration needed! The scripts use:
- Repository APIs (public, no authentication required)
- Optional: `GEO_API_KEY` for higher NCBI rate limits (set in `.env`)

## Tips

1. **Start Small**: Test with `--dry-run` first
2. **Use Filters**: Narrow down results with `--disease`, `--organism`, etc.
3. **Save Results**: Use `--output` to save discovery results for batch processing
4. **Check Logs**: Watch for rate limit warnings and errors

## Troubleshooting

### "Repository not found"
- Check repository name spelling (case-sensitive)
- Use `--list-repositories` to see available names

### "Study not found"
- Verify study ID format (GSE, PXD, MTBLS, ST)
- Check repository website to confirm study exists

### Rate Limiting
- Scripts automatically respect rate limits
- For GEO, set `GEO_API_KEY` in `.env` for higher limits

### Ingestion Errors
- Check logs for specific errors
- Verify data files are accessible
- Ensure Postgres connection is configured

## Next Steps

After harvesting:
1. Check dashboard for imported datasets
2. Verify feature linking worked
3. Run signature matching
4. Query via RAG

## See Also

- `docs/POSTGRES_MIGRATION_GUIDE.md` - Postgres integration details
- `docs/NEXT_STEPS.md` - Roadmap and priorities

