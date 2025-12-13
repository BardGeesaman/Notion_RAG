# Postgres-only Ingestion

**Status**: Active Default Ingestion Method
**Last Updated**: December 9, 2025

As of Phase 3, the Amprenta RAG system has fully transitioned to a Postgres-first architecture. Direct Notion harvesting is deprecated.

## Active Ingestion Scripts

Use these scripts for all data ingestion going forward:

### 1. Public Repository Ingestion
**Script**: `scripts/import_all_omics_repositories.py`

Use this for discovering and importing data from public repositories (Metabolomics Workbench, GEO, PRIDE, MetaboLights).

```bash
# Example: Discover and import ALS studies
python scripts/import_all_omics_repositories.py --disease "ALS" --import --ingest
```

### 2. Local File Ingestion
**Script**: `scripts/batch_ingest_omics.py`

Use this for ingesting local CSV/TSV files (Lipidomics, Metabolomics, Proteomics, Transcriptomics).

```bash
# Example: Batch ingest a directory of files
python scripts/batch_ingest_omics.py --directory ./data/incoming --create-pages
```

## Legacy Harvesters

Legacy Notion-only harvesters (e.g., `harvest_lipidomics.py` older versions) have been archived.

**Archive Location**: `md_archive/scripts/legacy/`

These scripts are preserved for historical reference only and should not be used for new data ingestion.

