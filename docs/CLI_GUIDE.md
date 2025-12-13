# Comprehensive User Guide

Complete guide to using all features of the Amprenta RAG System.

## Table of Contents

1. [Introduction](#introduction)
2. [Data Ingestion](#data-ingestion)
3. [Signature Management](#signature-management)
4. [RAG Queries](#rag-queries)
5. [Cross-Omics Reasoning](#cross-omics-reasoning)
6. [Batch Operations](#batch-operations)
7. [Analysis & Reports](#analysis--reports)
8. [Best Practices](#best-practices)

## Introduction

The Amprenta RAG System is a comprehensive multi-omics knowledge management platform that:
- Ingests data from multiple omics types (lipidomics, metabolomics, proteomics, transcriptomics)
- Manages multi-omics signatures
- Provides semantic search across all data sources
- Generates evidence-based reports
- Discovers patterns automatically

## Data Ingestion

### Multi-Omics Ingestion

The system supports ingestion of four omics types:

#### Lipidomics Ingestion

```bash
# Basic ingestion
python scripts/ingest_lipidomics.py \
  --file data/lipidomics.csv \
  --create-page

# Link to programs and experiments
python scripts/ingest_lipidomics.py \
  --file data/lipidomics.csv \
  --create-page \
  --program-id PROGRAM_PAGE_ID \
  --experiment-id EXPERIMENT_PAGE_ID

# Attach to existing dataset page
python scripts/ingest_lipidomics.py \
  --file data/lipidomics.csv \
  --dataset-page-id DATASET_PAGE_ID
```

**File Format**: CSV/TSV with lipid species in various formats:
- `Cer(d18:1/16:0)` (canonical)
- `CER 16:0`
- `Ceramide(d18:1/16:0)`

The system automatically normalizes to canonical format.

#### Metabolomics Ingestion

```bash
python scripts/ingest_metabolomics.py \
  --file data/metabolomics.csv \
  --create-page \
  --program-id PROGRAM_PAGE_ID
```

**File Format**: CSV/TSV with metabolite names (automatically normalized).

#### Proteomics Ingestion

```bash
python scripts/ingest_proteomics.py \
  --file data/proteomics.csv \
  --create-page
```

**File Format**: CSV/TSV with protein identifiers (UniProt, gene symbols, etc.).

#### Transcriptomics Ingestion

```bash
python scripts/ingest_transcriptomics.py \
  --file data/rna_seq.csv \
  --create-page
```

**File Format**: CSV/TSV with gene identifiers (gene symbols, Ensembl IDs, etc.).

### Public Repository Ingestion

#### Discover Studies

```bash
# Search Metabolomics Workbench
python scripts/discover_omics_studies.py \
  --repository MW \
  --keyword "ALS" \
  --keyword "ceramide" \
  --max-results 20

# Search GEO (transcriptomics)
python scripts/discover_omics_studies.py \
  --repository GEO \
  --keyword "neurodegeneration" \
  --max-results 10
```

#### Harvest Studies

```bash
# Harvest from Metabolomics Workbench
python scripts/harvest_repository_study.py \
  --repository MW \
  --study-id ST004396 \
  --create-notion \
  --ingest

# Harvest from GEO
python scripts/harvest_repository_study.py \
  --repository GEO \
  --study-id GSE12345 \
  --create-notion \
  --ingest
```

### Batch Ingestion

#### Basic Batch Processing

```bash
# Auto-detect omics type and ingest all files
python scripts/batch_ingest_omics.py \
  --directory /path/to/omics/data

# Force specific omics type
python scripts/batch_ingest_omics.py \
  --directory /path/to/data \
  --omics-type lipidomics
```

#### Parallel Processing

```bash
# Process files in parallel (faster)
python scripts/batch_ingest_omics.py \
  --directory /path/to/data \
  --parallel \
  --max-workers 5

# Export results to JSON
python scripts/batch_ingest_omics.py \
  --directory /path/to/data \
  --export-results batch_results.json
```

#### Filtering Files

```bash
# Filter by pattern
python scripts/batch_ingest_omics.py \
  --directory /path/to/data \
  --pattern "*lipid*.csv"
```

## Signature Management

### Creating Signatures

#### Manual Signature Creation

Create a TSV file with signature components:

```tsv
feature_type	feature_name	direction	weight
lipid	Cer(d18:1/16:0)	↑	1.0
lipid	Cer(d18:1/18:0)	↑	1.0
gene	TP53	↓	0.8
protein	P04637	↓	0.8
metabolite	Glutamate	↑	1.0
```

Then ingest:

```bash
python scripts/ingest_signature.py \
  --file signature.tsv \
  --signature-type "Literature-derived" \
  --description "Multi-omics signature for ALS"
```

#### Bulk Signature Ingestion

```bash
# Ingest multiple signatures from a directory
python scripts/bulk_ingest_signatures.py \
  --directory /path/to/signatures \
  --pattern "*.tsv"
```

### Automated Signature Discovery

#### Discover Patterns from Datasets

```bash
# Discover signatures from all datasets
python scripts/discover_signatures.py \
  --all-datasets \
  --min-support 3 \
  --min-confidence 0.7

# Discover from specific datasets
python scripts/discover_signatures.py \
  --dataset-ids DATASET_ID_1 DATASET_ID_2 \
  --min-support 2
```

#### Tune Discovery Parameters

```bash
python scripts/discover_signatures.py \
  --all-datasets \
  --min-support 5 \          # Minimum datasets with pattern
  --min-features 5 \         # Minimum features in signature
  --max-features 20 \        # Maximum features in signature
  --min-co-occurrence 3 \    # Minimum co-occurrence count
  --min-confidence 0.8       # Minimum confidence score
```

#### Export and Ingest Discovered Signatures

```bash
# Export to JSON for review
python scripts/discover_signatures.py \
  --all-datasets \
  --output discovered_signatures.json

# Auto-ingest to Notion
python scripts/discover_signatures.py \
  --all-datasets \
  --ingest
```

### Signature Scoring

#### Score a Dataset Against Signatures

```bash
# Score a specific dataset
python scripts/score_signature.py \
  --dataset-id DATASET_PAGE_ID

# Score against specific signature
python scripts/score_signature.py \
  --dataset-id DATASET_PAGE_ID \
  --signature-id SIGNATURE_PAGE_ID
```

## RAG Queries

### Basic Queries

```bash
# Simple semantic search
python scripts/rag_query.py \
  --query "What is the role of ceramides in ALS?"

# With filters
python scripts/rag_query.py \
  --query "ceramide signature" \
  --source-type Dataset \
  --disease ALS \
  --top-k 10
```

### Advanced Filtering

```bash
# Multiple source types
python scripts/rag_query.py \
  --query "neurodegeneration" \
  --source-type Literature Dataset Experiment \
  --top-k 20

# Filter by signature
python scripts/rag_query.py \
  --query "signature patterns" \
  --signature "ALS-CSF-Core-6Ceramides"

# Filter by multiple criteria
python scripts/rag_query.py \
  --query "lipid changes" \
  --disease ALS \
  --source-type Dataset \
  --signature "ALS-CSF-Core-6Ceramides"
```

### Signature Similarity Queries

```bash
# Find similar signatures to a dataset
python scripts/rag_query.py \
  --signature-score DATASET_PAGE_ID

# Find datasets matching a signature
python scripts/rag_query.py \
  --signature-score-from-signature SIGNATURE_PAGE_ID
```

### Cross-Omics Reasoning

```bash
# Program-level summary
python scripts/rag_query.py \
  --cross-omics-program PROGRAM_PAGE_ID

# Signature-level summary
python scripts/rag_query.py \
  --cross-omics-signature SIGNATURE_PAGE_ID

# Feature-level summary
python scripts/rag_query.py \
  --cross-omics-feature "Cer(d18:1/16:0)" lipid

# Dataset-level summary
python scripts/rag_query.py \
  --cross-omics-dataset DATASET_PAGE_ID
```

## Cross-Omics Reasoning

### Program Summaries

```bash
# Generate comprehensive program summary
python scripts/rag_query.py \
  --cross-omics-program PROGRAM_PAGE_ID \
  --top-k-per-omics 20
```

This analyzes:
- All datasets linked to the program
- Multi-omics patterns
- Cross-omics convergence
- Disease/model system context

### Signature Analysis

```bash
# Analyze a signature across all omics
python scripts/rag_query.py \
  --cross-omics-signature SIGNATURE_PAGE_ID \
  --top-k-datasets 20 \
  --top-k-chunks 100
```

### Feature Analysis

```bash
# Analyze a specific feature (e.g., a lipid species)
python scripts/rag_query.py \
  --cross-omics-feature "Cer(d18:1/16:0)" lipid \
  --top-k-datasets 20
```

## Batch Operations

### Batch Feature Extraction

```bash
# Extract features from multiple datasets
python scripts/batch_feature_extraction.py \
  --dataset-ids ID1 ID2 ID3 \
  --output features.json
```

### Batch Signature Scoring

```bash
# Score all datasets against all signatures
python scripts/batch_signature_scoring.py \
  --all-datasets \
  --all-signatures \
  --output scores.json
```

## Analysis & Reports

### Evidence Reports

```bash
# Generate program evidence report
python scripts/generate_evidence_report.py \
  --program-id PROGRAM_PAGE_ID \
  --output program_evidence.md

# Generate dataset report
python scripts/generate_evidence_report.py \
  --dataset-id DATASET_PAGE_ID \
  --output dataset_evidence.md
```

### Dataset Comparison

```bash
# Compare two datasets
python scripts/compare_datasets.py \
  --dataset-id-1 DATASET_ID_1 \
  --dataset-id-2 DATASET_ID_2 \
  --output comparison.json
```

### Program Signature Maps

```bash
# Generate signature map for a program
python scripts/generate_program_signature_map.py \
  --program-id PROGRAM_PAGE_ID \
  --output signature_map.json
```

### Pathway Enrichment

```bash
# Enrichment analysis for a dataset
python scripts/pathway_enrichment.py \
  --dataset-id DATASET_PAGE_ID \
  --output enrichment_results.json
```

## Best Practices

### Data Organization

1. **Use Consistent Naming**: Follow naming conventions for files and datasets
2. **Link to Programs/Experiments**: Always link datasets to programs and experiments
3. **Add Metadata**: Fill in all relevant metadata fields in Notion

### Performance Optimization

1. **Use Feature Caching**: The system automatically caches feature extractions
2. **Parallel Processing**: Use `--parallel` for batch operations
3. **Batch Operations**: Use batch scripts for multiple files

### Signature Management

1. **Review Discovered Signatures**: Always review auto-discovered signatures before ingesting
2. **Document Signatures**: Add descriptions and source information
3. **Version Control**: Use version fields for signature updates

### Query Optimization

1. **Use Filters**: Narrow results with appropriate filters
2. **Limit Results**: Use `--top-k` to limit result size
3. **Use Cross-Omics**: Leverage cross-omics reasoning for comprehensive insights

## Troubleshooting

See the [Troubleshooting Guide](TROUBLESHOOTING.md) for:
- Common errors and solutions
- Performance issues
- Configuration problems
- API errors

## Next Steps

- [Architecture Overview](ARCHITECTURE.md) - Understand the system design
- [API Reference](API_REFERENCE.md) - Detailed API documentation
- [Configuration Guide](CONFIGURATION.md) - Advanced configuration options

