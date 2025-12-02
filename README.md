# Amprenta RAG System

A modular RAG (Retrieval-Augmented Generation) system for scientific literature and email/note ingestion, designed specifically for ceramide/sphingolipid neurodegeneration research.

## Overview

This system ingests content from **Zotero** and **Notion** (emails/notes), embeds it using OpenAI, stores vectors in **Pinecone**, and provides semantic search with answer synthesis. It's designed to help researchers quickly find and understand relevant scientific literature and internal notes.

## Architecture

- **Data Sources**: Zotero (literature PDFs/notes), Notion (emails/notes/experiments/datasets), Metabolomics Workbench (lipidomics studies)
- **Vector Store**: Pinecone
- **Embeddings**: OpenAI `text-embedding-3-large`
- **Query/Answer**: OpenAI GPT models
- **Metadata**: Rich semantic metadata (diseases, targets, lipid species, signatures)

## Project Structure

```
amprenta_rag/
├── clients/          # API client wrappers (OpenAI, Pinecone, Notion)
├── config.py         # Configuration (uses environment variables)
├── ingestion/        # Data ingestion pipelines
│   ├── zotero_ingest.py      # Zotero → Notion → Pinecone
│   ├── email_ingestion.py    # Notion emails → Pinecone
│   ├── notion_pages.py       # Notion API helpers
│   ├── metadata_semantic.py  # Semantic metadata extraction
│   └── ...
├── query/            # RAG query engine
│   ├── pinecone_query.py    # Low-level Pinecone operations
│   ├── rag_engine.py        # High-level RAG orchestration
│   └── rag_query_engine.py  # Compatibility wrapper
├── metadata/         # Metadata classification
├── maintenance/      # Cleanup and sync utilities
└── tests/           # Unit tests

scripts/              # Command-line scripts
├── ingest_collection.py         # Ingest Zotero collection
├── ingest_email.py              # Ingest Notion emails/notes
├── ingest_experiment.py         # Ingest Notion experiments
├── ingest_dataset.py            # Ingest Notion datasets
├── harvest_mw_studies.py        # Harvest MW studies → Notion
├── convert_mwtab_to_csv.py      # Convert mwTab to CSV
├── scan_ceramides_in_mwtab_csv.py  # Scan CSV for ceramides
├── rag_query.py                 # Query the RAG system
└── ...
```

## Setup

### 1. Prerequisites

- Python 3.10+
- API keys for:
  - OpenAI
  - Pinecone
  - Notion
  - Zotero

### 2. Installation

```bash
# Clone the repository
git clone https://github.com/BardGeesaman/Notion_RAG.git
cd Notion_RAG

# Create virtual environment (recommended)
python -m venv venv
source venv/bin/activate  # On Windows: venv\Scripts\activate

# Install dependencies
pip install openai pinecone requests pypdf
```

### 3. Configuration

#### Option 1: Use .env file (Recommended)

Copy the example file and fill in your API keys:

```bash
cp .env.example .env
# Edit .env with your actual API keys
```

The `.env` file is automatically loaded by `amprenta_rag.config` if `python-dotenv` is installed. Install it with:

```bash
pip install python-dotenv
```

#### Option 2: Environment Variables

Alternatively, set environment variables:

```bash
export OPENAI_API_KEY="your_openai_key"
export PINECONE_API_KEY="your_pinecone_key"
export NOTION_API_KEY="your_notion_key"
export ZOTERO_API_KEY="your_zotero_key"
```

**Note**: The `.env` file is already in `.gitignore` and will not be committed. The `config.py` file also contains hardcoded database IDs and other constants that may need to be updated for your workspace.

## Usage

### Ingest Zotero Collection

```bash
python scripts/ingest_collection.py \
  --collection-key YOUR_COLLECTION_KEY \
  --parent-type Literature
```

This will:
1. Fetch items from the Zotero collection
2. Download PDFs and extract text
3. Chunk and embed using OpenAI
4. Create Notion pages for literature items
5. Create RAG chunk pages in Notion
6. Upsert vectors to Pinecone with rich metadata

### Ingest Emails/Notes from Notion

```bash
python scripts/ingest_email.py
```

This processes all emails/notes in your Notion Email DB that have `Embedding Status = "Not Embedded"`.

### Ingest Experiments from Notion

```bash
python scripts/ingest_experiment.py --experiment-page-id YOUR_EXPERIMENT_PAGE_ID
```

Ingests a single experiment page from the Notion "Experiments" database into Pinecone with full metadata and signature awareness.

### Ingest Datasets from Notion

```bash
python scripts/ingest_dataset.py --dataset-page-id YOUR_DATASET_PAGE_ID
```

Ingests a single dataset page from the Notion "Experimental Data Assets" database into Pinecone.

### Harvest Metabolomics Workbench Studies

The MW harvester fetches lipidomics studies from Metabolomics Workbench and creates/updates Dataset pages in Notion.

#### Discover Studies by Keyword

```bash
python scripts/harvest_mw_studies.py \
  --search-keyword ALS \
  --search-keyword amyotrophic \
  --max-search-results 50 \
  --dry-run
```

Searches MW for studies matching keywords in title/summary/disease fields. Use `--dry-run` to preview results without creating pages.

#### Harvest Specific Studies

```bash
python scripts/harvest_mw_studies.py \
  --study-id ST004396 \
  --create-notion \
  --ingest
```

This will:
1. Fetch study metadata from MW
2. Fetch mwTab data
3. Create/update a Dataset page in Notion
4. Embed mwTab content as code blocks
5. Optionally trigger dataset ingestion into Pinecone (`--ingest`)

#### Filter by Lipid Content (Future)

```bash
python scripts/harvest_mw_studies.py \
  --search-keyword ALS \
  --lipid-filter Ceramide \
  --create-notion \
  --ingest
```

Note: Lipid filtering requires specific RefMet names. Generic names like "Ceramide" may not match MW's metabolite database.

### Convert mwTab to CSV

```bash
python scripts/convert_mwtab_to_csv.py --study-id ST004396
```

Downloads mwTab content from MW and extracts tabular metabolite data to CSV format. Saves to `data/mwtab/<study_id>.csv` by default.

```bash
# Custom output directory
python scripts/convert_mwtab_to_csv.py \
  --study-id ST004396 \
  --output-dir data/mwtab_export
```

**Note**: Some studies may have JSON-only mwTab format without tab-separated sections. The script will report if no tabular data is found.

### Scan CSV for Ceramides

```bash
python scripts/scan_ceramides_in_mwtab_csv.py --study-id ST004396
```

Scans a CSV file (from `convert_mwtab_to_csv.py`) for ceramide-like metabolites and produces a summary of hits.

```bash
# Scan specific CSV file
python scripts/scan_ceramides_in_mwtab_csv.py \
  --csv-path data/mwtab/ST004396.csv

# Generate a report file
python scripts/scan_ceramides_in_mwtab_csv.py \
  --study-id ST004396 \
  --output-report data/mwtab/ST004396_ceramide_report.txt
```

### Query the RAG System

```bash
# Single source type
python scripts/rag_query.py \
  --query "ceramide dysregulation and neurodegeneration" \
  --top-k 5 \
  --source-type Literature \
  --disease ALS \
  --signature "ALS-CSF-Core-6Ceramides"

# Multiple source types
python scripts/rag_query.py \
  --query "ALS ceramide signature" \
  --top-k 10 \
  --source-type Literature Experiment Dataset \
  --signature "ALS-CSF-Core-6Ceramides"
```

Options:
- `--query`: Your search query
- `--top-k`: Number of results (default: 10)
- `--source-type`: Filter by source type(s). Can specify multiple: `Literature`, `Email`, `Experiment`, `Dataset` (default: `Literature`)
- `--disease`: Filter by disease (e.g., "ALS", "AD")
- `--target`: Filter by molecular target (e.g., "SPTLC1")
- `--lipid`: Filter by lipid species
- `--signature`: Filter by lipid signature (Short ID)
- `--tag`: Filter by tag
- `--show-context`: Show context chunks used
- `--raw-json`: Output raw JSON response
- `--no-answer`: Skip answer synthesis (faster)

### Classify Literature Metadata

```bash
python scripts/classify_literature_metadata.py
```

Uses OpenAI to classify papers into semantic metadata (diseases, targets, lipid signatures, etc.).

### Maintenance Scripts

- `scripts/cleanup_deleted_items.py` - Remove deleted Zotero items
- `scripts/sync_collection_state.py` - Sync collection state
- `scripts/verify_rag_metadata.py` - Verify metadata consistency
- `scripts/reset_all.py` - Full system reset (use with caution!)

## Features

### Idempotent Ingestion
- Zotero attachments: Uses MD5 hash to skip unchanged files
- Zotero notes: Uses content hash to skip unchanged notes
- Safe to re-run ingestion scripts

### Rich Metadata
- **Document-level**: title, journal, DOI, year, importance
- **Semantic**: diseases, molecular targets, modality, stage
- **Lipid-specific**: lipid species, lipid signatures, phenotype axes
- **Source tracking**: Zotero item keys, Notion page IDs
- **Dataset-specific**: MW Study ID, data origin, sample type, matrix, model systems

### Multi-Source RAG
- Query across multiple source types: Literature, Email, Experiments, Datasets
- Unified filtering by disease, target, lipid, signature across all sources
- Provenance tracking showing which sources contributed to results

### Structured Logging
All modules use consistent logging prefixes:
- `[INGEST][ZOTERO]` - Zotero ingestion
- `[INGEST][EMAIL]` - Email/note ingestion
- `[INGEST][DATASET]` - Dataset ingestion
- `[MW]` - Metabolomics Workbench operations
- `[NOTION]` - Notion API operations
- `[PINECONE]` - Pinecone operations
- `[RAG]` - Query/RAG operations

## Development

### Running Tests

```bash
pytest amprenta_rag/tests -v
```

### Code Structure

The codebase follows a modular architecture:

- **Clients** (`clients/`): Thin wrappers around external APIs
- **Ingestion** (`ingestion/`): Data ingestion pipelines
- **Query** (`query/`): RAG query engine (separated into low-level and high-level)
- **Metadata** (`metadata/`): Semantic metadata classification
- **Maintenance** (`maintenance/`): Cleanup and sync utilities

### Refactoring Status

All major refactoring phases are complete:
- ✅ Phase 1-2: Ingestion pipeline modularization
- ✅ Phase 3: Query engine refactoring
- ✅ Phase 4: Logging and error handling
- ✅ Phase 5: Code hygiene and cleanup

See `amprenta_rag/amprenta_rag_roadmap.md` for details.

## Troubleshooting

### Common Issues

1. **ModuleNotFoundError**: Make sure you're in a virtual environment with all dependencies installed.

2. **API Key Errors**: Verify environment variables are set correctly.

3. **Notion API Errors**: Check that database IDs in `config.py` are correct (no dashes in IDs).

4. **Pinecone Namespace Issues**: Ensure you're using the correct namespace (configured in `config.py`).

## License

[Your license here]

## Contributing

[Your contributing guidelines here]

