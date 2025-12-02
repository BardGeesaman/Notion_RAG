# Amprenta RAG System

A modular RAG (Retrieval-Augmented Generation) system for scientific literature and email/note ingestion, designed specifically for ceramide/sphingolipid neurodegeneration research.

## Overview

This system ingests content from **Zotero** and **Notion** (emails/notes), embeds it using OpenAI, stores vectors in **Pinecone**, and provides semantic search with answer synthesis. It's designed to help researchers quickly find and understand relevant scientific literature and internal notes.

## Architecture

- **Data Sources**: Zotero (literature PDFs/notes), Notion (emails/notes)
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
├── ingest_collection.py    # Ingest Zotero collection
├── ingest_email.py         # Ingest Notion emails/notes
├── rag_query.py            # Query the RAG system
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

Set environment variables for API keys:

```bash
export OPENAI_API_KEY="your_openai_key"
export PINECONE_API_KEY="your_pinecone_key"
export NOTION_API_KEY="your_notion_key"
export ZOTERO_API_KEY="your_zotero_key"
```

Or create a `.env` file (already in `.gitignore`):

```bash
OPENAI_API_KEY=your_openai_key
PINECONE_API_KEY=your_pinecone_key
NOTION_API_KEY=your_notion_key
ZOTERO_API_KEY=your_zotero_key
```

Then load it with `python-dotenv` (optional):

```python
from dotenv import load_dotenv
load_dotenv()
```

**Note**: The `config.py` file also contains hardcoded database IDs and other constants that may need to be updated for your workspace.

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

### Query the RAG System

```bash
python scripts/rag_query.py \
  --query "ceramide dysregulation and neurodegeneration" \
  --top-k 5 \
  --source-type Literature \
  --disease ALS \
  --signature "ALS-CSF-Core-6Ceramides"
```

Options:
- `--query`: Your search query
- `--top-k`: Number of results (default: 10)
- `--source-type`: Filter by source ("Literature", "Email", "Note")
- `--disease`: Filter by disease (e.g., "ALS", "AD")
- `--target`: Filter by molecular target (e.g., "SPTLC1")
- `--lipid`: Filter by lipid species
- `--signature`: Filter by lipid signature
- `--tag`: Filter by tag
- `--show-context`: Show context chunks used
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

### Structured Logging
All modules use consistent logging prefixes:
- `[INGEST][ZOTERO]` - Zotero ingestion
- `[INGEST][EMAIL]` - Email/note ingestion
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

