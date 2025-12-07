# Amprenta RAG System

A comprehensive multi-omics knowledge management platform with RAG (Retrieval-Augmented Generation) capabilities for scientific research, designed specifically for ceramide/sphingolipid neurodegeneration research and extensible to any omics domain.

[![Python 3.10+](https://img.shields.io/badge/python-3.10+-blue.svg)](https://www.python.org/downloads/)
[![License](https://img.shields.io/badge/license-MIT-green.svg)](LICENSE)

## Overview

The Amprenta RAG System is a production-ready platform that:
- **Ingests multi-omics data** from multiple sources (lipidomics, metabolomics, proteomics, transcriptomics)
- **Manages multi-omics signatures** with automatic discovery and scoring
- **Provides semantic search** across all data sources with RAG
- **Generates evidence-based reports** using cross-omics reasoning
- **Discovers patterns automatically** using statistical analysis
- **Integrates with public repositories** (Metabolomics Workbench, GEO, PRIDE, MetaboLights)

## Quick Start

```bash
# Clone and setup
git clone https://github.com/BardGeesaman/Notion_RAG.git
cd Notion_RAG
python -m venv venv && source venv/bin/activate
pip install -r requirements.txt

# Configure (see docs/CONFIGURATION.md)
cp .env.example .env
# Edit .env with your API keys

# Verify setup
python scripts/validate_configuration.py

# Your first ingestion
python scripts/ingest_lipidomics.py --file data.csv --create-page
```

📖 **See [Quick Start Guide](docs/QUICK_START.md) for detailed setup instructions**

## Key Features

### 🔬 Multi-Omics Ingestion

- **Lipidomics**: Automatic species normalization, canonical format conversion
- **Metabolomics**: Metabolite name normalization and linking
- **Proteomics**: Protein identifier normalization (UniProt, gene symbols)
- **Transcriptomics**: Gene identifier normalization (Ensembl, gene symbols)
- **Batch Processing**: Auto-detect omics type from filename and content, 4x speedup with parallel workers
- **Postgres-First Architecture**: All omics pipelines use Postgres as primary database for fast, scalable ingestion
- **Automatic Linking**: Intelligent program/experiment linking based on metadata confidence scoring

### 📊 Signature Management

- **Multi-Omics Signatures**: Support for genes, proteins, metabolites, lipids
- **Automatic Discovery**: Statistical pattern detection across datasets
- **Signature Scoring**: Match datasets against signatures with direction consistency
- **Feature Linking**: Automatic linking to canonical feature pages

### 🔍 Advanced RAG Queries

- **Semantic Search**: Query across all data sources (Literature, Experiments, Datasets)
- **Advanced Filtering**: Filter by disease, target, signature, omics type
- **Cross-Omics Reasoning**: LLM-powered multi-omics evidence summaries
- **Signature Similarity**: Find matching datasets and signatures

### 📈 Analysis & Reports

- **Evidence Reports**: Automated cross-omics evidence summaries
- **Dataset Comparison**: Jaccard similarity, shared/differential features
- **Program Signature Maps**: Program × Signature matrices
- **Pathway Enrichment**: KEGG/Reactome pathway analysis

### 🚀 Performance & Production

- **Feature Caching**: LRU cache with persistence for 10-100x performance gains (see [Feature Caching Guide](docs/FEATURE_CACHING.md))
- **Parallel Processing**: Configurable worker pools for batch operations
- **Production Hardening**: Error handling, retry logic, circuit breakers
- **Health Monitoring**: Comprehensive health checks and performance metrics
- **REST API**: FastAPI endpoints for programmatic access (compounds, screening, multi-omics data)

## Architecture

```
┌─────────────────────────────────────────────────────────────┐
│                    Data Sources                              │
├─────────────────────────────────────────────────────────────┤
│  Zotero │ Notion │ Public Repos (MW, GEO, PRIDE, etc.)     │
└──────────────────────┬──────────────────────────────────────┘
                       │
                       ▼
┌─────────────────────────────────────────────────────────────┐
│              Ingestion Pipelines                             │
├─────────────────────────────────────────────────────────────┤
│  Multi-Omics │ Signatures │ Literature │ Experiments        │
│  • Normalization                                            │
│  • Feature Extraction                                       │
│  • Signature Matching                                       │
└──────────────────────┬──────────────────────────────────────┘
                       │
                       ▼
┌─────────────────────────────────────────────────────────────┐
│              Knowledge Graph (Notion)                        │
├─────────────────────────────────────────────────────────────┤
│  Programs │ Experiments │ Datasets │ Signatures │ Features  │
└──────────────────────┬──────────────────────────────────────┘
                       │
                       ▼
┌─────────────────────────────────────────────────────────────┐
│              Vector Store (Pinecone)                         │
├─────────────────────────────────────────────────────────────┤
│  Embedded Chunks │ Rich Metadata │ Multi-source tracking    │
└──────────────────────┬──────────────────────────────────────┘
                       │
                       ▼
┌─────────────────────────────────────────────────────────────┐
│              RAG Query Engine                                │
├─────────────────────────────────────────────────────────────┤
│  Semantic Search │ Cross-Omics Reasoning │ Evidence Reports │
└─────────────────────────────────────────────────────────────┘
```

**Key Components**:
- **Data Sources**: Zotero (literature), Notion (internal), Public repositories
- **Vector Store**: Pinecone with OpenAI embeddings
- **Knowledge Graph**: Notion databases for structured metadata
- **Query Engine**: RAG with cross-omics reasoning

📖 **See [Architecture Overview](docs/ARCHITECTURE.md) for detailed system design**

## Documentation

| Document | Description |
|----------|-------------|
| [📚 Quick Start Guide](docs/QUICK_START.md) | Get up and running in minutes |
| [📖 Comprehensive User Guide](docs/USER_GUIDE.md) | Complete feature documentation |
| [🔧 Troubleshooting Guide](docs/TROUBLESHOOTING.md) | Common issues and solutions |
| [⚙️ Configuration Guide](docs/CONFIGURATION.md) | Detailed configuration options |
| [🏗️ Architecture Overview](docs/ARCHITECTURE.md) | System design and components |
| [📋 API Reference](docs/API_REFERENCE.md) | REST API & module documentation |
| [💡 Usage Examples](docs/USAGE_EXAMPLES.md) | Practical code examples |
| [⚡ Feature Caching Guide](docs/FEATURE_CACHING.md) | Performance optimization with caching |
| [🔗 Auto-Linking Guide](docs/AUTO_LINKING.md) | Automatic program/experiment linking |
| [🗄️ Notion Database Setup](docs/NOTION_DATABASE_SETUP.md) | Database configuration guide |
| [🛡️ Production Hardening](docs/PRODUCTION_HARDENING.md) | Production deployment guide |
| [📧 Gmail Setup](docs/setup/GMAIL_SETUP.md) | Gmail API integration guide |
| [🔐 OAuth Setup](docs/setup/OAUTH_SETUP.md) | OAuth2 authentication guide |

## Usage Examples

### Ingest Multi-Omics Data

```bash
# Single dataset ingestion
python scripts/ingest_lipidomics.py --file data.csv --create-page
python scripts/ingest_metabolomics.py --file data.csv --create-page

# Batch ingestion (auto-detect type, 4x faster with parallel processing)
python scripts/batch_ingest_omics.py --directory /path/to/data --create-pages

# Batch with program linking
python scripts/batch_ingest_omics.py \
    --directory /path/to/data \
    --create-pages \
    --program-id PROGRAM_UUID \
    --parallel 4

# Override type detection
python scripts/batch_ingest_omics.py \
    --directory /path/to/data \
    --type metabolomics \
    --create-pages
```

### Discover and Manage Signatures

```bash
# Discover patterns from datasets
python scripts/discover_signatures.py --all-datasets --min-confidence 0.7

# Ingest a signature
python scripts/ingest_signature.py --file signature.tsv

# Score dataset against signatures
python scripts/score_signature.py --dataset-id DATASET_ID
```

### Query with RAG

```bash
# Simple query
python scripts/rag_query.py --query "ceramide dysregulation in ALS"

# Cross-omics reasoning
python scripts/rag_query.py --cross-omics-program PROGRAM_ID

# Signature similarity
python scripts/rag_query.py --signature-score DATASET_ID
```

### Generate Reports

```bash
# Evidence report
python scripts/generate_evidence_report.py --program-id PROGRAM_ID

# Dataset comparison
python scripts/compare_datasets.py --dataset-id-1 ID1 --dataset-id-2 ID2

# Program signature map
python scripts/generate_program_signature_map.py --program-id PROGRAM_ID
```

### REST API Access

```bash
# Start FastAPI server
uvicorn amprenta_rag.api.main:app --reload --port 8000

# Access chemistry compounds
curl http://localhost:8000/api/v1/compounds/

# Get HTS screening campaigns
curl http://localhost:8000/api/v1/screening/campaigns

# Interactive API docs
open http://localhost:8000/docs
```

📖 **See [Usage Examples](docs/USAGE_EXAMPLES.md) for more examples**

## Project Structure

```
amprenta_rag/
├── clients/              # API client wrappers (OpenAI, Pinecone, Notion)
├── config.py             # Configuration management
├── ingestion/            # Data ingestion pipelines
│   ├── lipidomics/      # Lipidomics ingestion
│   ├── metabolomics/    # Metabolomics ingestion
│   ├── proteomics/      # Proteomics ingestion
│   ├── transcriptomics/ # Transcriptomics ingestion
│   ├── signature_ingestion.py
│   └── ...
├── query/                # RAG query engine
│   ├── rag_engine.py
│   └── cross_omics/     # Cross-omics reasoning
├── signatures/           # Signature management
│   ├── discovery.py     # Automatic signature discovery
│   ├── signature_loader.py
│   └── ...
├── analysis/             # Analysis tools
│   ├── dataset_comparison.py
│   ├── pathway_analysis.py
│   └── ...
└── utils/                # Utilities
    ├── error_handling.py
    ├── performance.py
    └── ...

scripts/                  # Command-line scripts
├── ingest_lipidomics.py
├── ingest_metabolomics.py
├── batch_ingest_omics.py
├── discover_signatures.py
├── rag_query.py
└── ...
```

## Setup

### Prerequisites

- Python 3.10 or higher
- API keys for:
  - OpenAI (embeddings and LLM)
  - Pinecone (vector database)
  - Notion (knowledge graph)
  - Zotero (optional, for literature)

### Installation

```bash
# Clone repository
git clone https://github.com/BardGeesaman/Notion_RAG.git
cd Notion_RAG

# Create virtual environment
python -m venv venv
source venv/bin/activate  # On Windows: venv\Scripts\activate

# Install dependencies
pip install -r requirements.txt
```

### Configuration

1. **Copy environment template**:
   ```bash
   cp .env.example .env
   ```

2. **Fill in API keys** (see `.env.example` for all required keys)

3. **Set up Notion databases** (see [Notion Database Setup](docs/NOTION_DATABASE_SETUP.md))

4. **Verify configuration**:
   ```bash
   python scripts/validate_configuration.py
   ```

📖 **See [Configuration Guide](docs/CONFIGURATION.md) for detailed setup**

## Features in Detail

### Multi-Omics Support

- **Automatic Type Detection**: Detects omics type from file names and headers
- **Feature Normalization**: Converts various formats to canonical forms
- **Feature Linking**: Automatic creation/linking to feature pages in Notion
- **Signature Scoring**: Scores datasets against multi-omics signatures

### Signature Discovery

- **Pattern Detection**: Statistical analysis of feature co-occurrence
- **Direction Consistency**: Verifies directional consistency across datasets
- **Clustering**: Groups features into signature candidates
- **Confidence Scoring**: Ranks candidates by statistical significance

### Performance Optimization

- **Feature Caching**: LRU cache with file persistence (10-100x faster) - see [Feature Caching Guide](docs/FEATURE_CACHING.md)
- **Parallel Processing**: Configurable worker pools for batch operations
- **Progress Tracking**: Real-time progress bars with tqdm
- **Error Aggregation**: Comprehensive error reporting

**Cache Management**:

```bash
# Warm cache before batch operations
python scripts/warm_feature_cache.py --all-datasets

# Monitor cache performance
python scripts/manage_feature_cache.py --stats

# Export/import cache for session continuity
python scripts/manage_feature_cache.py --export cache_backup.json
python scripts/manage_feature_cache.py --import cache_backup.json
```

### Production Ready

- **Error Handling**: Retry logic with exponential backoff
- **Circuit Breakers**: Prevents cascade failures
- **Health Checks**: Comprehensive system health monitoring
- **Configuration Validation**: Startup validation of all components

## Development

### Running Tests

```bash
# Run all tests
pytest amprenta_rag/tests -v

# Run specific test
pytest amprenta_rag/tests/test_signature_scoring.py -v

# Run with coverage
pytest --cov=amprenta_rag --cov-report=html
```

### Code Quality

The codebase follows best practices:
- Modular architecture
- Comprehensive error handling
- Structured logging
- Type hints throughout
- Comprehensive documentation

## Troubleshooting

Common issues and quick fixes:

| Issue | Solution |
|-------|----------|
| ModuleNotFoundError | `pip install -r requirements.txt` |
| API Key Missing | Check `.env` file exists and contains keys |
| Database Access Error | Run `python scripts/verify_notion_setup.py` |
| No Query Results | Verify data is ingested: check Pinecone index |
| Performance Issues | Enable feature caching and parallel processing |

📖 **See [Troubleshooting Guide](docs/TROUBLESHOOTING.md) for detailed solutions**

## Roadmap

### Completed ✅

- Multi-omics ingestion (lipidomics, metabolomics, proteomics, transcriptomics)
- Multi-omics signature support
- Automated signature discovery
- Cross-omics reasoning
- Batch ingestion framework
- Performance optimization (feature caching)
- Production hardening

### In Progress 🚧

- Enhanced cross-omics reasoning
- Visualization dashboards
- Advanced analytics

### Planned 📋

- Architecture evolution (Postgres + FastAPI)
- Enhanced batch ingestion features
- Extended testing suite

📖 **See [Strategic Roadmap](context/UNIFIED_STRATEGIC_ROADMAP.md) for full roadmap**

## Contributing

Contributions are welcome! Please:
1. Check existing issues and documentation
2. Follow code style and conventions
3. Add tests for new features
4. Update documentation as needed

## License

[Your license here]

## Support

- 📖 **Documentation**: See `docs/` directory
- 🐛 **Issues**: Report on GitHub
- 💬 **Questions**: Check [Troubleshooting Guide](docs/TROUBLESHOOTING.md)

## Acknowledgments

Built for ceramide/sphingolipid neurodegeneration research, with extensibility to any omics domain.

---

**Last Updated**: 2025-01-XX | **Version**: 2.0.0
