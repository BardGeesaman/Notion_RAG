# Quick Start Guide

Get up and running with the Amprenta RAG System in minutes.

## Prerequisites

- Python 3.10 or higher
- PostgreSQL 13 or higher
- API keys for:
  - OpenAI (for embeddings and LLM)
  - Pinecone (vector database)
  - Zotero (optional, for literature ingestion)

## Installation

### 1. Clone and Setup

```bash
# Clone the repository
git clone https://github.com/BardGeesaman/Notion_RAG.git
cd Notion_RAG

# Create virtual environment
python -m venv venv
source venv/bin/activate  # On Windows: venv\Scripts\activate

# Install dependencies
pip install -r requirements.txt
```

### 2. Configure Environment

Create a `.env` file in the project root:

```bash
# Required API Keys
OPENAI_API_KEY=your_openai_key_here
PINECONE_API_KEY=your_pinecone_key_here

# Postgres Database
POSTGRES_HOST=localhost
POSTGRES_PORT=5432
POSTGRES_USER=your_username
POSTGRES_PASSWORD=your_password
POSTGRES_DB=amprenta

# Optional
ZOTERO_API_KEY=your_zotero_key_here
ZOTERO_LIBRARY_ID=your_library_id
ZOTERO_LIBRARY_TYPE=user  # or "group"
```

### 3. Setup Postgres Database

```bash
# Create the database
createdb amprenta

# Run migrations
alembic upgrade head
```

### 4. Verify Configuration

```bash
python scripts/validate_configuration.py
```

This will check:
- All API keys are present
- Postgres database is accessible
- Pinecone connection works
- Configuration is valid

## Your First Ingestion

### Option 1: Ingest a Lipidomics Dataset

```bash
# Create a simple test dataset file
echo "Species,Intensity
Cer(d18:1/16:0),1000
Cer(d18:1/18:0),1500
SM(d18:1/16:0),2000" > test_lipidomics.csv

# Ingest it
python scripts/ingest_lipidomics.py --file test_lipidomics.csv --create-page
```

### Option 2: Discover and Harvest a Public Study

```bash
# Discover studies from Metabolomics Workbench
python scripts/discover_omics_studies.py --repository MW --keyword "ceramide" --max-results 5

# Harvest a specific study
python scripts/harvest_repository_study.py --repository MW --study-id ST004396 --ingest
```

### Option 3: Ingest a Signature

```bash
# Create a signature file
cat > test_signature.tsv << EOF
feature_type	feature_name	direction	weight
lipid	Cer(d18:1/16:0)	↑	1.0
lipid	Cer(d18:1/18:0)	↑	1.0
lipid	SM(d18:1/16:0)	↓	1.0
EOF

# Ingest the signature
python scripts/ingest_signature.py --file test_signature.tsv
```

## Your First Query

```bash
# Simple query
python scripts/rag_query.py --query "What is the role of ceramides in neurodegeneration?"

# Query with filters
python scripts/rag_query.py \
  --query "ALS ceramide signature" \
  --source-type Dataset \
  --top-k 10

# Cross-omics reasoning
python scripts/rag_query.py \
  --cross-omics-program YOUR_PROGRAM_PAGE_ID
```

## Common Workflows

### 1. Batch Ingest Multiple Files

```bash
# Ingest all files in a directory (auto-detect omics type)
python scripts/batch_ingest_omics.py \
  --directory /path/to/omics/data \
  --parallel \
  --max-workers 4

# Export results
python scripts/batch_ingest_omics.py \
  --directory /path/to/data \
  --export-results results.json
```

### 2. Discover Signatures from Datasets

```bash
# Discover patterns across all datasets
python scripts/discover_signatures.py \
  --all-datasets \
  --min-support 3 \
  --min-confidence 0.7 \
  --output discovered_signatures.json

# Ingest discovered signatures
python scripts/discover_signatures.py \
  --all-datasets \
  --ingest
```

### 3. Generate Evidence Reports

```bash
# Generate report for a program
python scripts/generate_evidence_report.py \
  --program-id YOUR_PROGRAM_ID \
  --output program_report.md

# Generate report for a dataset
python scripts/generate_evidence_report.py \
  --dataset-id YOUR_DATASET_ID \
  --output dataset_report.md
```

## Health Check

Run a comprehensive health check:

```bash
python scripts/health_check.py
```

This checks:
- ✅ Postgres database connectivity
- ✅ Pinecone connectivity
- ✅ OpenAI API access
- ✅ Configuration validity

## Next Steps

1. **Learn More**: Read the [Comprehensive User Guide](USER_GUIDE.md)
2. **Understand Architecture**: See [Architecture Overview](ARCHITECTURE.md)
3. **Troubleshoot Issues**: Check [Troubleshooting Guide](TROUBLESHOOTING.md)
4. **API Reference**: See [API Reference](API_REFERENCE.md)

## Getting Help

- Check the [Troubleshooting Guide](TROUBLESHOOTING.md) for common issues
- Review [Usage Examples](USAGE_EXAMPLES.md) for advanced workflows
- See [Configuration Guide](CONFIGURATION.md) for detailed setup options

## System Status

Your system is ready when:
- ✅ Configuration validation passes
- ✅ Health check shows all green
- ✅ You can successfully ingest at least one file
- ✅ You can run at least one query

Congratulations! You're ready to use the Amprenta RAG System.

