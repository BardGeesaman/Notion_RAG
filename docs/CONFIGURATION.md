# Configuration Guide

**Last Updated**: 2025-12-03

This document describes how to configure the Amprenta RAG platform.

---

## Environment Variables

All configuration is managed via environment variables in a `.env` file in the project root.

### Required Configuration

#### Notion API

```bash
NOTION_API_KEY=secret_xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
```

Get your Notion integration token from:
1. Go to https://www.notion.so/my-integrations
2. Create a new integration
3. Copy the "Internal Integration Token"

#### Notion Database IDs

```bash
# Core databases
NOTION_DATASET_DB_ID=xxxxxxxx-xxxx-xxxx-xxxx-xxxxxxxxxxxx
NOTION_SIGNATURE_DB_ID=xxxxxxxx-xxxx-xxxx-xxxx-xxxxxxxxxxxx
NOTION_SIGNATURE_COMPONENT_DB_ID=xxxxxxxx-xxxx-xxxx-xxxx-xxxxxxxxxxxx
NOTION_LIPID_SPECIES_DB_ID=xxxxxxxx-xxxx-xxxx-xxxx-xxxxxxxxxxxx

# Feature databases
NOTION_GENE_FEATURES_DB_ID=xxxxxxxx-xxxx-xxxx-xxxx-xxxxxxxxxxxx
NOTION_PROTEIN_FEATURES_DB_ID=xxxxxxxx-xxxx-xxxx-xxxx-xxxxxxxxxxxx
NOTION_METABOLITE_FEATURES_DB_ID=xxxxxxxx-xxxx-xxxx-xxxx-xxxxxxxxxxxx

# Optional databases
NOTION_PROGRAMS_DB_ID=xxxxxxxx-xxxx-xxxx-xxxx-xxxxxxxxxxxx
NOTION_EXPERIMENTS_DB_ID=xxxxxxxx-xxxx-xxxx-xxxx-xxxxxxxxxxxx
```

**How to find Database IDs**:
1. Open the database in Notion
2. Copy the URL: `https://www.notion.so/workspace/xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx?v=...`
3. Extract the 32-character hex ID (with dashes)

#### Pinecone

```bash
PINECONE_API_KEY=xxxxxxxx-xxxx-xxxx-xxxx-xxxxxxxxxxxx
PINECONE_INDEX_NAME=amprenta-rag
PINECONE_ENVIRONMENT=us-east-1-aws
```

Get your Pinecone API key from:
1. Go to https://app.pinecone.io/
2. Navigate to API Keys
3. Copy your API key

#### OpenAI

```bash
OPENAI_API_KEY=sk-xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
```

Get your OpenAI API key from:
1. Go to https://platform.openai.com/api-keys
2. Create a new secret key
3. Copy the key (starts with `sk-`)

---

### Optional Configuration

#### Signature Ingestion

```bash
# Directory containing signature TSV/CSV files
SIGNATURES_DIR=data/signatures

# Signature overlap threshold (0.0 to 1.0)
SIGNATURE_OVERLAP_THRESHOLD=0.3

# Enable/disable signature scoring
ENABLE_SIGNATURE_SCORING=True

# Enable/disable lipid mapping
ENABLE_LIPID_MAPPING=True
```

#### Feature Cache

The feature cache provides 10-100x performance improvements for signature scoring operations by caching dataset features in memory with LRU eviction.

```bash
# Enable/disable feature cache (default: true)
FEATURE_CACHE_ENABLED=true

# Cache TTL in seconds (default: 3600 = 1 hour)
FEATURE_CACHE_TTL=3600

# Maximum number of cached datasets (default: 1000)
FEATURE_CACHE_MAX_SIZE=1000

# Enable disk-based persistence (default: true)
FEATURE_CACHE_ENABLE_PERSISTENCE=true

# Custom cache directory (optional, defaults to temp directory)
FEATURE_CACHE_DIR=/path/to/cache

# Parallel workers for cache warming (default: 5)
FEATURE_CACHE_PARALLEL_WORKERS=5
```

**Performance Impact**:
- Cold cache: 100-500ms per dataset (Postgres queries)
- Warm cache: 1-5ms per dataset (memory access)
- Eliminates 90-99% of Postgres queries during signature scoring

**Cache Management**:

```bash
# Warm cache before batch operations
python scripts/warm_feature_cache.py --all-datasets

# Monitor cache statistics
python scripts/manage_feature_cache.py --stats

# Export cache for backup
python scripts/manage_feature_cache.py --export cache_backup.json
```

📖 **See [Feature Caching Guide](FEATURE_CACHING.md) for detailed usage**

#### Logging

```bash
# Log level: DEBUG, INFO, WARNING, ERROR
LOG_LEVEL=INFO
```

#### OpenAI Models

```bash
# Chat model (default: gpt-4o-mini)
OPENAI_CHAT_MODEL=gpt-4o-mini

# Embedding model (default: text-embedding-3-small)
OPENAI_EMBEDDING_MODEL=text-embedding-3-small
```

---

## Example .env File

```bash
# Notion
NOTION_API_KEY=secret_xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
NOTION_DATASET_DB_ID=xxxxxxxx-xxxx-xxxx-xxxx-xxxxxxxxxxxx
NOTION_SIGNATURE_DB_ID=xxxxxxxx-xxxx-xxxx-xxxx-xxxxxxxxxxxx
NOTION_SIGNATURE_COMPONENT_DB_ID=xxxxxxxx-xxxx-xxxx-xxxx-xxxxxxxxxxxx
NOTION_LIPID_SPECIES_DB_ID=xxxxxxxx-xxxx-xxxx-xxxx-xxxxxxxxxxxx
NOTION_GENE_FEATURES_DB_ID=xxxxxxxx-xxxx-xxxx-xxxx-xxxxxxxxxxxx
NOTION_PROTEIN_FEATURES_DB_ID=xxxxxxxx-xxxx-xxxx-xxxx-xxxxxxxxxxxx
NOTION_METABOLITE_FEATURES_DB_ID=xxxxxxxx-xxxx-xxxx-xxxx-xxxxxxxxxxxx
NOTION_PROGRAMS_DB_ID=xxxxxxxx-xxxx-xxxx-xxxx-xxxxxxxxxxxx
NOTION_EXPERIMENTS_DB_ID=xxxxxxxx-xxxx-xxxx-xxxx-xxxxxxxxxxxx

# Pinecone
PINECONE_API_KEY=xxxxxxxx-xxxx-xxxx-xxxx-xxxxxxxxxxxx
PINECONE_INDEX_NAME=amprenta-rag
PINECONE_ENVIRONMENT=us-east-1-aws

# OpenAI
OPENAI_API_KEY=sk-xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
OPENAI_CHAT_MODEL=gpt-4o-mini
OPENAI_EMBEDDING_MODEL=text-embedding-3-small

# Signature Configuration
SIGNATURES_DIR=data/signatures
SIGNATURE_OVERLAP_THRESHOLD=0.3
ENABLE_SIGNATURE_SCORING=True
ENABLE_LIPID_MAPPING=True

# Feature Cache (10-100x performance boost)
FEATURE_CACHE_ENABLED=true
FEATURE_CACHE_TTL=3600
FEATURE_CACHE_MAX_SIZE=1000
FEATURE_CACHE_ENABLE_PERSISTENCE=true
FEATURE_CACHE_PARALLEL_WORKERS=5

# Logging
LOG_LEVEL=INFO
```

---

## Notion Database Setup

### Required Databases

#### 1. Experimental Data Assets

**Properties**:
- `Title` (title): Dataset name
- `Omics Type` (select): Lipidomics, Metabolomics, Proteomics, Transcriptomics
- `Summary` (rich_text): Dataset description
- `Signature Match Score` (number): Highest signature match score
- `Related Signature(s)` (relation): Links to Lipid Signatures
- `Related Programs` (relation): Links to Programs (optional)
- `Related Experiments` (relation): Links to Experiments (optional)
- `Embedding IDs` (rich_text): Pinecone embedding IDs
- `Last Embedded` (date): Last embedding timestamp

#### 2. Lipid Signatures

**Properties**:
- `Name` (title): Signature name
- `Short ID` (rich_text): Unique short identifier
- `Type` (select): Literature-derived, Consortium, Open Dataset
- `Modalities` (multi_select): Gene, Protein, Metabolite, Lipid
- `Description` (rich_text): Signature description
- `Related Components` (relation): Links to Lipid Signature Components
- `Related Datasets` (relation): Links to Experimental Data Assets

#### 3. Lipid Signature Components

**Properties**:
- `Component Name` (title): Feature name
- `Feature Type` (select): gene, protein, metabolite, lipid
- `Direction` (select): ↑ (up), ↓ (down)
- `Weight` (number): Component weight
- `Related Signature` (relation): Links to Lipid Signatures
- `Related Feature` (relation): Links to Gene/Protein/Metabolite/Lipid Species DB

#### 4. Lipid Species

**Properties**:
- `Species Name` (title): Canonical lipid species name
- `Synonyms` (rich_text): Alternative names
- `Class` (select): Cer, SM, HexCer, etc.
- `Related Datasets` (relation): Links to Experimental Data Assets

#### 5. Gene Features

**Properties**:
- `Gene Name` (title): Gene symbol or Ensembl ID
- `Related Datasets` (relation): Links to Experimental Data Assets

#### 6. Protein Features

**Properties**:
- `Protein Name` (title): UniProt ID or canonical name
- `Related Datasets` (relation): Links to Experimental Data Assets

#### 7. Metabolite Features

**Properties**:
- `Metabolite Name` (title): Normalized metabolite name
- `Related Datasets` (relation): Links to Experimental Data Assets

---

### Optional Databases

#### Programs

**Properties**:
- `Program` (title): Program name
- `Description` (rich_text): Program description
- `Experiments` (relation): Links to Experiments
- `Related Datasets` (relation): Links to Experimental Data Assets

#### Experiments

**Properties**:
- `Title` (title): Experiment name
- `Type` (select): Experiment type
- `Disease` (multi_select): Disease context
- `Matrix` (multi_select): Sample matrix (CSF, Plasma, Serum, etc.)
- `Model Systems` (multi_select): Model systems
- `Related Programs` (relation): Links to Programs
- `Related Datasets` (relation): Links to Experimental Data Assets
- `Readout Signatures` (relation): Links to Lipid Signatures

---

## Pinecone Index Setup

### Create Index

```python
import pinecone

pinecone.init(api_key="your-api-key", environment="us-east-1-aws")

pinecone.create_index(
    name="amprenta-rag",
    dimension=1536,  # text-embedding-3-small dimension
    metric="cosine"
)
```

### Index Configuration

- **Dimension**: 1536 (for `text-embedding-3-small`)
- **Metric**: cosine
- **Pod Type**: s1.x1 (or higher for production)

---

## Verification

### Test Configuration

```python
from amprenta_rag.config import get_config

cfg = get_config()

# Verify Notion
print(f"Notion API Key: {cfg.notion.api_key[:10]}...")
print(f"Dataset DB ID: {cfg.notion.dataset_db_id}")

# Verify Pinecone
print(f"Pinecone Index: {cfg.pinecone.index_name}")

# Verify OpenAI
print(f"OpenAI Chat Model: {cfg.openai.chat_model}")
```

### Test Notion Connection

```python
from amprenta_rag.clients.notion_client import notion_headers
import requests

headers = notion_headers()
url = f"https://api.notion.com/v1/databases/{cfg.notion.dataset_db_id}"
resp = requests.get(url, headers=headers)
print(f"Notion Status: {resp.status_code}")
```

### Test Pinecone Connection

```python
from amprenta_rag.clients.pinecone_client import get_pinecone_index

index = get_pinecone_index()
stats = index.describe_index_stats()
print(f"Pinecone Stats: {stats}")
```

---

## Troubleshooting

### Common Issues

1. **"Database ID not found"**
   - Verify the database ID is correct (32 hex characters with dashes)
   - Ensure the Notion integration has access to the database
   - Check database permissions in Notion

2. **"Pinecone index not found"**
   - Verify the index name matches exactly
   - Check Pinecone environment matches
   - Ensure API key has access to the index

3. **"OpenAI API error"**
   - Verify API key is valid
   - Check API key has sufficient credits
   - Verify model names are correct

4. **"Feature cache not working"**
   - Check `FEATURE_CACHE_ENABLED=true` in `.env`
   - Verify `FEATURE_CACHE_TTL` is set correctly
   - Check `FEATURE_CACHE_MAX_SIZE` is appropriate for your workload
   - Monitor cache stats: `python scripts/manage_feature_cache.py --stats`
   - See [Feature Caching Guide](FEATURE_CACHING.md) for detailed troubleshooting

---

## Security Best Practices

1. **Never commit `.env` file**
   - Add `.env` to `.gitignore`
   - Use environment variables in production
   - Rotate API keys regularly

2. **Limit Notion permissions**
   - Grant integration access only to required databases
   - Use read-only access where possible

3. **Monitor API usage**
   - Set up usage alerts for OpenAI
   - Monitor Pinecone index size
   - Track Notion API rate limits

---

## See Also

- [Feature Caching Guide](FEATURE_CACHING.md) - Performance optimization
- [API Reference](API_REFERENCE.md) - Module documentation
- [Architecture Overview](ARCHITECTURE.md) - System design
- [Usage Examples](USAGE_EXAMPLES.md) - Practical examples

