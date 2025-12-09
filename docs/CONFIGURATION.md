# Configuration Guide

**Last Updated**: 2025-12-03

This document describes how to configure the Amprenta RAG platform.

---

## Environment Variables

All configuration is managed via environment variables in a `.env` file in the project root.

### Required Configuration

#### Postgres Database

```bash
# Postgres connection settings
POSTGRES_HOST=localhost
POSTGRES_PORT=5432
POSTGRES_USER=your_username
POSTGRES_PASSWORD=your_password
POSTGRES_DB=amprenta

# Full connection URL (alternative to individual settings)
# DATABASE_URL=postgresql://user:password@localhost:5432/amprenta
```

**Setup**:
1. Install PostgreSQL (version 13+)
2. Create a database: `createdb amprenta`
3. Run migrations: `alembic upgrade head`

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
# Postgres Database (required)
POSTGRES_HOST=localhost
POSTGRES_PORT=5432
POSTGRES_USER=your_username
POSTGRES_PASSWORD=your_password
POSTGRES_DB=amprenta

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

## Database Schema

The Amprenta RAG platform uses PostgreSQL as its primary data store. The schema is managed via Alembic migrations.

### Core Tables

- **datasets**: Experimental data assets (omics datasets)
- **features**: Molecular features (genes, proteins, metabolites, lipids)
- **signatures**: Lipid signatures for pattern matching
- **programs**: Research programs grouping related experiments
- **experiments**: Individual experiments within programs
- **literature**: Literature references from Zotero

Run migrations to create the schema:
```bash
alembic upgrade head
```
- `Related Signature` (relation): Links to Lipid Signatures
- `Related Feature` (relation): Links to Gene/Protein/Metabolite/Lipid Species DB

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

# Verify Postgres
print(f"Postgres Host: {cfg.postgres.host}")
print(f"Postgres Database: {cfg.postgres.database}")

# Verify Pinecone
print(f"Pinecone Index: {cfg.pinecone.index_name}")

# Verify OpenAI
print(f"OpenAI Chat Model: {cfg.openai.chat_model}")
```

### Test Postgres Connection

```python
from amprenta_rag.database.base import get_db

db = next(get_db())
result = db.execute("SELECT 1")
print(f"Postgres Status: Connected")
db.close()
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

1. **"Database connection failed"**
   - Verify Postgres is running: `pg_isready`
   - Check database exists: `psql -l | grep amprenta`
   - Verify credentials in `.env`

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

2. **Secure Postgres access**
   - Use strong passwords
   - Limit network access to trusted hosts
   - Use read-only users where possible

3. **Monitor API usage**
   - Set up usage alerts for OpenAI
   - Monitor Pinecone index size
   - Track database query performance

---

## See Also

- [Feature Caching Guide](FEATURE_CACHING.md) - Performance optimization
- [API Reference](API_REFERENCE.md) - Module documentation
- [Architecture Overview](ARCHITECTURE.md) - System design
- [Usage Examples](USAGE_EXAMPLES.md) - Practical examples

