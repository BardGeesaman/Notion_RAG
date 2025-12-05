# Troubleshooting Guide

Common issues and solutions for the Amprenta RAG System.

## Table of Contents

1. [Installation Issues](#installation-issues)
2. [Configuration Problems](#configuration-problems)
3. [API Errors](#api-errors)
4. [Ingestion Problems](#ingestion-problems)
5. [Query Issues](#query-issues)
6. [Performance Problems](#performance-problems)
7. [Data Issues](#data-issues)

## Installation Issues

### ModuleNotFoundError

**Error**: `ModuleNotFoundError: No module named 'amprenta_rag'`

**Solution**:
```bash
# Make sure you're in the project root directory
cd /path/to/Notion_RAG

# Verify Python path includes project root
python -c "import sys; print(sys.path)"

# Install in development mode
pip install -e .
```

### RDKit Installation Issues

**Error**: `ModuleNotFoundError: No module named 'rdkit'`

**Solution**:
RDKit is best installed via conda:
```bash
conda install -c conda-forge rdkit
```

Or if using pip:
```bash
pip install rdkit-pypi
```

### Missing Dependencies

**Error**: Various import errors

**Solution**:
```bash
# Install all dependencies
pip install -r requirements.txt

# Verify installation
python scripts/validate_configuration.py
```

## Configuration Problems

### API Key Not Found

**Error**: `Configuration validation failed: Missing OPENAI_API_KEY`

**Solution**:
1. Check your `.env` file exists and contains the key:
   ```bash
   cat .env | grep OPENAI_API_KEY
   ```

2. Verify the key is set (no quotes needed):
   ```bash
   OPENAI_API_KEY=sk-...
   ```

3. Test configuration:
   ```bash
   python scripts/validate_configuration.py
   ```

### Database ID Errors

**Error**: `400 Bad Request` when accessing Notion databases

**Solution**:
1. Verify database IDs in `.env` are correct (no dashes):
   ```bash
   # Correct format
   NOTION_EXP_DATA_DB_ID=1234567890abcdef1234567890abcdef
   
   # Wrong format (has dashes)
   NOTION_EXP_DATA_DB_ID=12345678-90ab-cdef-1234-567890abcdef
   ```

2. Run database setup verification:
   ```bash
   python scripts/verify_notion_setup.py
   ```

3. See [Notion Database Setup Guide](NOTION_DATABASE_SETUP.md)

### Configuration Validation Fails

**Error**: `Configuration validation failed`

**Solution**:
```bash
# Run validation with detailed output
python scripts/validate_configuration.py --verbose

# Check specific components
python scripts/health_check.py
```

## API Errors

### Notion API Rate Limits

**Error**: `429 Too Many Requests`

**Solution**:
1. Reduce parallel workers:
   ```bash
   # Use fewer workers
   python scripts/batch_ingest_omics.py \
     --directory /path/to/data \
     --parallel \
     --max-workers 2  # Reduced from 5
   ```

2. Add delays between requests (handled automatically, but can be tuned)

3. Use caching to reduce API calls

### Notion API Authentication

**Error**: `401 Unauthorized`

**Solution**:
1. Verify API key is correct and active
2. Check integration has access to databases
3. Verify database is shared with integration (see [SHARE_DATABASES_WITH_INTEGRATION.md](SHARE_DATABASES_WITH_INTEGRATION.md))

### Pinecone Connection Issues

**Error**: `Connection timeout` or `Unable to connect to Pinecone`

**Solution**:
1. Verify API key:
   ```bash
   echo $PINECONE_API_KEY
   ```

2. Check index name in configuration:
   ```python
   # In amprenta_rag/config.py
   PINECONE_INDEX_NAME = "your-index-name"
   ```

3. Test connection:
   ```bash
   python -c "from amprenta_rag.clients.pinecone_client import get_pinecone_index; print(get_pinecone_index())"
   ```

### OpenAI API Errors

**Error**: `Rate limit exceeded` or `Insufficient quota`

**Solution**:
1. Check API quota in OpenAI dashboard
2. Reduce batch sizes for embeddings
3. Use cached embeddings when possible
4. Consider upgrading OpenAI plan

## Ingestion Problems

### File Not Found

**Error**: `FileNotFoundError: Lipidomics file not found`

**Solution**:
1. Check file path (use absolute paths if needed):
   ```bash
   python scripts/ingest_lipidomics.py \
     --file /absolute/path/to/file.csv \
     --create-page
   ```

2. Verify file exists:
   ```bash
   ls -lh /path/to/file.csv
   ```

### No Features Extracted

**Error**: `No species extracted from file` or `No features found`

**Solution**:
1. Check file format (should be CSV/TSV)
2. Verify column names match expected patterns
3. Check file encoding (should be UTF-8):
   ```bash
   file -i your_file.csv
   ```

4. Try manual inspection:
   ```python
   import pandas as pd
   df = pd.read_csv('your_file.csv')
   print(df.head())
   print(df.columns.tolist())
   ```

### Feature Normalization Fails

**Error**: Features not matching expected format

**Solution**:
1. Check normalization logs:
   ```bash
   python scripts/ingest_lipidomics.py --file your_file.csv --create-page 2>&1 | grep -i "normalize\|map"
   ```

2. Verify feature names are recognized
3. Check normalization function documentation

### Signature Ingestion Errors

**Error**: `Invalid signature format` or `Missing required columns`

**Solution**:
1. Verify TSV format has required columns:
   ```
   feature_type	feature_name	direction	weight
   ```

2. Check feature names are valid
3. Verify directions are valid: `↑`, `↓`, `neutral`, or empty

## Query Issues

### No Results Returned

**Error**: Query returns empty results

**Solution**:
1. Check if data is ingested:
   ```bash
   python scripts/rag_query.py \
     --query "test" \
     --source-type Dataset \
     --show-context
   ```

2. Verify filters aren't too restrictive:
   ```bash
   # Try without filters first
   python scripts/rag_query.py --query "your query" --top-k 10
   ```

3. Check Pinecone index has data:
   ```bash
   python -c "from amprenta_rag.clients.pinecone_client import get_pinecone_index; idx = get_pinecone_index(); print(idx.describe_index_stats())"
   ```

### Query Timeout

**Error**: Query takes too long or times out

**Solution**:
1. Reduce `--top-k` value:
   ```bash
   python scripts/rag_query.py \
     --query "your query" \
     --top-k 5  # Reduced from 10
   ```

2. Use more specific filters
3. Check system resources (CPU, memory)

### Cross-Omics Reasoning Errors

**Error**: `Error generating cross-omics summary`

**Solution**:
1. Verify linked datasets exist:
   ```bash
   python scripts/list_cross_omics_test_ids.py
   ```

2. Check database IDs are configured:
   ```bash
   python scripts/verify_notion_setup.py
   ```

3. Verify Program/Experiment pages exist and are linked

## Performance Problems

### Slow Ingestion

**Symptoms**: Ingestion takes very long

**Solution**:
1. Enable parallel processing:
   ```bash
   python scripts/batch_ingest_omics.py \
     --directory /path/to/data \
     --parallel \
     --max-workers 5
   ```

2. Use feature caching (enabled by default)
3. Check network speed for API calls
4. Verify system resources aren't constrained

### Memory Issues

**Error**: `MemoryError` or system becomes unresponsive

**Solution**:
1. Process files in smaller batches
2. Reduce parallel workers:
   ```bash
   --max-workers 2  # Instead of 5
   ```

3. Clear cache if needed:
   ```bash
   rm -rf .cache/feature_cache/*
   ```

### High API Costs

**Symptoms**: Unexpected API costs

**Solution**:
1. Use caching to avoid redundant API calls
2. Batch operations to reduce overhead
3. Review embedding usage (largest cost)
4. Consider using smaller embedding models for testing

## Data Issues

### Feature Linking Fails

**Error**: Features not linking to feature pages

**Solution**:
1. Verify feature database IDs are configured:
   ```bash
   python scripts/verify_notion_setup.py
   ```

2. Check feature pages exist in Notion
3. Review linking logs for errors

### Signature Scoring Issues

**Error**: Signature scores are 0 or unexpected

**Solution**:
1. Verify dataset has features:
   ```bash
   python scripts/find_datasets_with_features.py
   ```

2. Check feature names match signature components
3. Verify feature normalization is consistent
4. Review scoring logs:
   ```bash
   python scripts/score_signature.py \
     --dataset-id DATASET_ID \
     2>&1 | grep -i "score\|match"
   ```

### Metadata Inconsistencies

**Error**: Metadata doesn't match expected format

**Solution**:
1. Run metadata verification:
   ```bash
   python scripts/rag_verify_metadata.py
   ```

2. Check database schema matches expected format
3. Review metadata extraction logs

## Getting More Help

### Enable Debug Logging

```bash
# Set logging level
export LOG_LEVEL=DEBUG

# Or in Python
import logging
logging.basicConfig(level=logging.DEBUG)
```

### Check System Health

```bash
# Comprehensive health check
python scripts/health_check.py

# Configuration validation
python scripts/validate_configuration.py --verbose
```

### Review Logs

Check log files for detailed error messages:
```bash
# Find log files
find . -name "*.log" -type f

# View recent errors
grep -i error *.log | tail -20
```

### Common Solutions Summary

| Issue | Quick Fix |
|-------|-----------|
| ModuleNotFoundError | `pip install -r requirements.txt` |
| API Key Missing | Check `.env` file |
| Database Access | Run `verify_notion_setup.py` |
| No Query Results | Check data is ingested |
| Slow Performance | Enable parallel processing |
| Memory Issues | Reduce parallel workers |

## Reporting Issues

When reporting issues, include:
1. Error message (full traceback)
2. Command that triggered the error
3. Configuration (without sensitive keys)
4. System information (OS, Python version)
5. Relevant log output

Example:
```bash
# Capture full error
python scripts/your_script.py --args 2>&1 | tee error.log

# System info
python --version
uname -a
```

## Next Steps

- [Configuration Guide](CONFIGURATION.md) - Detailed configuration options
- [User Guide](USER_GUIDE.md) - Comprehensive usage guide
- [Architecture Overview](ARCHITECTURE.md) - System design details

