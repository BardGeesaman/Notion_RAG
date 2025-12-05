# How to Get API Keys for GEO and PRIDE Repositories

## Overview

Yes, authentication (API keys) could definitely be the issue! Both GEO and PRIDE repositories may require or benefit from API keys for reliable access.

---

## 🔑 GEO/NCBI API Key

### Why You Need It

- **Higher Rate Limits**: 10 requests/second vs 3 requests/second without key
- **More Reliable Access**: Better success rate for API calls
- **No IP Blocking**: Reduces risk of temporary IP blocks

### How to Get NCBI API Key

1. **Create NCBI Account** (if you don't have one)
   - Go to: https://www.ncbi.nlm.nih.gov/account/
   - Click "Register for an NCBI account"
   - Fill in your information and verify email

2. **Generate API Key**
   - Go to: https://www.ncbi.nlm.nih.gov/account/settings/
   - Scroll down to "API Key Management"
   - Click "Create an API Key"
   - Copy your API key (starts with something like: `abc123def456...`)

3. **Add to Configuration**
   - Open your `.env` file
   - Add:
     ```bash
     GEO_API_KEY=your_api_key_here
     ```
   - Save the file

### Testing API Key

```bash
# Test if API key works
python -c "
import os
os.environ['GEO_API_KEY'] = 'your_key_here'
from amprenta_rag.ingestion.repositories.geo import GEORepository
repo = GEORepository(api_key=os.environ['GEO_API_KEY'])
studies = repo.search_studies(['cancer'], max_results=5)
print(f'Found {len(studies)} studies')
"
```

---

## 🔑 PRIDE API Key (If Required)

### Current Status

PRIDE Archive API is currently publicly accessible and doesn't require authentication for public data. However, this may change in the future.

### If PRIDE Requires Authentication

1. **Check PRIDE Documentation**
   - Visit: https://www.ebi.ac.uk/pride/archive/api
   - Look for authentication requirements

2. **Register for API Access** (if needed)
   - Visit PRIDE Archive website
   - Create account
   - Generate API credentials

3. **Add to Configuration**
   ```bash
   PRIDE_API_KEY=your_pride_api_key_here
   PRIDE_API_SECRET=your_pride_secret_here  # If required
   ```

---

## ✅ Quick Setup Guide

### Step 1: Get NCBI API Key

```bash
# 1. Visit: https://www.ncbi.nlm.nih.gov/account/settings/
# 2. Create API key
# 3. Copy the key
```

### Step 2: Add to .env

```bash
# Add to your .env file
GEO_API_KEY=your_ncbi_api_key_here
```

### Step 3: Test

```bash
# Test GEO import with API key
python scripts/harvest_repository_study.py \
    --study-id GSE12345 \
    --repository GEO \
    --ingest
```

---

## 🧪 Verify API Keys Are Being Used

Check logs for API key usage:

```bash
# Look for rate limit info in logs
# With API key: faster rate limits (10 req/sec)
# Without API key: slower (3 req/sec)
```

---

## 📚 References

- **NCBI API Keys**: https://www.ncbi.nlm.nih.gov/account/settings/
- **NCBI API Documentation**: https://www.ncbi.nlm.nih.gov/books/NBK25497/
- **PRIDE API**: https://www.ebi.ac.uk/pride/archive/api

---

## 💡 Tips

1. **Free**: NCBI API keys are completely free
2. **Easy Setup**: Takes 5 minutes to get a key
3. **Highly Recommended**: Significantly improves GEO repository reliability
4. **Secure**: Store API keys in `.env` file (already in `.gitignore`)

---

## 🔧 Current Implementation

The system already supports API keys:

- ✅ GEO repository accepts API keys
- ✅ Configuration system ready (`GEO_API_KEY` in config)
- ✅ Discovery and harvest scripts will use API keys automatically once configured

Just add your API key to `.env` and it will be used automatically!

