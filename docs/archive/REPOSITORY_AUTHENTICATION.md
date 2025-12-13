# Repository Authentication Guide

## Overview

Some repositories may require authentication (API keys) for reliable access. This guide explains authentication requirements and how to configure them.

---

## 🔑 Authentication Status

### ✅ No Authentication Required
- **Metabolomics Workbench (MW)** - Public API, no auth needed
- **MetaboLights** - Public API, no auth needed

### ⚠️ Optional Authentication (Recommended)
- **GEO (NCBI)** - Optional API key for higher rate limits and more reliable access
- **PRIDE** - May benefit from authentication for some endpoints

---

## 🔧 Configuration

### GEO/NCBI API Key

GEO uses NCBI Entrez E-utilities API. While authentication is not strictly required, an API key provides:
- Higher rate limits (10 requests/second vs 3 requests/second)
- More reliable access
- Better error handling

**How to Get NCBI API Key:**
1. Go to https://ncbiinsights.ncbi.nlm.nih.gov/2017/11/02/new-api-keys-for-the-e-utilities/
2. Or visit: https://www.ncbi.nlm.nih.gov/account/settings/
3. Create an account (if you don't have one)
4. Generate an API key
5. Add to your `.env` file:

```bash
GEO_API_KEY=your_ncbi_api_key_here
```

### PRIDE API Authentication

PRIDE Archive API is publicly accessible but may have rate limits. Currently, PRIDE doesn't require authentication for public data access, but some endpoints might work better with proper headers.

**If PRIDE API changes and requires authentication:**
1. Check PRIDE API documentation: https://www.ebi.ac.uk/pride/archive/api
2. Obtain API credentials if needed
3. Add to `.env`:

```bash
PRIDE_API_KEY=your_pride_api_key_here
PRIDE_API_SECRET=your_pride_api_secret_here  # If required
```

---

## 📝 Current Implementation

### GEO Repository

The GEO repository already supports optional API keys:

```python
# In amprenta_rag/config.py
GEO_API_KEY = os.getenv("GEO_API_KEY", "")  # Optional
```

The repository uses the API key when available:
- Higher rate limits (0.34s delay without key, faster with key)
- Included in all API requests automatically

### PRIDE Repository

PRIDE repository currently doesn't require authentication, but the implementation supports adding it if needed.

---

## 🚀 Enabling Authentication

### Step 1: Get API Keys

**For GEO/NCBI:**
1. Visit https://www.ncbi.nlm.nih.gov/account/settings/
2. Create account or log in
3. Generate API key
4. Copy the key

**For PRIDE:**
- Currently not required
- Check PRIDE documentation if you encounter access issues

### Step 2: Add to Configuration

Add to your `.env` file:

```bash
# GEO/NCBI API Key (optional but recommended)
GEO_API_KEY=your_ncbi_api_key_here

# PRIDE API Key (if required in future)
PRIDE_API_KEY=your_pride_api_key_here
```

### Step 3: Restart Application

Restart your application/scripts to load the new API keys.

---

## 🧪 Testing Authentication

Test if API keys help:

```bash
# Test GEO with API key
python -c "
from amprenta_rag.ingestion.repositories.geo import GEORepository
from amprenta_rag.config import get_config

cfg = get_config()
repo = GEORepository(api_key=cfg.geo_api_key if hasattr(cfg, 'geo_api_key') else None)
studies = repo.search_studies(['cancer'], max_results=5)
print(f'Found {len(studies)} studies')
"
```

---

## ⚠️ Troubleshooting

### GEO Still Not Working

1. **Check API Key**: Verify `GEO_API_KEY` is set correctly in `.env`
2. **Rate Limits**: NCBI has strict rate limits, ensure delays are working
3. **Search Terms**: Try different search queries - GEO database may need specific formats
4. **Direct Access**: Try accessing studies directly by ID instead of searching

### PRIDE Still Not Working

1. **Check API Status**: PRIDE API might be temporarily down
2. **Endpoint Changes**: PRIDE may have changed API endpoints
3. **Project IDs**: Verify project IDs are valid and publicly accessible
4. **Headers**: Ensure proper Accept headers are sent (already implemented)

---

## 📚 References

- **NCBI API Documentation**: https://www.ncbi.nlm.nih.gov/books/NBK25497/
- **NCBI API Keys**: https://www.ncbi.nlm.nih.gov/account/settings/
- **PRIDE API Documentation**: https://www.ebi.ac.uk/pride/archive/api
- **Entrez E-utilities**: https://eutils.ncbi.nlm.nih.gov/

---

## 💡 Recommendations

1. **For GEO**: Get an NCBI API key - it's free and significantly improves reliability
2. **For PRIDE**: Monitor API documentation for authentication changes
3. **For Production**: Always use API keys when available for better reliability
4. **Rate Limiting**: Respect rate limits to avoid IP blocking

---

## 🔄 Current Status

- ✅ **GEO API Key Support**: Already implemented, just needs to be configured
- ✅ **PRIDE Authentication**: Framework ready, not currently required
- ✅ **Configuration System**: Ready to accept API keys
- ⚠️ **Working Without Keys**: Both repositories should work, but keys improve reliability

