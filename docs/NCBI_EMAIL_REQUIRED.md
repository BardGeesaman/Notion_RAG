# NCBI Email Configuration Required

## ⚠️ IMPORTANT: NCBI Requires Email Address

NCBI requires an email address for API access. This is a mandatory requirement by NCBI policy.

## 🔧 How to Configure

### Step 1: Get Your Email Address

Use your email address (any valid email is fine).

### Step 2: Add to .env File

Add this line to your `.env` file:

```bash
NCBI_EMAIL=your.email@example.com
```

### Example

```bash
# .env file
GEO_API_KEY=b625510afa75c7df47141f6b82baf70d4008
NCBI_EMAIL=your.email@example.com
```

## 📝 Why It's Required

NCBI uses email addresses to:
- Identify API users
- Contact users if there are issues
- Monitor API usage
- Follow their usage policies

This is a **mandatory requirement**, not optional.

## ✅ After Configuration

Once configured, the GEO repository will automatically use:
- Your email address (from NCBI_EMAIL)
- Your API key (from GEO_API_KEY)

Both are set globally in Bio.Entrez before any API calls.

