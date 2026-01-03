# Secrets Management

## Overview
Unified secrets abstraction with priority order:
1. AWS Secrets Manager (production)
2. Environment variables (fallback)
3. Dev defaults when DISABLE_AUTH=true (dev/test only)

## Usage

```python
from amprenta_rag.config_secrets import get_secret

jwt_key = get_secret("JWT_SECRET_KEY")
sig_key = get_secret("SIGNATURE_SECRET_KEY")
custom = get_secret("CUSTOM_KEY", default="fallback")
```

## Development Mode

Set `DISABLE_AUTH=true` to use dev defaults:
```bash
DISABLE_AUTH=true uvicorn amprenta_rag.api.main:app --port 8000
```

⚠️ **Security Warning**: Dev defaults are NOT FOR PRODUCTION.
Production deployments will raise RuntimeError if DISABLE_AUTH is set.

## AWS Secrets Manager (Future)

Naming convention: `amprenta/{env}/{key}`
- `amprenta/dev/jwt_secret_key`
- `amprenta/prod/signature_secret_key`

Phase 3 (Terraform integration) is deferred.

## Testing

Mock `get_secret()` in tests - never hit real AWS:
```python
from unittest.mock import patch

@patch("amprenta_rag.config_secrets.get_secret")
def test_example(mock_get_secret):
    mock_get_secret.return_value = "test-key"
    # ... test code
```