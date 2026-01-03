# Secrets Management Overhaul

## Problem

FastAPI crashes on startup without `SIGNATURE_SECRET_KEY` and `JWT_SECRET_KEY`, even when `DISABLE_AUTH=true`. This blocks:
- Local development
- E2E tests (Playwright)
- CI/CD pipelines without secrets configured

## Solution Architecture

```
Priority Order (per Reviewer feedback):
1. AWS Secrets Manager (production)
2. Environment variables (fallback)
3. Dev defaults when DISABLE_AUTH=true (dev/test only - LAST RESORT)
```

## Implementation

### Phase 1: Fix Local Dev (P0) - 15 minutes

**File:** `amprenta_rag/utils/config_check.py`

Update `validate_required_secrets()` to skip validation when `DISABLE_AUTH=true`:

```python
def validate_required_secrets() -> Tuple[bool, List[str]]:
    """Validate all required secrets are present at startup.
    
    Returns:
        Tuple of (all_valid, list_of_missing_secrets)
    """
    # Skip validation in dev/test mode (NEVER use in production)
    if os.getenv("DISABLE_AUTH", "").lower() in ("1", "true"):
        logger.info("DISABLE_AUTH=true - skipping secret validation (DEV MODE ONLY)")
        return True, []
    
    # ... existing validation logic unchanged
```

### Phase 2: Unified Secrets Abstraction (P1) - 1-2 hours

**New File:** `amprenta_rag/config/secrets.py`

```python
import os
import logging
from functools import lru_cache
from typing import Optional

logger = logging.getLogger(__name__)

# Dev-only defaults (NEVER use in production)
_DEV_DEFAULTS = {
    "SIGNATURE_SECRET_KEY": "dev-signature-key-not-for-production",
    "JWT_SECRET_KEY": "dev-jwt-key-not-for-production",
}


def _get_from_aws(name: str) -> Optional[str]:
    """Attempt to get secret from AWS Secrets Manager."""
    try:
        import boto3
        from botocore.exceptions import ClientError
        
        client = boto3.client("secretsmanager")
        env = os.getenv("ENVIRONMENT", "dev")
        aws_name = f"amprenta/{env}/{name.lower()}"
        
        response = client.get_secret_value(SecretId=aws_name)
        return response.get("SecretString")
    except Exception as e:
        logger.debug(f"AWS Secrets Manager not available for {name}: {e}")
        return None


@lru_cache(maxsize=128)
def get_secret(name: str, default: Optional[str] = None) -> Optional[str]:
    """Get secret with priority: AWS -> env var -> default.
    
    In dev mode (DISABLE_AUTH=true), also falls back to dev defaults.
    
    Args:
        name: Secret name (e.g., "JWT_SECRET_KEY")
        default: Default value if not found
        
    Returns:
        Secret value or default
    """
    # Priority 1: AWS Secrets Manager (production)
    value = _get_from_aws(name)
    if value:
        return value
    
    # Priority 2: Environment variable
    value = os.getenv(name)
    if value:
        return value
    
    # Priority 3: Dev defaults (ONLY when DISABLE_AUTH=true)
    if os.getenv("DISABLE_AUTH", "").lower() in ("1", "true"):
        if name in _DEV_DEFAULTS:
            logger.warning(f"Using dev default for {name} - NOT FOR PRODUCTION")
            return _DEV_DEFAULTS[name]
    
    # Priority 4: Provided default
    return default
```

### Phase 3: AWS Secrets Manager Setup (P1) - 1 day

**Terraform:** `terraform/secrets.tf`

```hcl
resource "aws_secretsmanager_secret" "jwt_secret" {
  name = "amprenta/${var.environment}/jwt_secret_key"
}

resource "aws_secretsmanager_secret" "signature_secret" {
  name = "amprenta/${var.environment}/signature_secret_key"
}
```

**Migration:**
1. Create AWS secrets via Terraform
2. Update key `os.getenv()` calls to use `get_secret()`
3. Test with AWS credentials configured locally
4. Deploy to ECS (secrets injected via task definition)

## Files Changed

| File | Change |
|------|--------|
| `amprenta_rag/utils/config_check.py` | Skip validation when DISABLE_AUTH=true |
| `amprenta_rag/config/secrets.py` | New unified secrets abstraction |
| `terraform/secrets.tf` | AWS Secrets Manager resources |
| `docs/SECRETS_MANAGEMENT.md` | Update documentation |

## Security Notes

- `DISABLE_AUTH=true` is DEV/TEST ONLY
- Production deployments MUST use AWS Secrets Manager or env vars
- Dev defaults are clearly marked as not for production
- AWS secrets use environment-scoped naming: `amprenta/{env}/{key}`

## Testing

- Unit tests mock `get_secret()` - never hit real AWS
- E2E tests use `DISABLE_AUTH=true` (dev defaults)
- Integration tests can use real AWS with test secrets

## Batches

| Batch | Content | Effort |
|-------|---------|--------|
| 1 | Phase 1: config_check.py fix | 15 min |
| 2 | Phase 2: secrets.py abstraction | 1-2 hrs |
| 3 | Phase 3: Terraform + migration | 1 day |

