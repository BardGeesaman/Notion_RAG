# Secrets Management Guide

## Overview

Amprenta RAG uses a dual-mode secrets management system that automatically adapts to your environment:

- **Local Development**: Environment variables (`.env` file)
- **Production (AWS)**: AWS Secrets Manager

The application automatically detects the environment and retrieves secrets accordingly, providing a seamless experience for developers and secure, centralized management in production.

## Quick Start (Local Development)

### 1. Copy the example file
```bash
cp .env.example .env
```

### 2. Fill in required values (minimum for basic functionality)
```bash
# Edit .env file
POSTGRES_PASSWORD=your_secure_password
OPENAI_API_KEY=sk-...  # For LLM features
NCBI_EMAIL=your.email@example.com  # For GEO queries
```

### 3. Start the application
```bash
# Activate conda environment
conda activate myenv

# Start database (if using Docker)
docker-compose up -d postgres

# Run database migrations
alembic upgrade head

# Start API server
python -m uvicorn amprenta_rag.api.main:app --reload

# Start dashboard (in another terminal)
streamlit run scripts/dashboard/main.py
```

### 4. Verify setup
```bash
# Test database connection
python -c "from amprenta_rag.utils.secrets import get_database_url; print('Database:', get_database_url()[:50] + '...')"

# Test API keys
python -c "from amprenta_rag.utils.secrets import get_api_key; print('OpenAI configured:', bool(get_api_key('openai')))"
```

## Secret Categories

The secrets are organized into logical groups for better management:

| Category | Secrets | Required For | AWS Secret Group |
|----------|---------|--------------|------------------|
| **Database** | `POSTGRES_PASSWORD`, `POSTGRES_URL` | Core functionality | `database` |
| **API Keys** | `OPENAI_API_KEY`, `ZOTERO_API_KEY`, `NOTION_API_KEY`, `GEO_API_KEY` | External integrations | `api-keys` |
| **Auth** | `NCBI_EMAIL` | NCBI/GEO queries | `auth` |
| **Integrations** | 13× `NOTION_*_DB_ID` | Notion sync (legacy) | `integrations` |
| **Backup** | `BACKUP_KMS_KEY_ID` | Encrypted backups | `backup` |

### Required vs Optional Secrets

**Absolutely Required** (for core functionality):
- `POSTGRES_PASSWORD` - Database access
- `OPENAI_API_KEY` - LLM features (RAG, chat, summarization)

**Recommended** (for full functionality):
- `NCBI_EMAIL` - Required by NCBI for API access
- `ZOTERO_API_KEY` - Bibliography management

**Optional** (for specific features):
- `NOTION_API_KEY` + database IDs - Legacy Notion sync
- `GEO_API_KEY` - Higher rate limits for GEO queries
- `BACKUP_KMS_KEY_ID` - Encrypted backups

## Local Development Setup

### Using the Interactive Setup Helper

For guided setup, use the interactive helper:

```bash
python scripts/setup_env.py
```

This will:
1. Check for existing `.env` file
2. Copy `.env.example` to `.env`
3. Prompt for required values
4. Provide next steps

### Manual Setup

1. **Copy environment template**:
   ```bash
   cp .env.example .env
   ```

2. **Set database credentials**:
   ```bash
   # Required for local database access
   POSTGRES_PASSWORD=your_secure_password
   
   # Optional: Override other database settings
   POSTGRES_USER=postgres
   POSTGRES_HOST=localhost
   POSTGRES_PORT=5432
   POSTGRES_DB=amprenta
   ```

3. **Configure API keys** (as needed):
   ```bash
   # OpenAI for LLM features
   OPENAI_API_KEY=sk-proj-...
   
   # Zotero for bibliography
   ZOTERO_API_KEY=your_zotero_key
   ZOTERO_LIBRARY_ID=your_library_id
   
   # NCBI for public data
   NCBI_EMAIL=your.email@example.com
   GEO_API_KEY=your_geo_key  # Optional
   ```

4. **Set feature flags** (optional):
   ```bash
   # Development settings
   DEBUG=true
   LOG_LEVEL=DEBUG
   DISABLE_AUTH=true  # Skip authentication for local dev
   
   # Environment designation
   ENVIRONMENT=dev
   ```

### Environment Detection

The application automatically detects where it's running:

```python
from amprenta_rag.utils.secrets import get_environment_info
print(get_environment_info())
```

**Local Environment** (IS_AWS=False):
- Detected when: No AWS environment variables present
- Secrets source: Environment variables (`.env` file)
- Fallback: Direct `os.getenv()` calls

**AWS Environment** (IS_AWS=True):
- Detected when: `AWS_EXECUTION_ENV` or `ECS_CONTAINER_METADATA_URI` present
- Secrets source: AWS Secrets Manager via boto3
- Fallback: Environment variables if AWS unavailable

## AWS Secrets Manager (Production)

### Secret Groups Architecture

In AWS, secrets are organized into 5 groups under `amprenta/{environment}/`:

```
amprenta/
├── prod/
│   ├── database          # PostgreSQL credentials
│   ├── api-keys          # External service API keys
│   ├── auth              # Authentication credentials
│   ├── integrations      # Third-party integration IDs
│   └── backup            # Backup encryption keys
├── staging/
│   └── ... (same structure)
└── dev/
    └── ... (same structure)
```

### Secret Group Contents

#### 1. Database (`amprenta/{env}/database`)
```json
{
  "POSTGRES_PASSWORD": "secure_generated_password",
  "POSTGRES_URL": "postgresql://user:password@host:5432/db"
}
```

#### 2. API Keys (`amprenta/{env}/api-keys`)
```json
{
  "OPENAI_API_KEY": "sk-proj-...",
  "ZOTERO_API_KEY": "your_zotero_key",
  "NOTION_API_KEY": "secret_...",
  "GEO_API_KEY": "your_geo_key"
}
```

#### 3. Auth (`amprenta/{env}/auth`)
```json
{
  "NCBI_EMAIL": "your.email@example.com"
}
```

#### 4. Integrations (`amprenta/{env}/integrations`)
```json
{
  "NOTION_EXP_DATA_DB_ID": "2b7adf6142ab...",
  "NOTION_METABOLITE_FEATURES_DB_ID": "3c8bef7253bc...",
  "NOTION_PROTEIN_FEATURES_DB_ID": "4d9cfg8364cd...",
  "NOTION_GENE_FEATURES_DB_ID": "5e0dhg9475de...",
  "NOTION_SIGNATURE_DB_ID": "6f1ehi0586ef...",
  "NOTION_SIGNATURE_COMPONENT_DB_ID": "7g2fij1697fg...",
  "NOTION_LIPID_SPECIES_DB_ID": "8h3gjk2708gh...",
  "NOTION_PROGRAMS_DB_ID": "9i4hkl3819hi...",
  "NOTION_EXPERIMENTS_DB_ID": "0j5ilm4920ij...",
  "NOTION_COMPOUND_FEATURES_DB_ID": "1k6jmn5031jk...",
  "NOTION_HTS_CAMPAIGNS_DB_ID": "2l7kno6142kl...",
  "NOTION_BIOCHEMICAL_HITS_DB_ID": "3m8lop7253lm...",
  "NOTION_PATHWAYS_DB_ID": "4n9mpq8364mn..."
}
```

#### 5. Backup (`amprenta/{env}/backup`)
```json
{
  "BACKUP_KMS_KEY_ID": "arn:aws:kms:us-east-1:123456789:key/..."
}
```

### ECS Integration

ECS tasks automatically inject secrets as environment variables through the task definition:

```hcl
secrets = [
  { name = "POSTGRES_URL", valueFrom = "${aws_secretsmanager_secret.database.arn}:POSTGRES_URL::" },
  { name = "OPENAI_API_KEY", valueFrom = "${aws_secretsmanager_secret.api_keys.arn}:OPENAI_API_KEY::" },
  # ... all 21 secrets
]
```

**Benefits**:
- No code changes between local and production
- Secrets never appear in logs or process lists
- Automatic rotation support
- Fine-grained IAM permissions

### Populating Secrets in AWS

After running `terraform apply`, manually populate secrets in the AWS Console:

1. **Navigate to AWS Secrets Manager**
2. **Find secrets** under `amprenta/{environment}/`
3. **Edit each secret** and add key-value pairs as JSON
4. **Verify** ECS tasks can access secrets (check IAM permissions)

#### Example: Populating Database Secret

1. Go to AWS Secrets Manager → `amprenta/prod/database`
2. Click "Retrieve secret value" → "Edit"
3. Replace with:
   ```json
   {
     "POSTGRES_PASSWORD": "your_secure_production_password",
     "POSTGRES_URL": "postgresql://amprenta_user:your_secure_production_password@your-rds-endpoint:5432/amprenta"
   }
   ```
4. Save changes

## Security Best Practices

### Development Environment

1. **Never commit `.env` files**
   - Already in `.gitignore`
   - Contains sensitive credentials

2. **Use strong passwords**
   - Minimum 16 characters for database passwords
   - Mix of letters, numbers, symbols

3. **Rotate API keys regularly**
   - Especially if potentially compromised
   - Update both local `.env` and AWS secrets

4. **Limit local permissions**
   - Don't share `.env` files
   - Use separate API keys for development vs production

### Production Environment

1. **Use AWS IAM roles**
   - No hardcoded AWS credentials
   - Principle of least privilege

2. **Enable AWS CloudTrail**
   - Audit secret access
   - Monitor for unusual patterns

3. **Regular secret rotation**
   - Quarterly for high-value secrets
   - Automated where possible

4. **Network security**
   - VPC-only access to secrets
   - No internet gateway for production subnets

## Troubleshooting

### "Secret not found" errors

**Check environment detection**:
```python
from amprenta_rag.utils.secrets import get_environment_info
print(get_environment_info())
```

**Expected output (local)**:
```json
{
  "is_aws": false,
  "environment": "dev",
  "aws_region": null,
  "execution_env": null,
  "ecs_metadata_uri": null,
  "cache_size": 0
}
```

**Expected output (AWS)**:
```json
{
  "is_aws": true,
  "environment": "prod",
  "aws_region": "us-east-1",
  "execution_env": "AWS_ECS_FARGATE",
  "ecs_metadata_uri": "http://169.254.170.2/v4/...",
  "cache_size": 3
}
```

### Local environment issues

**Verify `.env` file is loaded**:
```bash
# Should show your password
echo $POSTGRES_PASSWORD

# If empty, check file location and syntax
cat .env | grep POSTGRES_PASSWORD
```

**Test database connection**:
```python
from amprenta_rag.utils.secrets import get_database_url
print("Database URL:", get_database_url())

# Test actual connection
from amprenta_rag.database.session import db_session
with db_session() as db:
    result = db.execute("SELECT 1").scalar()
    print("Database connection:", "OK" if result == 1 else "FAILED")
```

**Check dotenv loading**:
```python
import os
from dotenv import load_dotenv

# Manual load
load_dotenv()
print("POSTGRES_PASSWORD loaded:", bool(os.getenv("POSTGRES_PASSWORD")))
```

### AWS permission errors

**Common IAM permission issues**:

1. **ECS task role missing permissions**:
   ```json
   {
     "Effect": "Allow",
     "Action": [
       "secretsmanager:GetSecretValue",
       "secretsmanager:DescribeSecret"
     ],
     "Resource": [
       "arn:aws:secretsmanager:region:account:secret:amprenta/prod/*"
     ]
   }
   ```

2. **Wrong secret ARN in task definition**:
   - Check secret exists in AWS console
   - Verify ARN format: `arn:aws:secretsmanager:region:account:secret:name`

3. **Cross-region access**:
   - Secrets and ECS tasks must be in same region
   - Check `AWS_REGION` environment variable

**Debug AWS secrets access**:
```python
# Test boto3 access
import boto3
from botocore.exceptions import ClientError

try:
    client = boto3.client("secretsmanager")
    response = client.get_secret_value(SecretId="amprenta/prod/database")
    print("AWS access: OK")
except ClientError as e:
    print("AWS error:", e.response['Error']['Code'])
except Exception as e:
    print("Other error:", e)
```

### Application startup issues

**Missing required secrets**:
```python
# Check all required secrets
from amprenta_rag.utils.secrets import get_database_url, get_api_key

required_checks = [
    ("Database URL", get_database_url()),
    ("OpenAI API Key", get_api_key("openai")),
]

for name, value in required_checks:
    status = "✅ OK" if value else "❌ MISSING"
    print(f"{name}: {status}")
```

**Config loading errors**:
```python
try:
    from amprenta_rag.config import get_config
    config = get_config()
    print("✅ Config loaded successfully")
    
    # Check key components
    print(f"Database configured: {bool(config.postgres.url)}")
    print(f"OpenAI configured: {bool(config.openai.api_key)}")
    
except Exception as e:
    print("❌ Config loading failed:", e)
    import traceback
    traceback.print_exc()
```

## Migration Guide

### From hardcoded values to secrets

1. **Identify the variable** in `docs/SECRETS_INVENTORY.md`
2. **Add to AWS Secrets Manager** (production)
3. **Add to `.env.example` and your local `.env`** (development)
4. **Update code** to use `get_secret()` or appropriate helper
5. **Test** both local and AWS environments

### From Terraform variables to Secrets Manager

**Before** (insecure):
```hcl
variable "db_password" {
  description = "Database password"
  type        = string
  sensitive   = true
}

resource "aws_db_instance" "postgres" {
  password = var.db_password  # Visible in state!
}
```

**After** (secure):
```hcl
# No password variable needed
resource "aws_db_instance" "postgres" {
  manage_master_user_password = true  # AWS-generated
}

# Application gets password from Secrets Manager
resource "aws_secretsmanager_secret" "database" {
  name = "amprenta/${var.environment}/database"
}
```

**Migration steps**:
1. Update Terraform to remove password variables
2. `terraform apply` creates AWS-managed password
3. Retrieve password from RDS console or Secrets Manager
4. Update application secrets with the new password
5. Deploy updated application code

## API Reference

### Core Functions

```python
from amprenta_rag.utils.secrets import (
    get_secret,
    get_database_url,
    get_api_key,
    get_auth_credential,
    get_integration_id,
    get_backup_config,
    clear_cache,
    get_environment_info
)
```

#### `get_secret(secret_name, key=None)`

Get secret from AWS Secrets Manager or environment variable.

**Parameters**:
- `secret_name` (str): AWS secret name (e.g., "amprenta/prod/database")
- `key` (str, optional): JSON key within the secret

**Returns**: Secret value or None if not found

**Examples**:
```python
# Get entire secret (if it's a simple string)
password = get_secret("amprenta/prod/simple-password")

# Get specific key from JSON secret
db_password = get_secret("amprenta/prod/database", "POSTGRES_PASSWORD")
```

#### `get_database_url()`

Get database connection URL.

**Returns**: PostgreSQL connection URL

**Local**: Constructs from individual `POSTGRES_*` environment variables
**AWS**: Uses `POSTGRES_URL` from database secret group

#### `get_api_key(service)`

Get API key for external service.

**Parameters**:
- `service` (str): Service name ("openai", "zotero", "notion", "geo")

**Returns**: API key or None if not found

#### `get_auth_credential(credential)`

Get authentication credential.

**Parameters**:
- `credential` (str): Credential name ("ncbi_email")

**Returns**: Credential value or None if not found

#### `get_integration_id(integration)`

Get integration identifier (e.g., Notion database IDs).

**Parameters**:
- `integration` (str): Integration name ("experiments", "compounds", etc.)

**Returns**: Integration ID or None if not found

#### `get_backup_config(config)`

Get backup configuration value.

**Parameters**:
- `config` (str): Configuration name ("kms_key_id")

**Returns**: Configuration value or None if not found

#### `clear_cache()`

Clear the secrets cache. Useful for testing or when secrets are rotated.

#### `get_environment_info()`

Get information about the current environment and secrets configuration.

**Returns**: Dictionary with environment information

### Error Handling

All functions handle errors gracefully:
- Missing secrets return `None` (not exceptions)
- AWS errors are logged and fall back to environment variables
- Network issues are retried with exponential backoff
- Invalid JSON in secrets is logged as error

### Caching

AWS Secrets Manager calls are cached in memory to avoid repeated API calls:
- Cache is per-secret, per-key
- Cache persists for the lifetime of the process
- Use `clear_cache()` to force refresh
- Cache size is included in `get_environment_info()`

## Best Practices Summary

### For Developers

1. **Use the setup helper**: `python scripts/setup_env.py`
2. **Start minimal**: Only set required secrets initially
3. **Test incrementally**: Verify each secret as you add it
4. **Keep `.env` private**: Never commit or share

### For DevOps

1. **Use AWS IAM roles**: No hardcoded credentials
2. **Organize by environment**: Separate secrets for dev/staging/prod
3. **Monitor access**: Enable CloudTrail logging
4. **Rotate regularly**: Quarterly for high-value secrets

### For Security

1. **Principle of least privilege**: Minimal required permissions
2. **Audit regularly**: Review IAM policies and access logs
3. **Encrypt in transit**: Use TLS for all secret transfers
4. **Backup secrets**: Secure backup of critical secrets

---

## Additional Resources

- [AWS Secrets Manager Documentation](https://docs.aws.amazon.com/secretsmanager/)
- [Terraform AWS Provider - Secrets Manager](https://registry.terraform.io/providers/hashicorp/aws/latest/docs/resources/secretsmanager_secret)
- [ECS Secrets from Parameter Store and Secrets Manager](https://docs.aws.amazon.com/AmazonECS/latest/developerguide/specifying-sensitive-data-secrets.html)
- [Project Secrets Inventory](./SECRETS_INVENTORY.md)
