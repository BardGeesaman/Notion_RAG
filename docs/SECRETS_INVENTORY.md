# Secrets Inventory

**Generated**: 2025-01-01  
**Source**: `amprenta_rag/config.py`  
**Total Environment Variables**: 57

## Summary

- **Total environment variables**: 57
- **Required secrets**: 19
- **Optional secrets**: 6
- **Required configuration**: 7  
- **Optional configuration**: 24

## AWS Secrets Manager Groups

### Group: `amprenta/{env}/database`

| Variable | Category | Default | Description |
|----------|----------|---------|-------------|
| POSTGRES_PASSWORD | REQUIRED_SECRET | - | Database password |
| POSTGRES_URL | OPTIONAL_SECRET | - | Full database connection string |

### Group: `amprenta/{env}/api-keys`

| Variable | Category | Default | Description |
|----------|----------|---------|-------------|
| OPENAI_API_KEY | REQUIRED_SECRET | - | OpenAI API key for LLM operations |
| ZOTERO_API_KEY | REQUIRED_SECRET | - | Zotero API key for literature management |
| NOTION_API_KEY | OPTIONAL_SECRET | - | Notion API key (only if ENABLE_NOTION_SYNC=true) |
| GEO_API_KEY | OPTIONAL_SECRET | - | NCBI GEO API key for higher rate limits |

### Group: `amprenta/{env}/auth`

| Variable | Category | Default | Description |
|----------|----------|---------|-------------|
| NCBI_EMAIL | REQUIRED_SECRET | - | Email required by NCBI for API access |
| SIGNATURE_SECRET_KEY | REQUIRED_SECRET | - | HMAC key for electronic signatures (32+ hex chars) |
| JWT_SECRET_KEY | REQUIRED_SECRET | - | JWT signing key for authentication (32+ hex chars) |

### Group: `amprenta/{env}/integrations`

| Variable | Category | Default | Description |
|----------|----------|---------|-------------|
| NOTION_EXP_DATA_DB_ID | REQUIRED_SECRET | - | Notion Experimental Data Assets DB ID |
| NOTION_METABOLITE_FEATURES_DB_ID | REQUIRED_SECRET | - | Notion Metabolite Features DB ID |
| NOTION_PROTEIN_FEATURES_DB_ID | REQUIRED_SECRET | - | Notion Protein Features DB ID |
| NOTION_GENE_FEATURES_DB_ID | REQUIRED_SECRET | - | Notion Gene Features DB ID |
| NOTION_SIGNATURE_DB_ID | REQUIRED_SECRET | - | Notion Lipid Signatures DB ID |
| NOTION_SIGNATURE_COMPONENT_DB_ID | REQUIRED_SECRET | - | Notion Signature Components DB ID |
| NOTION_LIPID_SPECIES_DB_ID | REQUIRED_SECRET | - | Notion Lipid Species DB ID |
| NOTION_PROGRAMS_DB_ID | REQUIRED_SECRET | - | Notion Programs DB ID |
| NOTION_EXPERIMENTS_DB_ID | REQUIRED_SECRET | - | Notion Experiments DB ID |
| NOTION_COMPOUND_FEATURES_DB_ID | REQUIRED_SECRET | - | Notion Compound Features DB ID |
| NOTION_HTS_CAMPAIGNS_DB_ID | REQUIRED_SECRET | - | Notion HTS Campaigns DB ID |
| NOTION_BIOCHEMICAL_HITS_DB_ID | REQUIRED_SECRET | - | Notion Biochemical Hits DB ID |
| NOTION_PATHWAYS_DB_ID | REQUIRED_SECRET | - | Notion Pathways DB ID |

### Group: `amprenta/{env}/backup`

| Variable | Category | Default | Description |
|----------|----------|---------|-------------|
| BACKUP_KMS_KEY_ID | OPTIONAL_SECRET | - | AWS KMS key ID for backup encryption |

### Group: `amprenta/{env}/email`

| Variable | Category | Default | Description |
|----------|----------|---------|-------------|
| SMTP_HOST | OPTIONAL_CONFIG | smtp.gmail.com | SMTP server hostname |
| SMTP_PORT | OPTIONAL_CONFIG | 587 | SMTP server port |
| SMTP_USER | OPTIONAL_SECRET | - | SMTP username/email |
| SMTP_PASSWORD | OPTIONAL_SECRET | - | SMTP password or app password |
| FROM_EMAIL | OPTIONAL_CONFIG | - | Override sender email address |

## Environment Variables (Non-Secret)

| Variable | Category | Default | Description |
|----------|----------|---------|-------------|
| POSTGRES_HOST | REQUIRED_CONFIG | localhost | Database host |
| POSTGRES_PORT | REQUIRED_CONFIG | 5432 | Database port |
| POSTGRES_DB | REQUIRED_CONFIG | amprenta | Database name |
| POSTGRES_USER | REQUIRED_CONFIG | bard | Database username |
| POSTGRES_ECHO | OPTIONAL_CONFIG | false | Enable SQLAlchemy query logging |
| VECTOR_BACKEND | REQUIRED_CONFIG | pgvector | Vector database backend |
| CORS_ORIGINS | REQUIRED_CONFIG | http://localhost:8501 | CORS allowed origins |
| SIGNATURES_DIR | REQUIRED_CONFIG | - | Directory for signature files |
| ENABLE_NOTION_SYNC | OPTIONAL_CONFIG | false | Enable Notion synchronization |
| ENABLE_DUAL_WRITE | OPTIONAL_CONFIG | false | Enable dual-write mode |
| USE_POSTGRES_AS_SOT | OPTIONAL_CONFIG | true | Use Postgres as source of truth |
| ENABLE_SIGNATURE_SCORING | OPTIONAL_CONFIG | true | Enable signature scoring |
| ENABLE_LIPID_MAPPING | OPTIONAL_CONFIG | true | Enable lipid mapping |
| ENABLE_FEATURE_LINKING | OPTIONAL_CONFIG | true | Enable feature linking |
| ENABLE_LLM_SEMANTIC_EXTRACTION | OPTIONAL_CONFIG | false | Enable LLM semantic extraction |
| FEATURE_LINKING_MAX_WORKERS | OPTIONAL_CONFIG | 10 | Max workers for feature linking |
| FEATURE_CACHE_ENABLED | OPTIONAL_CONFIG | true | Enable feature caching |
| FEATURE_CACHE_TTL | OPTIONAL_CONFIG | 3600 | Feature cache TTL in seconds |
| FEATURE_CACHE_TTL_SECONDS | OPTIONAL_CONFIG | 3600 | Legacy feature cache TTL |
| FEATURE_CACHE_MAX_SIZE | OPTIONAL_CONFIG | 1000 | Max feature cache size |
| FEATURE_CACHE_ENABLE_PERSISTENCE | OPTIONAL_CONFIG | true | Enable persistent feature cache |
| FEATURE_CACHE_DIR | OPTIONAL_CONFIG | - | Feature cache directory |
| FEATURE_CACHE_PARALLEL_WORKERS | OPTIONAL_CONFIG | 5 | Parallel workers for cache |
| SIGNATURE_OVERLAP_THRESHOLD | OPTIONAL_CONFIG | 0.3 | Signature overlap threshold |
| AUTO_LINK_ENABLED | OPTIONAL_CONFIG | true | Enable auto-linking |
| AUTO_LINK_MIN_CONFIDENCE | OPTIONAL_CONFIG | 0.8 | Auto-link minimum confidence |
| BACKUP_S3_ENABLED | OPTIONAL_CONFIG | false | Enable S3 backup |
| BACKUP_S3_BUCKET | OPTIONAL_CONFIG | - | S3 bucket for backups |
| BACKUP_LOCAL_DIR | OPTIONAL_CONFIG | ./backups | Local backup directory |
| BACKUP_RETENTION_DAYS | OPTIONAL_CONFIG | 365 | Backup retention period |

## Decision Matrix

### AWS Secrets Manager (Secrets)

**Required Secrets (11 variables)**:
- `OPENAI_API_KEY` - Critical for LLM operations
- `ZOTERO_API_KEY` - Required for literature management
- `POSTGRES_PASSWORD` - Database access credential
- `NCBI_EMAIL` - Required by NCBI API terms
- All Notion DB IDs (10 variables) - Sensitive integration identifiers

**Optional Secrets (3 variables)**:
- `NOTION_API_KEY` - Only needed if Notion sync enabled
- `GEO_API_KEY` - Optional performance enhancement
- `POSTGRES_URL` - Alternative to individual DB params
- `BACKUP_KMS_KEY_ID` - Optional backup encryption

### ECS Environment Variables (Non-Secrets)

**Required Configuration (8 variables)**:
- Database connection parameters (host, port, db, user)
- Service configuration (vector backend, CORS origins)
- Required directory paths

**Optional Configuration (30 variables)**:
- Feature flags and toggles
- Performance tuning parameters
- Cache configuration
- Backup settings

### Application Constants (Hardcoded)

**Should be hardcoded in application**:
- Default model names (`OPENAI_CHAT_MODEL`, `OPENAI_EMBED_MODEL`)
- API versions (`NOTION_VERSION`)
- Static IDs that don't change per environment
- Library identifiers (`ZOTERO_LIBRARY_ID`)

## Security Considerations

### High Priority for AWS Secrets Manager

1. **API Keys**: All external service API keys must be in Secrets Manager
2. **Database Credentials**: Passwords and connection strings
3. **Integration IDs**: Notion database IDs are sensitive and environment-specific

### Medium Priority

1. **Optional API Keys**: GEO API key can remain in ECS environment variables
2. **Email Addresses**: NCBI email could be ECS env var but treating as secret for compliance

### Low Priority (Keep in ECS)

1. **Feature Flags**: Boolean toggles for application behavior
2. **Performance Parameters**: Cache sizes, worker counts, thresholds
3. **Directory Paths**: File system paths and backup locations

## Migration Strategy

### Phase 1: Critical Secrets
- Move all REQUIRED_SECRET variables to AWS Secrets Manager
- Update application to read from Secrets Manager

### Phase 2: Optional Secrets  
- Move OPTIONAL_SECRET variables to AWS Secrets Manager
- Maintain fallback to environment variables

### Phase 3: Configuration Review
- Review all OPTIONAL_CONFIG variables
- Consider hardcoding stable values
- Consolidate related parameters

## Ambiguous Items Requiring Decision

### 1. NCBI_EMAIL
- **Current**: Treated as REQUIRED_SECRET
- **Consideration**: Could be REQUIRED_CONFIG (not sensitive)
- **Recommendation**: Keep as secret for compliance/audit trail

### 2. Notion Database IDs
- **Current**: Treated as REQUIRED_SECRET
- **Consideration**: Could be REQUIRED_CONFIG (environment-specific but not sensitive)
- **Recommendation**: Keep as secrets (contain sensitive workspace identifiers)

### 3. POSTGRES_URL vs Individual Parameters
- **Current**: Both approaches supported
- **Consideration**: Redundant configuration options
- **Recommendation**: Standardize on individual parameters, keep URL as optional override

### 4. Feature Cache Configuration
- **Current**: 6 separate variables for cache tuning
- **Consideration**: Could be consolidated into fewer parameters
- **Recommendation**: Keep granular control for performance tuning

## Environment-Specific Values

Variables that will have different values per environment:

### Development
- `POSTGRES_HOST`: localhost
- `POSTGRES_DB`: amprenta_dev
- `CORS_ORIGINS`: http://localhost:8501,http://localhost:3000

### Staging  
- `POSTGRES_HOST`: staging-db.internal
- `POSTGRES_DB`: amprenta_staging
- `CORS_ORIGINS`: https://staging.amprenta.com

### Production
- `POSTGRES_HOST`: prod-db.internal  
- `POSTGRES_DB`: amprenta_prod
- `CORS_ORIGINS`: https://app.amprenta.com

## Validation Requirements

### Required for Application Startup
These variables must have values or the application will fail to start:
- `OPENAI_API_KEY`
- `ZOTERO_API_KEY`
- `POSTGRES_HOST`, `POSTGRES_PORT`, `POSTGRES_DB`, `POSTGRES_USER`
- `NOTION_API_KEY` (only if `ENABLE_NOTION_SYNC=true`)

### Optional with Graceful Degradation
These variables can be empty and the application will function with reduced capabilities:
- `GEO_API_KEY` (lower NCBI rate limits)
- `BACKUP_S3_BUCKET` (no S3 backup)
- All feature flags (use defaults)

## Next Steps

1. **Create AWS Secrets Manager secrets** for each group
2. **Update application code** to read from Secrets Manager
3. **Test configuration** in development environment
4. **Deploy to staging** with new secrets management
5. **Migrate production** after staging validation
6. **Remove old environment variables** after successful migration
