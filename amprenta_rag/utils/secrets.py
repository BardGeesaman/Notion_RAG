"""AWS Secrets Manager integration with local fallback.

This module provides seamless integration between AWS Secrets Manager (for production)
and environment variables (for local development). It automatically detects the
runtime environment and uses the appropriate secret retrieval method.

Usage:
    # Get database URL (automatically constructed)
    db_url = get_database_url()
    
    # Get API keys
    openai_key = get_api_key("openai")
    zotero_key = get_api_key("zotero")
    
    # Get specific secret values
    email = get_secret("amprenta/prod/auth", "NCBI_EMAIL")
"""

import os
import json
import logging
from typing import Optional, Dict, Any

logger = logging.getLogger(__name__)

# Check if running in AWS (ECS sets these environment variables)
IS_AWS = (
    os.getenv("AWS_EXECUTION_ENV") is not None or 
    os.getenv("ECS_CONTAINER_METADATA_URI") is not None or
    os.getenv("AWS_REGION") is not None
)

# Global cache for secrets to avoid repeated AWS API calls
_secrets_cache: Dict[str, str] = {}


def get_secret(secret_name: str, key: Optional[str] = None) -> Optional[str]:
    """
    Get secret from AWS Secrets Manager or fall back to environment variable.
    
    Args:
        secret_name: AWS secret name (e.g., "amprenta/prod/database")
        key: Optional JSON key within the secret
    
    Returns:
        Secret value or None if not found
        
    Examples:
        # Get entire secret (if it's a simple string)
        password = get_secret("amprenta/prod/simple-password")
        
        # Get specific key from JSON secret
        db_password = get_secret("amprenta/prod/database", "POSTGRES_PASSWORD")
    """
    if IS_AWS:
        return _get_aws_secret(secret_name, key)
    else:
        # Local development: use environment variables
        # Convert secret name and key to environment variable format
        if key:
            env_key = key.upper().replace("-", "_")
        else:
            env_key = secret_name.upper().replace("/", "_").replace("-", "_")
        
        value = os.getenv(env_key)
        if value:
            logger.debug(f"Retrieved local env var: {env_key}")
        return value


def _get_aws_secret(secret_name: str, key: Optional[str] = None) -> Optional[str]:
    """Fetch secret from AWS Secrets Manager with caching."""
    global _secrets_cache
    
    # Check cache first to avoid repeated API calls
    cache_key = f"{secret_name}:{key or 'full'}"
    if cache_key in _secrets_cache:
        return _secrets_cache[cache_key]
    
    # Fetch from AWS Secrets Manager
    try:
        import boto3
        from botocore.exceptions import ClientError, NoCredentialsError
        
        client = boto3.client("secretsmanager")
        response = client.get_secret_value(SecretId=secret_name)
        secret_value = response.get("SecretString")
        
        if not secret_value:
            logger.warning(f"Secret {secret_name} is empty")
            return None
            
        logger.debug(f"Retrieved AWS secret: {secret_name}")
        
        # Parse JSON if key is specified
        if key and secret_value:
            try:
                parsed = json.loads(secret_value)
                value = parsed.get(key)
                _secrets_cache[cache_key] = value
                return value
            except json.JSONDecodeError as e:
                logger.error(f"Failed to parse JSON in secret {secret_name}: {e}")
                return None
        else:
            # Return full secret value
            _secrets_cache[cache_key] = secret_value
            return secret_value
            
    except ClientError as e:
        error_code = e.response.get('Error', {}).get('Code')
        if error_code == 'ResourceNotFoundException':
            logger.error(f"Secret not found: {secret_name}")
        elif error_code == 'UnauthorizedOperation':
            logger.error(f"Insufficient permissions for secret: {secret_name}")
        else:
            logger.error(f"AWS error retrieving secret {secret_name}: {e}")
        return None
        
    except NoCredentialsError:
        logger.error("AWS credentials not configured")
        return None
        
    except ImportError:
        logger.warning("boto3 not available, falling back to env vars")
        # Fallback to environment variable
        env_key = key.upper() if key else secret_name.upper().replace("/", "_").replace("-", "_")
        return os.getenv(env_key)
    
    except Exception as e:
        logger.error(f"Unexpected error retrieving secret {secret_name}: {e}")
        return None


def get_database_url() -> str:
    """
    Get database connection URL.
    
    In AWS: Uses POSTGRES_URL from database secret group
    Locally: Constructs URL from individual environment variables
    
    Returns:
        PostgreSQL connection URL
    """
    # Determine environment for secret name
    environment = os.getenv("ENVIRONMENT", "dev")
    
    if IS_AWS:
        # Try to get full URL first (preferred)
        postgres_url = get_secret(f"amprenta/{environment}/database", "POSTGRES_URL")
        if postgres_url:
            return postgres_url
            
        # Fallback: construct from individual components
        user = os.getenv("POSTGRES_USER", "postgres")
        password = get_secret(f"amprenta/{environment}/database", "POSTGRES_PASSWORD")
        host = os.getenv("POSTGRES_HOST", "localhost")
        port = os.getenv("POSTGRES_PORT", "5432")
        db = os.getenv("POSTGRES_DB", "amprenta")
        
        if not password:
            logger.warning("POSTGRES_PASSWORD not found in AWS secrets")
            password = ""
            
        return f"postgresql://{user}:{password}@{host}:{port}/{db}"
    else:
        # Local development: construct from environment variables
        user = os.getenv("POSTGRES_USER", "postgres")
        password = os.getenv("POSTGRES_PASSWORD", "")
        host = os.getenv("POSTGRES_HOST", "localhost")
        port = os.getenv("POSTGRES_PORT", "5432")
        db = os.getenv("POSTGRES_DB", "amprenta")
        
        return f"postgresql://{user}:{password}@{host}:{port}/{db}"


def get_api_key(service: str) -> Optional[str]:
    """
    Get API key for external service.
    
    Args:
        service: Service name (openai, zotero, notion, geo)
        
    Returns:
        API key or None if not found
        
    Examples:
        openai_key = get_api_key("openai")
        zotero_key = get_api_key("zotero")
    """
    key_mapping = {
        "openai": "OPENAI_API_KEY",
        "zotero": "ZOTERO_API_KEY", 
        "notion": "NOTION_API_KEY",
        "geo": "GEO_API_KEY",
    }
    
    key_name = key_mapping.get(service.lower())
    if not key_name:
        logger.warning(f"Unknown service: {service}")
        return None
    
    # Determine environment for secret name
    environment = os.getenv("ENVIRONMENT", "dev")
    
    if IS_AWS:
        return get_secret(f"amprenta/{environment}/api-keys", key_name)
    else:
        return os.getenv(key_name)


def get_auth_credential(credential: str) -> Optional[str]:
    """
    Get authentication credential.
    
    Args:
        credential: Credential name (ncbi_email, etc.)
        
    Returns:
        Credential value or None if not found
    """
    credential_mapping = {
        "ncbi_email": "NCBI_EMAIL",
    }
    
    key_name = credential_mapping.get(credential.lower())
    if not key_name:
        logger.warning(f"Unknown credential: {credential}")
        return None
    
    # Determine environment for secret name
    environment = os.getenv("ENVIRONMENT", "dev")
    
    if IS_AWS:
        return get_secret(f"amprenta/{environment}/auth", key_name)
    else:
        return os.getenv(key_name)


def get_integration_id(integration: str) -> Optional[str]:
    """
    Get integration identifier (e.g., Notion database IDs).
    
    Args:
        integration: Integration name (e.g., "experiments", "compounds")
        
    Returns:
        Integration ID or None if not found
    """
    # Map integration names to environment variable names
    integration_mapping = {
        "exp_data": "NOTION_EXP_DATA_DB_ID",
        "metabolite_features": "NOTION_METABOLITE_FEATURES_DB_ID",
        "protein_features": "NOTION_PROTEIN_FEATURES_DB_ID",
        "gene_features": "NOTION_GENE_FEATURES_DB_ID",
        "signatures": "NOTION_SIGNATURE_DB_ID",
        "signature_components": "NOTION_SIGNATURE_COMPONENT_DB_ID",
        "lipid_species": "NOTION_LIPID_SPECIES_DB_ID",
        "programs": "NOTION_PROGRAMS_DB_ID",
        "experiments": "NOTION_EXPERIMENTS_DB_ID",
        "compound_features": "NOTION_COMPOUND_FEATURES_DB_ID",
        "hts_campaigns": "NOTION_HTS_CAMPAIGNS_DB_ID",
        "biochemical_hits": "NOTION_BIOCHEMICAL_HITS_DB_ID",
        "pathways": "NOTION_PATHWAYS_DB_ID",
    }
    
    key_name = integration_mapping.get(integration.lower())
    if not key_name:
        logger.warning(f"Unknown integration: {integration}")
        return None
    
    # Determine environment for secret name
    environment = os.getenv("ENVIRONMENT", "dev")
    
    if IS_AWS:
        return get_secret(f"amprenta/{environment}/integrations", key_name)
    else:
        return os.getenv(key_name)


def get_backup_config(config: str) -> Optional[str]:
    """
    Get backup configuration value.
    
    Args:
        config: Configuration name (kms_key_id, etc.)
        
    Returns:
        Configuration value or None if not found
    """
    config_mapping = {
        "kms_key_id": "BACKUP_KMS_KEY_ID",
    }
    
    key_name = config_mapping.get(config.lower())
    if not key_name:
        logger.warning(f"Unknown backup config: {config}")
        return None
    
    # Determine environment for secret name
    environment = os.getenv("ENVIRONMENT", "dev")
    
    if IS_AWS:
        return get_secret(f"amprenta/{environment}/backup", key_name)
    else:
        return os.getenv(key_name)


def clear_cache() -> None:
    """
    Clear the secrets cache.
    
    Useful for testing or when secrets are rotated and need to be refreshed.
    """
    global _secrets_cache
    _secrets_cache.clear()
    logger.debug("Secrets cache cleared")


def get_environment_info() -> Dict[str, Any]:
    """
    Get information about the current environment and secrets configuration.
    
    Returns:
        Dictionary with environment information
    """
    return {
        "is_aws": IS_AWS,
        "environment": os.getenv("ENVIRONMENT", "dev"),
        "aws_region": os.getenv("AWS_REGION"),
        "execution_env": os.getenv("AWS_EXECUTION_ENV"),
        "ecs_metadata_uri": os.getenv("ECS_CONTAINER_METADATA_URI"),
        "cache_size": len(_secrets_cache),
    }


# Backward compatibility aliases
def get_database_password() -> Optional[str]:
    """Get database password (deprecated - use get_database_url instead)."""
    environment = os.getenv("ENVIRONMENT", "dev")
    if IS_AWS:
        return get_secret(f"amprenta/{environment}/database", "POSTGRES_PASSWORD")
    else:
        return os.getenv("POSTGRES_PASSWORD")


def get_openai_key() -> Optional[str]:
    """Get OpenAI API key (deprecated - use get_api_key('openai') instead)."""
    return get_api_key("openai")


def get_zotero_key() -> Optional[str]:
    """Get Zotero API key (deprecated - use get_api_key('zotero') instead)."""
    return get_api_key("zotero")
