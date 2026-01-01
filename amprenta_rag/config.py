# amprenta_rag/config.py

"""
Central configuration for the Amprenta RAG engine.
All IDs and constants are hard-coded for stability.
API keys are loaded from environment variables or .env file.

This module automatically loads .env file if python-dotenv is installed.
Otherwise, it falls back to environment variables.
"""

from __future__ import annotations

import os
from dataclasses import dataclass, field
from typing import List

# Import secrets management utilities
from amprenta_rag.utils.secrets import (
    get_api_key,
    get_auth_credential,
    get_integration_id,
    get_backup_config,
    get_database_url
)

# Load .env file if python-dotenv is installed (optional dependency)
try:
    from dotenv import load_dotenv

    try:
        load_dotenv()
    except Exception:
        # Skip silently if .env is not readable in sandbox/test environments
        pass
except ImportError:
    # If dotenv is not installed, silently ignore; env vars still work.
    pass


# ---------------------------------------------------------
#  🔐 REQUIRED: API Keys (from environment or .env file)
#     (validated lazily at get_config() time)
# ---------------------------------------------------------

# ---------------------------------------------------------
#  🔐 SECRETS: API Keys and Credentials (from AWS Secrets Manager or .env)
# ---------------------------------------------------------

OPENAI_API_KEY = get_api_key("openai") or ""
ZOTERO_API_KEY = get_api_key("zotero") or ""
NOTION_API_KEY = get_api_key("notion") or ""
GEO_API_KEY = get_api_key("geo") or ""
NCBI_EMAIL = get_auth_credential("ncbi_email") or ""

# ---------------------------------------------------------
#  ⚙️ NON-SENSITIVE CONFIG: Feature flags and settings
# ---------------------------------------------------------

# Pinecone deprecated - using pgvector
VECTOR_BACKEND = os.getenv("VECTOR_BACKEND", "pgvector")

# Notion API key - OPTIONAL (only required if ENABLE_NOTION_SYNC=true)
ENABLE_NOTION_SYNC_ENV = os.getenv("ENABLE_NOTION_SYNC", "false").lower() == "true"

# API server configuration
def _parse_cors_origins() -> List[str]:
    raw = os.getenv("CORS_ORIGINS", "http://localhost:8501")
    return [o.strip() for o in raw.split(",") if o.strip()]

# ---------------------------------------------------------
#  🗄️ DATABASE CONFIG: Postgres connection (secrets + config)
# ---------------------------------------------------------

# Database URL (constructed securely from secrets)
POSTGRES_URL = get_database_url()

# Individual connection parameters (for non-URL usage)
POSTGRES_HOST = os.getenv("POSTGRES_HOST", "localhost")
POSTGRES_PORT = int(os.getenv("POSTGRES_PORT", "5432"))
POSTGRES_DB = os.getenv("POSTGRES_DB", "amprenta")
POSTGRES_USER = os.getenv("POSTGRES_USER", "bard")
POSTGRES_ECHO = os.getenv("POSTGRES_ECHO", "false").lower() == "true"

# Note: POSTGRES_PASSWORD is handled by secrets.py - not exposed here

# ---------------------------------------------------------
#  🔗 INTEGRATION IDs: Notion workspace identifiers (from secrets)
# ---------------------------------------------------------

# Static database IDs (hardcoded for stability)
NOTION_EMAIL_DB_ID = "2b7adf6142ab80878ccce09b0067db60"  # Email & Notes Inbox
NOTION_RAG_DB_ID = "2bbadf6142ab8076a5fbed30f2cfcbfb"  # RAG Engine

# Dynamic database IDs (from secrets for environment-specific workspaces)
NOTION_EXP_DATA_DB_ID = get_integration_id("exp_data") or ""
NOTION_METABOLITE_FEATURES_DB_ID = get_integration_id("metabolite_features") or ""
NOTION_PROTEIN_FEATURES_DB_ID = get_integration_id("protein_features") or ""
NOTION_GENE_FEATURES_DB_ID = get_integration_id("gene_features") or ""
NOTION_SIGNATURE_DB_ID = get_integration_id("signatures") or ""
NOTION_SIGNATURE_COMPONENT_DB_ID = get_integration_id("signature_components") or ""
NOTION_LIPID_SPECIES_DB_ID = get_integration_id("lipid_species") or ""
NOTION_PROGRAMS_DB_ID = get_integration_id("programs") or ""
NOTION_EXPERIMENTS_DB_ID = get_integration_id("experiments") or ""
NOTION_COMPOUND_FEATURES_DB_ID = get_integration_id("compound_features") or ""
NOTION_HTS_CAMPAIGNS_DB_ID = get_integration_id("hts_campaigns") or ""
NOTION_BIOCHEMICAL_HITS_DB_ID = get_integration_id("biochemical_hits") or ""
NOTION_PATHWAYS_DB_ID = get_integration_id("pathways") or ""

# Pipeline directories
SIGNATURES_DIR = os.getenv("SIGNATURES_DIR", "")

# Signature scoring configuration
SIGNATURE_OVERLAP_THRESHOLD = float(os.getenv("SIGNATURE_OVERLAP_THRESHOLD", "0.3"))
ENABLE_SIGNATURE_SCORING = (
    os.getenv("ENABLE_SIGNATURE_SCORING", "true").lower() == "true"
)
ENABLE_LIPID_MAPPING = os.getenv("ENABLE_LIPID_MAPPING", "true").lower() == "true"
ENABLE_FEATURE_LINKING = os.getenv("ENABLE_FEATURE_LINKING", "true").lower() == "true"
FEATURE_LINKING_MAX_WORKERS = int(os.getenv("FEATURE_LINKING_MAX_WORKERS", "10"))
ENABLE_LLM_SEMANTIC_EXTRACTION = os.getenv("ENABLE_LLM_SEMANTIC_EXTRACTION", "false").lower() == "true"

# Feature cache configuration
FEATURE_CACHE_ENABLED = os.getenv("FEATURE_CACHE_ENABLED", "true").lower() == "true"

# Prefer FEATURE_CACHE_TTL but fall back to legacy FEATURE_CACHE_TTL_SECONDS
_FEATURE_CACHE_TTL_RAW = os.getenv("FEATURE_CACHE_TTL")
if _FEATURE_CACHE_TTL_RAW is None:
    _FEATURE_CACHE_TTL_RAW = os.getenv("FEATURE_CACHE_TTL_SECONDS", "3600")
FEATURE_CACHE_TTL = int(_FEATURE_CACHE_TTL_RAW)

# Max size defaults to 1000 datasets unless explicitly overridden
_FEATURE_CACHE_MAX_SIZE_RAW = os.getenv("FEATURE_CACHE_MAX_SIZE")
if _FEATURE_CACHE_MAX_SIZE_RAW is None or _FEATURE_CACHE_MAX_SIZE_RAW == "":
    FEATURE_CACHE_MAX_SIZE = 1000
else:
    FEATURE_CACHE_MAX_SIZE = int(_FEATURE_CACHE_MAX_SIZE_RAW)

# Additional cache tuning settings (used by enhanced cache paths)
FEATURE_CACHE_ENABLE_PERSISTENCE = os.getenv("FEATURE_CACHE_ENABLE_PERSISTENCE", "true").lower() == "true"
FEATURE_CACHE_DIR = os.getenv("FEATURE_CACHE_DIR", "")
FEATURE_CACHE_PARALLEL_WORKERS = int(os.getenv("FEATURE_CACHE_PARALLEL_WORKERS", "5"))

# Backward-compatibility alias
FEATURE_CACHE_TTL_SECONDS = FEATURE_CACHE_TTL

# Auto-linking configuration
AUTO_LINK_ENABLED = os.getenv("AUTO_LINK_ENABLED", "true").lower() == "true"
AUTO_LINK_MIN_CONFIDENCE = float(os.getenv("AUTO_LINK_MIN_CONFIDENCE", "0.8"))


# ---------------------------------------------------------
#  📚 CONSTANT VALUES (YOUR REAL IDs)
# ---------------------------------------------------------

# Notion database IDs (NO DASHES)
NOTION_LIT_DB_ID = "c232e83e9c894e0ab1ac8a2a2e9d20a4"
NOTION_EMAIL_DB_ID = "2b7adf6142ab80878ccce09b0067db60"
NOTION_RAG_DB_ID = "2bbadf6142ab8076a5fbed30f2cfcbfb"

# Notion parent page (NO DASHES)
NOTION_PARENT_PAGE_ID = "2b9adf6142ab8026853ef58f725665a6"

# Zotero library info
ZOTERO_LIBRARY_TYPE = "group"
ZOTERO_LIBRARY_ID = "6259755"

# Pinecone config (deprecated - kept for backwards compatibility)
PINECONE_INDEX_NAME = "amprenta-rag"
PINECONE_REGION = "us-east-1"
PINECONE_NAMESPACE = "default"

# OpenAI model defaults
OPENAI_CHAT_MODEL = "gpt-4.1-mini"
OPENAI_EMBED_MODEL = "text-embedding-3-large"

# Notion API version
NOTION_VERSION = "2022-06-28"

# ---------------------------------------------------------
#  📦 CONFIG STRUCTURES
# ---------------------------------------------------------


@dataclass(frozen=True)
class OpenAIConfig:
    api_key: str = OPENAI_API_KEY
    model: str = OPENAI_CHAT_MODEL
    embedding_model: str = OPENAI_EMBED_MODEL


@dataclass(frozen=True)
class PineconeConfig:
    """Deprecated: Pinecone support removed. Kept for backwards compatibility."""
    api_key: str = ""
    index_name: str = PINECONE_INDEX_NAME
    region: str = PINECONE_REGION
    namespace: str = PINECONE_NAMESPACE


@dataclass(frozen=True)
class NotionConfig:
    # Notion API key - may be empty if Notion sync is disabled
    api_key: str = NOTION_API_KEY
    version: str = NOTION_VERSION
    # Kept for historical compatibility; Notion is deprecated and not used.
    base_url: str = ""
    lit_db_id: str = NOTION_LIT_DB_ID
    email_db_id: str = NOTION_EMAIL_DB_ID
    rag_db_id: str = NOTION_RAG_DB_ID
    parent_page_id: str = NOTION_PARENT_PAGE_ID
    exp_data_db_id: str = NOTION_EXP_DATA_DB_ID
    metabolite_features_db_id: str = NOTION_METABOLITE_FEATURES_DB_ID
    protein_features_db_id: str = NOTION_PROTEIN_FEATURES_DB_ID
    gene_features_db_id: str = NOTION_GENE_FEATURES_DB_ID
    signature_db_id: str = NOTION_SIGNATURE_DB_ID
    signature_component_db_id: str = NOTION_SIGNATURE_COMPONENT_DB_ID
    lipid_species_db_id: str = NOTION_LIPID_SPECIES_DB_ID
    programs_db_id: str = NOTION_PROGRAMS_DB_ID
    experiments_db_id: str = NOTION_EXPERIMENTS_DB_ID
    compound_features_db_id: str = NOTION_COMPOUND_FEATURES_DB_ID
    hts_campaigns_db_id: str = NOTION_HTS_CAMPAIGNS_DB_ID
    biochemical_hits_db_id: str = NOTION_BIOCHEMICAL_HITS_DB_ID
    pathways_db_id: str = NOTION_PATHWAYS_DB_ID


@dataclass(frozen=True)
class ZoteroConfig:
    api_key: str = ZOTERO_API_KEY
    library_type: str = ZOTERO_LIBRARY_TYPE
    library_id: str = ZOTERO_LIBRARY_ID


@dataclass(frozen=True)
class PipelineConfig:
    signatures_dir: str = SIGNATURES_DIR
    signature_overlap_threshold: float = SIGNATURE_OVERLAP_THRESHOLD
    enable_signature_scoring: bool = ENABLE_SIGNATURE_SCORING
    enable_lipid_mapping: bool = ENABLE_LIPID_MAPPING
    enable_feature_linking: bool = ENABLE_FEATURE_LINKING
    enable_llm_semantic_extraction: bool = ENABLE_LLM_SEMANTIC_EXTRACTION
    feature_linking_max_workers: int = FEATURE_LINKING_MAX_WORKERS
    auto_link_enabled: bool = AUTO_LINK_ENABLED
    auto_link_min_confidence: float = AUTO_LINK_MIN_CONFIDENCE
    # TIER 3: Postgres as Source of Truth (PRIMARY for performance)
    # Set to true to use Postgres as primary database (recommended for bulk ingestion)
    use_postgres_as_sot: bool = os.getenv("USE_POSTGRES_AS_SOT", "true").lower() == "true"
    # Optional: Sync to Notion for documentation (set to false to disable Notion entirely)
    enable_notion_sync: bool = os.getenv("ENABLE_NOTION_SYNC", "false").lower() == "true"
    # Dual-write mode (transition period only - writes to both Postgres and Notion)
    enable_dual_write: bool = os.getenv("ENABLE_DUAL_WRITE", "false").lower() == "true"


@dataclass(frozen=True)
class FeatureCacheConfig:
    enabled: bool = FEATURE_CACHE_ENABLED
    ttl: int = FEATURE_CACHE_TTL
    max_size: int = FEATURE_CACHE_MAX_SIZE


@dataclass(frozen=True)
class ServerConfig:
    cors_origins: List[str] = field(default_factory=_parse_cors_origins)


@dataclass(frozen=True)
class PostgresConfig:
    """Postgres database configuration for TIER 3 architecture evolution."""
    url: str = POSTGRES_URL
    host: str = POSTGRES_HOST
    port: int = POSTGRES_PORT
    db: str = POSTGRES_DB
    user: str = POSTGRES_USER
    echo: bool = POSTGRES_ECHO
    # Note: password is handled securely via secrets.py and included in url


# ---------------------------------------------------------
#  💾 BACKUP CONFIG: Backup and encryption settings
# ---------------------------------------------------------

# Non-sensitive backup configuration
BACKUP_S3_ENABLED = os.getenv("BACKUP_S3_ENABLED", "false").lower() == "true"
BACKUP_S3_BUCKET = os.getenv("BACKUP_S3_BUCKET", "")
BACKUP_LOCAL_DIR = os.getenv("BACKUP_LOCAL_DIR", "./backups")
BACKUP_RETENTION_DAYS = int(os.getenv("BACKUP_RETENTION_DAYS", "365"))

# Sensitive backup configuration (from secrets)
BACKUP_KMS_KEY_ID = get_backup_config("kms_key_id") or ""


@dataclass(frozen=True)
class BackupConfig:
    """Backup and disaster recovery configuration."""
    s3_enabled: bool = BACKUP_S3_ENABLED
    s3_bucket: str = BACKUP_S3_BUCKET
    kms_key_id: str = BACKUP_KMS_KEY_ID
    local_dir: str = BACKUP_LOCAL_DIR
    retention_days: int = BACKUP_RETENTION_DAYS


@dataclass(frozen=True)
class AppConfig:
    openai: OpenAIConfig
    pinecone: PineconeConfig
    notion: NotionConfig
    zotero: ZoteroConfig
    pipeline: PipelineConfig
    feature_cache: FeatureCacheConfig
    postgres: PostgresConfig
    server: ServerConfig
    backup: BackupConfig
    vector_backend: str = VECTOR_BACKEND


# ---------------------------------------------------------
#  🔧 SINGLETON ACCESSOR
# ---------------------------------------------------------

_config_singleton: AppConfig | None = None


def _validate_required_keys():
    missing = []
    if not OPENAI_API_KEY:
        missing.append("OPENAI_API_KEY")
    # Pinecone deprecated - no longer required
    if not ZOTERO_API_KEY:
        missing.append("ZOTERO_API_KEY")
    if ENABLE_NOTION_SYNC_ENV and not NOTION_API_KEY:
        missing.append("NOTION_API_KEY (required because ENABLE_NOTION_SYNC=true)")
    if missing:
        raise RuntimeError(
            "Missing required configuration values: " + ", ".join(missing)
        )


def get_config() -> AppConfig:
    global _config_singleton
    if _config_singleton is None:
        _validate_required_keys()
        _config_singleton = AppConfig(
            openai=OpenAIConfig(),
            pinecone=PineconeConfig(),
            notion=NotionConfig(),
            zotero=ZoteroConfig(),
            pipeline=PipelineConfig(),
            feature_cache=FeatureCacheConfig(),
            postgres=PostgresConfig(),
            server=ServerConfig(),
            backup=BackupConfig(),
            vector_backend=VECTOR_BACKEND,
        )
    return _config_singleton
