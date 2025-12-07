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
from dataclasses import dataclass
from typing import Optional
from typing import Optional
from typing import Optional

# Load .env file if python-dotenv is installed (optional dependency)
try:
    from dotenv import load_dotenv

    load_dotenv()
except ImportError:
    # If dotenv is not installed, silently ignore; env vars still work.
    pass


# ---------------------------------------------------------
#  🔐 REQUIRED: API Keys (from environment or .env file)
# ---------------------------------------------------------

OPENAI_API_KEY = os.getenv("OPENAI_API_KEY", "")
if not OPENAI_API_KEY:
    raise RuntimeError(
        "OPENAI_API_KEY is not set. "
        "Please define it in your environment or in a .env file. "
        "See .env.example for a template."
    )

PINECONE_API_KEY = os.getenv("PINECONE_API_KEY", "")
if not PINECONE_API_KEY:
    raise RuntimeError(
        "PINECONE_API_KEY is not set. "
        "Please define it in your environment or in a .env file. "
        "See .env.example for a template."
    )

# Notion API key - OPTIONAL (only required if ENABLE_NOTION_SYNC=true)
# Check if Notion sync is enabled first, then validate API key if needed
ENABLE_NOTION_SYNC_ENV = os.getenv("ENABLE_NOTION_SYNC", "false").lower() == "true"
NOTION_API_KEY = os.getenv("NOTION_API_KEY", "")

# Only require Notion API key if Notion sync is explicitly enabled
if ENABLE_NOTION_SYNC_ENV and not NOTION_API_KEY:
    raise RuntimeError(
        "NOTION_API_KEY is not set but ENABLE_NOTION_SYNC is true. "
        "Either set NOTION_API_KEY in your environment/.env file, "
        "or set ENABLE_NOTION_SYNC=false to disable Notion entirely. "
        "Notion is now optional and disabled by default."
    )

ZOTERO_API_KEY = os.getenv("ZOTERO_API_KEY", "")
if not ZOTERO_API_KEY:
    raise RuntimeError(
        "ZOTERO_API_KEY is not set. "
        "Please define it in your environment or in a .env file. "
        "See .env.example for a template."
    )

# Optional API keys for public repositories
GEO_API_KEY = os.getenv("GEO_API_KEY", "")  # Optional, for higher NCBI rate limits
NCBI_EMAIL = os.getenv("NCBI_EMAIL", "")  # Required by NCBI for API access

# Postgres Database Configuration (optional - for TIER 3 architecture evolution)
POSTGRES_URL = os.getenv("POSTGRES_URL", "")  # Full connection string (optional)
POSTGRES_HOST = os.getenv("POSTGRES_HOST", "localhost")
POSTGRES_PORT = int(os.getenv("POSTGRES_PORT", "5432"))
POSTGRES_DB = os.getenv("POSTGRES_DB", "amprenta_rag")
POSTGRES_USER = os.getenv("POSTGRES_USER", "postgres")
POSTGRES_PASSWORD = os.getenv("POSTGRES_PASSWORD", "")
POSTGRES_ECHO = os.getenv("POSTGRES_ECHO", "false").lower() == "true"

# Database IDs (no dashes)
NOTION_EMAIL_DB_ID = "2b7adf6142ab80878ccce09b0067db60"  # Email & Notes Inbox
NOTION_RAG_DB_ID = "2bbadf6142ab8076a5fbed30f2cfcbfb"  # RAG Engine
NOTION_EXP_DATA_DB_ID = os.getenv(
    "NOTION_EXP_DATA_DB_ID", ""
)  # Experimental Data Assets
NOTION_METABOLITE_FEATURES_DB_ID = os.getenv(
    "NOTION_METABOLITE_FEATURES_DB_ID", ""
)  # Metabolite Features
NOTION_PROTEIN_FEATURES_DB_ID = os.getenv(
    "NOTION_PROTEIN_FEATURES_DB_ID", ""
)  # Protein Features
NOTION_GENE_FEATURES_DB_ID = os.getenv(
    "NOTION_GENE_FEATURES_DB_ID", ""
)  # Gene Features
NOTION_SIGNATURE_DB_ID = os.getenv("NOTION_SIGNATURE_DB_ID", "")  # Lipid Signatures
NOTION_SIGNATURE_COMPONENT_DB_ID = os.getenv(
    "NOTION_SIGNATURE_COMPONENT_DB_ID", ""
)  # Lipid Signature Components
NOTION_LIPID_SPECIES_DB_ID = os.getenv(
    "NOTION_LIPID_SPECIES_DB_ID", ""
)  # Lipid Species
NOTION_PROGRAMS_DB_ID = os.getenv(
    "NOTION_PROGRAMS_DB_ID", ""
)  # Programs
NOTION_EXPERIMENTS_DB_ID = os.getenv(
    "NOTION_EXPERIMENTS_DB_ID", ""
)  # Experiments
NOTION_COMPOUND_FEATURES_DB_ID = os.getenv(
    "NOTION_COMPOUND_FEATURES_DB_ID", ""
)  # Compound Features
NOTION_HTS_CAMPAIGNS_DB_ID = os.getenv(
    "NOTION_HTS_CAMPAIGNS_DB_ID", ""
)  # HTS Campaigns
NOTION_BIOCHEMICAL_HITS_DB_ID = os.getenv(
    "NOTION_BIOCHEMICAL_HITS_DB_ID", ""
)  # Biochemical Hits
NOTION_PATHWAYS_DB_ID = os.getenv(
    "NOTION_PATHWAYS_DB_ID", ""
)  # Pathways

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

# Pinecone index config
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
    api_key: str = PINECONE_API_KEY
    index_name: str = PINECONE_INDEX_NAME
    region: str = PINECONE_REGION
    namespace: str = PINECONE_NAMESPACE


@dataclass(frozen=True)
class NotionConfig:
    # Notion API key - may be empty if Notion sync is disabled
    api_key: str = NOTION_API_KEY
    version: str = NOTION_VERSION
    base_url: str = "https://api.notion.com/v1"
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
class PostgresConfig:
    """Postgres database configuration for TIER 3 architecture evolution."""
    url: str = POSTGRES_URL
    host: str = POSTGRES_HOST
    port: int = POSTGRES_PORT
    db: str = POSTGRES_DB
    user: str = POSTGRES_USER
    password: str = POSTGRES_PASSWORD
    echo: bool = POSTGRES_ECHO


@dataclass(frozen=True)
class AppConfig:
    openai: OpenAIConfig
    pinecone: PineconeConfig
    notion: NotionConfig
    zotero: ZoteroConfig
    pipeline: PipelineConfig
    feature_cache: FeatureCacheConfig
    postgres: PostgresConfig


# ---------------------------------------------------------
#  🔧 SINGLETON ACCESSOR
# ---------------------------------------------------------

_config_singleton: AppConfig | None = None


def get_config() -> AppConfig:
    global _config_singleton
    if _config_singleton is None:
        _config_singleton = AppConfig(
            openai=OpenAIConfig(),
            pinecone=PineconeConfig(),
            notion=NotionConfig(),
            zotero=ZoteroConfig(),
            pipeline=PipelineConfig(),
            feature_cache=FeatureCacheConfig(),
            postgres=PostgresConfig(),
        )
    return _config_singleton
