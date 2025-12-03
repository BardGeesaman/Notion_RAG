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

NOTION_API_KEY = os.getenv("NOTION_API_KEY", "")
if not NOTION_API_KEY:
    raise RuntimeError(
        "NOTION_API_KEY is not set. "
        "Please define it in your environment or in a .env file. "
        "See .env.example for a template."
    )

ZOTERO_API_KEY = os.getenv("ZOTERO_API_KEY", "")
if not ZOTERO_API_KEY:
    raise RuntimeError(
        "ZOTERO_API_KEY is not set. "
        "Please define it in your environment or in a .env file. "
        "See .env.example for a template."
    )

# Database IDs (no dashes)
NOTION_EMAIL_DB_ID = "2b7adf6142ab80878ccce09b0067db60"  # Email & Notes Inbox
NOTION_RAG_DB_ID = "2bbadf6142ab8076a5fbed30f2cfcbfb"  # RAG Engine
NOTION_EXP_DATA_DB_ID = os.getenv(
    "NOTION_EXP_DATA_DB_ID", ""
)  # Experimental Data Assets
NOTION_METABOLITE_FEATURES_DB_ID = os.getenv(
    "NOTION_METABOLITE_FEATURES_DB_ID", ""
)  # Metabolite Features
NOTION_SIGNATURE_DB_ID = os.getenv("NOTION_SIGNATURE_DB_ID", "")  # Lipid Signatures
NOTION_SIGNATURE_COMPONENT_DB_ID = os.getenv(
    "NOTION_SIGNATURE_COMPONENT_DB_ID", ""
)  # Lipid Signature Components
NOTION_LIPID_SPECIES_DB_ID = os.getenv(
    "NOTION_LIPID_SPECIES_DB_ID", ""
)  # Lipid Species

# Pipeline directories
SIGNATURES_DIR = os.getenv("SIGNATURES_DIR", "")

# Signature scoring configuration
SIGNATURE_OVERLAP_THRESHOLD = float(os.getenv("SIGNATURE_OVERLAP_THRESHOLD", "0.3"))
ENABLE_SIGNATURE_SCORING = (
    os.getenv("ENABLE_SIGNATURE_SCORING", "true").lower() == "true"
)
ENABLE_LIPID_MAPPING = os.getenv("ENABLE_LIPID_MAPPING", "true").lower() == "true"


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
    api_key: str = NOTION_API_KEY
    version: str = NOTION_VERSION
    base_url: str = "https://api.notion.com/v1"
    lit_db_id: str = NOTION_LIT_DB_ID
    email_db_id: str = NOTION_EMAIL_DB_ID
    rag_db_id: str = NOTION_RAG_DB_ID
    parent_page_id: str = NOTION_PARENT_PAGE_ID
    exp_data_db_id: str = NOTION_EXP_DATA_DB_ID
    metabolite_features_db_id: str = NOTION_METABOLITE_FEATURES_DB_ID
    signature_db_id: str = NOTION_SIGNATURE_DB_ID
    signature_component_db_id: str = NOTION_SIGNATURE_COMPONENT_DB_ID
    lipid_species_db_id: str = NOTION_LIPID_SPECIES_DB_ID


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


@dataclass(frozen=True)
class AppConfig:
    openai: OpenAIConfig
    pinecone: PineconeConfig
    notion: NotionConfig
    zotero: ZoteroConfig
    pipeline: PipelineConfig


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
        )
    return _config_singleton
