# amprenta_rag/__init__.py

"""
Core RAG engine for the Amprenta Notion workspace.

This package centralizes:
- config (API keys, DB IDs, index names)
- shared API clients (OpenAI, Pinecone, Notion)
- query logic and, later, ingestion/sync tools.
"""

__all__ = ["config"]