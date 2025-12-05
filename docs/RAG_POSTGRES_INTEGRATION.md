# RAG Integration with Postgres - Phase 5

**Status**: Complete  
**Roadmap**: `context/UNIFIED_STRATEGIC_ROADMAP.md` Section 3.1, Phase 5

## Overview

Phase 5 integrates Postgres into the RAG system, enabling hybrid Postgres + Notion retrieval. The system now uses Postgres as the source of truth for structured data while leveraging Notion for narrative content.

## Architecture

### Hybrid Approach

- **Postgres**: Structured data (relationships, metadata, IDs)
- **Notion**: Narrative content (descriptions, results, conclusions)
- **Pinecone**: Semantic index (unchanged)

### Metadata Structure

Pinecone metadata now includes both Postgres IDs and Notion IDs:

```python
{
    "source_type": "Dataset",
    "dataset_id": "uuid-here",  # Postgres ID (primary)
    "notion_page_id": "notion-id",  # Notion ID (for narrative)
    "dataset_name": "...",
    "omics_type": "lipidomics",
    # ... other structured fields
}
```

## Components

### 1. Postgres Builder (`amprenta_rag/rag/postgres_builder.py`)

Builds RAG metadata and text from Postgres entities:

- `build_dataset_rag_metadata()` - Metadata for datasets
- `build_program_rag_metadata()` - Metadata for programs
- `build_experiment_rag_metadata()` - Metadata for experiments
- `build_signature_rag_metadata()` - Metadata for signatures
- `build_feature_rag_metadata()` - Metadata for features
- `build_dataset_rag_text()` - Text representation from Postgres
- `build_program_rag_text()` - Text representation from Postgres

### 2. Postgres Resolver (`amprenta_rag/rag/postgres_resolver.py`)

Resolves Postgres IDs and fetches context:

- `resolve_postgres_id_from_metadata()` - Extract Postgres ID from metadata
- `get_entity_type_from_metadata()` - Determine entity type
- `fetch_postgres_context()` - Fetch RAG text from Postgres
- `get_notion_id_from_postgres()` - Get Notion ID from Postgres entity

### 3. Hybrid Chunk Collection (`amprenta_rag/rag/hybrid_chunk_collection.py`)

Collects chunks from both Postgres and Notion:

- `collect_hybrid_chunks()` - Prefer Postgres, fallback to Notion
- `collect_enhanced_chunks()` - Combine Postgres + Notion for rich context

## Usage

### In RAG Queries

The RAG query system now supports Postgres by default:

```python
from amprenta_rag.query.rag.query import query_rag

# Query with Postgres support (default)
result = query_rag(
    "What datasets match the ALS signature?",
    use_postgres=True,  # Use hybrid Postgres + Notion
)

# Query with Notion-only (legacy)
result = query_rag(
    "What datasets match the ALS signature?",
    use_postgres=False,  # Use only Notion
)
```

### Building Metadata for Embedding

When embedding new data, use Postgres builders:

```python
from amprenta_rag.rag.postgres_builder import build_dataset_rag_metadata
from uuid import UUID

# Build metadata with Postgres ID
metadata = build_dataset_rag_metadata(
    dataset_id=UUID("..."),
    include_notion_id=True,  # Include Notion ID if available
)
```

### Fetching Context

When retrieving context for RAG:

```python
from amprenta_rag.rag.postgres_resolver import fetch_postgres_context
from uuid import UUID

# Fetch structured context from Postgres
context = fetch_postgres_context(
    postgres_id=UUID("..."),
    entity_type="Dataset",
    include_notion_narrative=True,  # Also fetch Notion narrative
)
```

## Benefits

1. **Structured Data**: Postgres provides reliable structured relationships
2. **Rich Context**: Notion provides narrative and human-readable content
3. **Backward Compatible**: Falls back to Notion when Postgres unavailable
4. **Flexible**: Can use Postgres-only, Notion-only, or hybrid

## Migration Path

1. **Current**: All metadata uses Notion IDs
2. **Transition**: New embeddings include both Postgres and Notion IDs
3. **Future**: Postgres IDs become primary, Notion IDs optional

## Integration Points

- **Ingestion**: Use `build_*_rag_metadata()` when embedding
- **Query**: Use `collect_hybrid_chunks()` for retrieval
- **Fallback**: Automatic fallback to Notion if Postgres unavailable

## Next Steps

Phase 6 will integrate this into ingestion pipelines, making Postgres the primary source for new data while maintaining Notion compatibility.

