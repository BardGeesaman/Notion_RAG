# Architecture Overview

**Last Updated**: 2025-12-03

This document provides a high-level overview of the Amprenta RAG multi-omics platform architecture.

---

## System Overview

The Amprenta platform is a unified AI-native multi-omics knowledge system that integrates:

- **Notion**: Canonical knowledge graph and schema
- **Pinecone**: Vector index for semantic search (RAG)
- **OpenAI**: Embeddings and LLM reasoning
- **Python**: Ingestion pipelines and query engine

```
┌─────────────────────────────────────────────────────────────┐
│                    Amprenta RAG Platform                    │
└─────────────────────────────────────────────────────────────┘
                              │
        ┌─────────────────────┼─────────────────────┐
        │                     │                     │
        ▼                     ▼                     ▼
   ┌─────────┐          ┌─────────┐         ┌─────────┐
   │ Notion  │          │ Pinecone │         │ OpenAI  │
   │   DB    │          │  Vector  │         │   API   │
   │         │          │  Index   │         │         │
   └─────────┘          └─────────┘         └─────────┘
        │                     │                     │
        └─────────────────────┼─────────────────────┘
                              │
                    ┌─────────▼─────────┐
                    │  Python Pipeline │
                    │  (Ingestion +     │
                    │   Query Engine)  │
                    └──────────────────┘
```

---

## Core Components

### 1. Ingestion Layer

**Purpose**: Ingest multi-omics data and create knowledge graph in Notion.

**Modules**:

- **Omics Ingestion** (`lipidomics_ingestion.py`, `metabolomics_ingestion.py`, etc.)
  - Parse CSV/TSV files
  - Normalize feature names
  - Create/update dataset pages in Notion
  - Link features to datasets

- **Signature Ingestion** (`signature_ingestion.py`)
  - Load signature definitions (TSV/CSV)
  - Create signature, component, and species pages
  - Link components to features

- **Feature Linking** (`features/linking.py`)
  - Create/find feature pages (Gene, Protein, Metabolite, Lipid Species)
  - Link features to datasets via relations
  - Maintain idempotency

- **Signature Matching** (`signature_matching/`)
  - Score datasets against signatures
  - Find matches above threshold
  - Update dataset pages with matches

- **RAG Embedding** (`text_embedding_utils.py`)
  - Build text representations
  - Chunk text into ~2000 token blocks
  - Embed via OpenAI
  - Upsert to Pinecone

**Data Flow**:
```
CSV/TSV File
    │
    ▼
Parse & Normalize Features
    │
    ▼
Create/Update Dataset Page (Notion)
    │
    ├──► Link Features to Dataset
    ├──► Score Against Signatures
    └──► Embed & Upsert to Pinecone
```

---

### 2. Query Layer

**Purpose**: Retrieve and reason over multi-omics knowledge.

**Modules**:

- **RAG Engine** (`rag_engine.py`)
  - Semantic search in Pinecone
  - Retrieve relevant chunks
  - Synthesize answers with LLM

- **Cross-Omics Reasoning** (`cross_omics/`)
  - Generate high-level summaries
  - Analyze multi-omics patterns
  - Cross-reference features across omics types

- **Signature Similarity** (`signature_matching/matching.py`)
  - Find matching signatures for datasets
  - Rank by score and overlap

**Query Flow**:
```
User Query
    │
    ▼
Semantic Search (Pinecone)
    │
    ├──► Retrieve Chunks
    ├──► Fetch Notion Pages
    └──► Synthesize with LLM
```

---

### 3. Feature Management

**Purpose**: Normalize and link omics features across the knowledge graph.

**Modules**:

- **Normalization** (`features/normalization.py`)
  - Metabolite name normalization
  - Protein/gene identifier cleaning
  - Lipid species canonicalization

- **Extraction** (`features/extraction.py`)
  - Extract features from mwTab JSON
  - Scan text for metabolite mentions

- **Linking** (`features/linking.py`)
  - Create/find feature pages
  - Link features to datasets
  - Maintain bidirectional relations

**Feature Types**:

- **Genes**: HGNC symbols, Ensembl IDs
- **Proteins**: UniProt IDs, canonical gene names
- **Metabolites**: Normalized metabolite names (HMDB, KEGG)
- **Lipids**: Canonical species format (Cer(d18:1/16:0))

---

### 4. Signature System

**Purpose**: Define and match multi-omics signatures.

**Components**:

- **Signature**: Named collection of components
- **Component**: Individual feature with direction (↑/↓) and weight
- **Species**: Canonical lipid species representation

**Signature Structure**:
```
Signature: "ALS-CSF-Core-6Ceramides"
├── Component: Cer(d18:1/16:0) ↑ 1.0
├── Component: Cer(d18:1/18:0) ↑ 0.8
├── Component: SM(d18:1/16:0) ↓ 0.5
└── ...
```

**Scoring Algorithm**:
1. Extract features from dataset
2. Match signature components to dataset features
3. Calculate overlap fraction (matched / total)
4. Apply direction consistency penalty
5. Weight by component weights
6. Return score and detailed metrics

---

### 5. Caching System

**Purpose**: Improve performance by caching dataset features.

**Module**: `dataset_feature_cache.py`

**Features**:

- In-memory cache with TTL (default: 1 hour)
- Singleton pattern
- Preloading for batch operations
- Automatic expiry

**Cache Structure**:
```python
{
    "dataset-page-id": {
        "features_by_type": {
            "gene": {"TP53", "TNF", ...},
            "protein": {"P04637", ...},
            "metabolite": {"Glutamate", ...},
            "lipid": {"Cer(d18:1/16:0)", ...}
        },
        "timestamp": datetime,
        "omics_type": "Lipidomics"
    }
}
```

---

## Module Organization

### Refactored Structure

The codebase has been refactored into focused modules:

```
amprenta_rag/
├── ingestion/
│   ├── features/
│   │   ├── constants.py          # Metabolite synonyms, lists
│   │   ├── normalization.py      # Feature name normalization
│   │   ├── extraction.py         # Feature extraction from data
│   │   └── linking.py            # Feature linking to Notion
│   ├── signatures/
│   │   ├── short_id.py            # Short ID generation
│   │   ├── signature_crud.py      # Signature CRUD operations
│   │   ├── component_crud.py      # Component CRUD operations
│   │   └── species_crud.py        # Species CRUD operations
│   ├── signature_matching/
│   │   ├── models.py              # Data models
│   │   ├── species_mapping.py      # Lipid name mapping
│   │   ├── signature_loader.py    # Load signatures from Notion
│   │   ├── matching.py            # Scoring and matching logic
│   │   └── writeback.py           # Notion writeback
│   ├── metadata/
│   │   ├── helpers.py              # Shared helpers
│   │   ├── signature_metadata.py   # Signature metadata collection
│   │   ├── literature_extraction.py
│   │   ├── email_extraction.py
│   │   ├── experiment_extraction.py
│   │   └── dataset_extraction.py
│   └── ...
├── query/
│   ├── cross_omics/
│   │   ├── helpers.py              # Shared helpers
│   │   ├── prompt_templates.py     # LLM prompts
│   │   ├── synthesis.py            # LLM synthesis
│   │   ├── program_summary.py
│   │   ├── signature_summary.py
│   │   ├── feature_summary.py
│   │   └── dataset_summary.py
│   └── ...
└── ...
```

---

## Data Flow

### Ingestion Flow

```
1. User uploads CSV/TSV file
   │
2. Parse file → Extract features
   │
3. Normalize feature names
   │
4. Create/update dataset page (Notion)
   │
5. Link features to dataset
   │
6. Score against signatures
   │
7. Update dataset with matches
   │
8. Build text representation
   │
9. Chunk text → Embed → Upsert (Pinecone)
   │
10. Update dataset with embedding IDs
```

### Query Flow

```
1. User submits query
   │
2. Embed query (OpenAI)
   │
3. Search Pinecone (vector similarity)
   │
4. Retrieve top-k chunks
   │
5. Fetch full text from Notion (if needed)
   │
6. Synthesize answer with LLM
   │
7. Return answer + sources
```

---

## Notion Schema

### Databases

1. **Lipid Signatures**
   - Name, Short ID, Type, Modalities, Description
   - Relations: Components, Related Datasets

2. **Lipid Signature Components**
   - Component Name, Feature Type, Direction, Weight
   - Relations: Signature, Feature (Gene/Protein/Metabolite/Lipid Species)

3. **Experimental Data Assets**
   - Title, Omics Type, Summary, Signature Match Score
   - Relations: Related Signature(s), Programs, Experiments

4. **Feature Databases**
   - Gene Features, Protein Features, Metabolite Features, Lipid Species
   - Relations: Related Datasets

5. **Programs** (optional)
   - Program name, description
   - Relations: Experiments, Datasets

6. **Experiments** (optional)
   - Experiment name, type, disease, matrix
   - Relations: Programs, Datasets, Signatures

---

## Performance Optimizations

1. **Feature Caching**
   - In-memory cache with TTL
   - Reduces Notion API calls
   - Preloading for batch operations

2. **Batch Operations**
   - Batch signature scoring
   - Batch Pinecone upserts (100 vectors per batch)
   - Parallel feature extraction

3. **Idempotency**
   - Safe to retry operations
   - No duplicate creation
   - Graceful updates

---

## Error Handling

- **Non-blocking**: Errors logged but don't stop execution
- **Graceful degradation**: Missing data handled gracefully
- **Retry logic**: Automatic retries for transient failures
- **Clear logging**: Consistent prefixes for easy filtering

---

## Security

- **API Keys**: Stored in `.env` file (not committed)
- **Notion Permissions**: Integration token with specific database access
- **Pinecone**: Namespace-based isolation
- **OpenAI**: API key with usage limits

---

## See Also

- [API Reference](API_REFERENCE.md)
- [Usage Examples](USAGE_EXAMPLES.md)
- [Configuration Guide](CONFIGURATION.md)

