# API Reference

**Last Updated**: 2025-12-03

This document provides a comprehensive API reference for the Amprenta RAG multi-omics platform.

---

## Table of Contents

1. [REST API Endpoints](#rest-api-endpoints)
2. [Ingestion Modules](#ingestion-modules)
3. [Query Modules](#query-modules)
4. [Signature Modules](#signature-modules)
5. [Feature Modules](#feature-modules)
6. [Metadata Modules](#metadata-modules)
7. [Client Modules](#client-modules)

---

## REST API Endpoints

The Amprenta platform exposes a **FastAPI REST API** for programmatic access to multi-omics data, chemistry compounds, and HTS screening campaigns.

### Base URL

```
http://localhost:8000
```

### API Version

All endpoints are prefixed with `/api/v1/`

---

### Compounds API

#### List All Compounds

```http
GET /api/v1/compounds/
```

**Description**: Retrieve all compounds from the chemistry database.

**Response** (200 OK):
```json
[
  {
    "compound_id": "AMPR-001",
    "smiles": "CC(C)Cc1ccc(cc1)C(C)C(=O)O",
    "inchi_key": "HEFNNWSXXWATRW-UHFFFAOYSA-N",
    "canonical_smiles": "CC(C)Cc1ccc(cc1)C(C)C(=O)O",
    "molecular_formula": "C13H18O2",
    "molecular_weight": 206.28,
    "logp": 3.5,
    "hbd_count": 1,
    "hba_count": 2,
    "rotatable_bonds": 4
  },
  ...
]
```

**Example**:
```bash
curl http://localhost:8000/api/v1/compounds/
```

---

#### Get Compound by ID

```http
GET /api/v1/compounds/{compound_id}
```

**Description**: Retrieve a specific compound by its ID.

**Parameters**:
- `compound_id` (path, required): Compound identifier (e.g., "AMPR-001")

**Response** (200 OK):
```json
{
  "compound_id": "AMPR-001",
  "smiles": "CC(C)Cc1ccc(cc1)C(C)C(=O)O",
  "inchi_key": "HEFNNWSXXWATRW-UHFFFAOYSA-N",
  "canonical_smiles": "CC(C)Cc1ccc(cc1)C(C)C(=O)O",
  "molecular_formula": "C13H18O2",
  "molecular_weight": 206.28,
  "logp": 3.5,
  "hbd_count": 1,
  "hba_count": 2,
  "rotatable_bonds": 4
}
```

**Response** (404 Not Found):
```json
{
  "detail": "Compound not found"
}
```

**Example**:
```bash
curl http://localhost:8000/api/v1/compounds/AMPR-001
```

---

#### Get Programs Linked to Compound

```http
GET /api/v1/compounds/{compound_id}/programs
```

**Description**: Retrieve all programs linked to a compound.

**Parameters**:
- `compound_id` (path, required): Compound identifier

**Response** (200 OK):
```json
[
  {
    "program_id": "2b9adf6142ab8026853ef58f725665a6",
    "status": "Active",
    "notes": "Lead compound for ALS program",
    "created_at": "2025-01-15T10:30:00Z"
  },
  ...
]
```

**Example**:
```bash
curl http://localhost:8000/api/v1/compounds/AMPR-001/programs
```

---

### Screening API

#### List HTS Campaigns

```http
GET /api/v1/screening/campaigns
```

**Description**: List all high-throughput screening campaigns.

**Response** (200 OK):
```json
[
  {
    "campaign_id": "HTS-2024-001",
    "campaign_name": "Cer Synthesis Inhibitor Screen",
    "description": "384-well primary screen for ceramide synthesis inhibitors",
    "assay_type": "Enzymatic",
    "target": "CerS1",
    "library_id": "DIVERSet-CL",
    "total_wells": 5120,
    "hit_count": 47,
    "run_date": "2024-03-15"
  },
  ...
]
```

**Example**:
```bash
curl http://localhost:8000/api/v1/screening/campaigns
```

---

#### Get HTS Campaign by ID

```http
GET /api/v1/screening/campaigns/{campaign_id}
```

**Description**: Retrieve details for a specific HTS campaign.

**Parameters**:
- `campaign_id` (path, required): Campaign identifier

**Response** (200 OK):
```json
{
  "campaign_id": "HTS-2024-001",
  "campaign_name": "Cer Synthesis Inhibitor Screen",
  "description": "384-well primary screen for ceramide synthesis inhibitors",
  "assay_type": "Enzymatic",
  "target": "CerS1",
  "library_id": "DIVERSet-CL",
  "total_wells": 5120,
  "hit_count": 47,
  "run_date": "2024-03-15"
}
```

**Response** (404 Not Found):
```json
{
  "detail": "Campaign not found"
}
```

**Example**:
```bash
curl http://localhost:8000/api/v1/screening/campaigns/HTS-2024-001
```

---

#### Get HTS Hits for Campaign

```http
GET /api/v1/screening/campaigns/{campaign_id}/hits
```

**Description**: Retrieve all hits (positive results) for a screening campaign.

**Parameters**:
- `campaign_id` (path, required): Campaign identifier

**Response** (200 OK):
```json
[
  {
    "result_id": "RES-001",
    "compound_id": "AMPR-123",
    "well_position": "A01",
    "raw_value": 1250.5,
    "normalized_value": 85.3,
    "z_score": 3.2,
    "hit_flag": 1,
    "hit_category": "Strong"
  },
  {
    "result_id": "RES-002",
    "compound_id": "AMPR-456",
    "well_position": "B03",
    "raw_value": 980.2,
    "normalized_value": 72.1,
    "z_score": 2.5,
    "hit_flag": 1,
    "hit_category": "Moderate"
  },
  ...
]
```

**Example**:
```bash
curl http://localhost:8000/api/v1/screening/campaigns/HTS-2024-001/hits
```

---

### Other API Endpoints

The platform also exposes REST APIs for:

- **Programs** (`/api/v1/programs/`)
- **Experiments** (`/api/v1/experiments/`)
- **Datasets** (`/api/v1/datasets/`)
- **Features** (`/api/v1/features/`)
- **Signatures** (`/api/v1/signatures/`)

See FastAPI interactive documentation at `/docs` when server is running.

---

### Running the API Server

#### Development Server

```bash
# Start FastAPI development server
cd /path/to/RAG
uvicorn amprenta_rag.api.main:app --reload --port 8000

# API available at: http://localhost:8000
# Interactive docs: http://localhost:8000/docs
# OpenAPI spec: http://localhost:8000/openapi.json
```

#### Production Server

```bash
# Use gunicorn with uvicorn workers
gunicorn amprenta_rag.api.main:app \
    --workers 4 \
    --worker-class uvicorn.workers.UvicornWorker \
    --bind 0.0.0.0:8000
```

---

### API Error Handling

All endpoints follow consistent error handling:

#### HTTP Status Codes

- **200 OK**: Successful request
- **201 Created**: Resource created successfully
- **204 No Content**: Resource deleted successfully
- **400 Bad Request**: Invalid request parameters
- **404 Not Found**: Resource not found
- **500 Internal Server Error**: Server error (check logs)

#### Error Response Format

```json
{
  "detail": "Human-readable error message"
}
```

**Example Error**:
```bash
curl http://localhost:8000/api/v1/compounds/INVALID-ID

# Response (404):
{
  "detail": "Compound not found"
}
```

---

### API Authentication

**Current**: No authentication (development mode)

**Future**: OAuth2/JWT authentication for production deployments

---

## Ingestion Modules

### Feature Extraction and Linking

#### `amprenta_rag.ingestion.features`

**Purpose**: Extract, normalize, and link omics features (genes, proteins, metabolites, lipids) to Notion.

**Key Functions**:

- `normalize_metabolite_name(raw: str) -> str`
  - Normalizes metabolite names from various formats (mwTab, CSV, vendor formats)
  - Handles database prefixes (HMDB, KEGG, CHEBI), synonyms, and case normalization
  - Returns canonical metabolite name

- `extract_features_from_mwtab(mwtab_json: Dict[str, Any]) -> List[str]`
  - Extracts metabolite features from Metabolomics Workbench mwTab JSON
  - Returns list of normalized metabolite names

- `extract_features_from_text(text: str) -> List[str]`
  - Scans text content for metabolite mentions (amino acids, nucleotides, common metabolites)
  - Returns list of detected metabolite names

- `link_feature(feature_type: str, feature_name: str, dataset_page_id: str) -> None`
  - Creates or finds a feature page in Notion (Gene, Protein, Metabolite, or Lipid Species DB)
  - Links the feature to the dataset page via relations
  - Idempotent: safe to call multiple times

**Example**:
```python
from amprenta_rag.ingestion.features import normalize_metabolite_name, link_feature

# Normalize a metabolite name
normalized = normalize_metabolite_name("HMDB:12345 L-glutamate")
# Returns: "Glutamate"

# Link a feature to a dataset
link_feature("metabolite", "Glutamate", "dataset-page-id-123")
```

---

### Signature Ingestion

#### `amprenta_rag.ingestion.signatures`

**Purpose**: Create and manage multi-omics signatures in Notion.

**Key Functions**:

- `generate_signature_short_id(signature_name: str, version: Optional[str] = None) -> str`
  - Generates a deterministic short ID for a signature
  - Handles special characters, truncation, and versioning
  - Returns normalized short ID (max 50 characters)

- `find_or_create_signature_page(signature: Signature, ...) -> Optional[str]`
  - Creates or finds a signature page in Notion "Lipid Signatures" database
  - Handles idempotency via name/Short ID lookup
  - Returns Notion page ID

- `find_or_create_component_page(component: SignatureComponent, ...) -> Optional[str]`
  - Creates or finds a signature component page
  - Links component to signature and feature pages
  - Returns Notion page ID

- `find_or_create_lipid_species_page(lipid_name: str) -> Optional[str]`
  - Creates or finds a lipid species page in Notion
  - Handles species normalization and synonym management
  - Returns Notion page ID

**Example**:
```python
from amprenta_rag.ingestion.signatures import find_or_create_signature_page
from amprenta_rag.signatures.signature_loader import Signature

signature = Signature(
    name="ALS-CSF-Core-6Ceramides",
    components=[...]
)

page_id = find_or_create_signature_page(signature)
```

---

### Signature Matching and Scoring

#### `amprenta_rag.ingestion.signature_matching`

**Purpose**: Score datasets against signatures and find matches.

**Key Functions**:

- `map_raw_lipid_to_canonical_species(raw_name: str) -> Optional[str]`
  - Maps vendor lipid names (e.g., "CER 16:0") to canonical format (e.g., "Cer(d18:1/16:0)")
  - Returns canonical species name or None if invalid

- `score_signature_against_dataset(signature: Signature, dataset_species: Set[str], ...) -> SignatureScoreResult`
  - Computes overlap score between a signature and dataset
  - Considers direction consistency and component weights
  - Returns detailed scoring result

- `find_matching_signatures_for_dataset(dataset_page_id: str, ...) -> List[SignatureMatchResult]`
  - Finds all signatures matching a dataset above threshold
  - Supports multi-omics scoring (genes, proteins, metabolites, lipids)
  - Returns ranked list of matches

- `update_dataset_with_signature_matches(dataset_page_id: str, matches: List[SignatureMatchResult]) -> None`
  - Updates Notion dataset page with signature matches
  - Writes relations, scores, and summary text
  - Idempotent and non-destructive

**Example**:
```python
from amprenta_rag.ingestion.signature_matching import find_matching_signatures_for_dataset

matches = find_matching_signatures_for_dataset(
    dataset_page_id="abc-123-def",
    overlap_threshold=0.3
)

for match in matches:
    print(f"{match.signature_name}: score={match.score}, overlap={match.overlap_fraction}")
```

---

### Multi-Omics Scoring

#### `amprenta_rag.ingestion.multi_omics_scoring`

**Purpose**: Extract and score multi-omics features from datasets.

**Key Functions**:

- `extract_dataset_features_by_type(dataset_page_id: str, ...) -> Dict[str, Set[str]]`
  - Extracts features from a dataset grouped by type (gene, protein, metabolite, lipid)
  - Queries Notion feature databases for linked features
  - Uses feature cache for performance
  - Returns dictionary mapping feature_type → set of feature names

- `score_multi_omics_signature_against_dataset(signature: Signature, dataset_page_id: str, ...) -> SignatureScoreResult`
  - Scores a multi-omics signature against a dataset
  - Handles mixed feature types (genes, proteins, metabolites, lipids)
  - Returns detailed scoring result

**Example**:
```python
from amprenta_rag.ingestion.multi_omics_scoring import extract_dataset_features_by_type

features_by_type = extract_dataset_features_by_type("dataset-page-id-123")

print(f"Genes: {len(features_by_type['gene'])}")
print(f"Proteins: {len(features_by_type['protein'])}")
print(f"Metabolites: {len(features_by_type['metabolite'])}")
print(f"Lipids: {len(features_by_type['lipid'])}")
```

---

### Dataset Feature Cache

#### `amprenta_rag.ingestion.dataset_feature_cache`

**Purpose**: In-memory caching of dataset features to improve performance.

**Key Functions**:

- `get_feature_cache(ttl_seconds: int = 3600) -> DatasetFeatureCache`
  - Returns singleton feature cache instance
  - Default TTL: 1 hour

- `DatasetFeatureCache.get_features(dataset_page_id: str, ...) -> Optional[Dict[str, Set[str]]]`
  - Retrieves cached features for a dataset
  - Returns None if not cached or stale

- `DatasetFeatureCache.set_features(dataset_page_id: str, features_by_type: Dict[str, Set[str]], ...) -> None`
  - Stores features in cache with timestamp
  - Automatically expires after TTL

- `DatasetFeatureCache.preload_datasets(dataset_page_ids: List[str], extract_fn: Callable, ...) -> Dict[str, bool]`
  - Preloads features for multiple datasets
  - Useful for batch scoring operations

**Example**:
```python
from amprenta_rag.ingestion.dataset_feature_cache import get_feature_cache

cache = get_feature_cache()
features = cache.get_features("dataset-page-id-123")

if features is None:
    # Cache miss - extract from Notion
    features = extract_dataset_features_by_type("dataset-page-id-123")
    cache.set_features("dataset-page-id-123", features)
```

---

## Query Modules

### Cross-Omics Reasoning

#### `amprenta_rag.query.cross_omics`

**Purpose**: Generate high-level multi-omics summaries using LLM reasoning.

**Key Functions**:

- `cross_omics_program_summary(program_page_id: str, top_k_per_omics: int = 20) -> str`
  - Summarizes multi-omics evidence for a Program
  - Retrieves related experiments, datasets, and signatures
  - Returns LLM-generated summary

- `cross_omics_signature_summary(signature_page_id: str, top_k_datasets: int = 20, top_k_chunks: int = 100) -> str`
  - Summarizes cross-omics patterns for a Signature
  - Analyzes datasets matching the signature
  - Returns LLM-generated summary

- `cross_omics_feature_summary(feature_name: str, feature_type: str, top_k_datasets: int = 20, top_k_chunks: int = 100) -> str`
  - Summarizes all omics evidence related to a specific Feature
  - Cross-references across datasets and signatures
  - Returns LLM-generated summary

- `cross_omics_dataset_summary(dataset_page_id: str, top_k_chunks: int = 100) -> str`
  - Summarizes multi-omics evidence for a Dataset
  - Includes related experiments, signatures, and features
  - Returns LLM-generated summary

**Example**:
```python
from amprenta_rag.query.cross_omics import cross_omics_program_summary

summary = cross_omics_program_summary(
    program_page_id="program-abc-123",
    top_k_per_omics=20
)
print(summary)
```

---

### RAG Engine

#### `amprenta_rag.query.rag_engine`

**Purpose**: Orchestrate RAG queries with Pinecone and OpenAI.

**Key Functions**:

- `query_rag(query: str, top_k: int = 10, source_type: Optional[str] = None, ...) -> List[Dict[str, Any]]`
  - Performs semantic search in Pinecone
  - Retrieves relevant chunks with metadata
  - Returns ranked list of chunks

- `signature_similarity_query(dataset_page_id: str, top_k: int = 10) -> List[Dict[str, Any]]`
  - Finds top-k signatures matching a dataset
  - Uses signature scoring engine
  - Returns ranked list with scores and overlap metrics

**Example**:
```python
from amprenta_rag.query.rag_engine import query_rag

results = query_rag(
    query="ceramide signature in ALS",
    top_k=10,
    source_type="Signature"
)

for result in results:
    print(f"{result['metadata']['source_page_id']}: {result['score']}")
```

---

## Metadata Modules

### Semantic Metadata Extraction

#### `amprenta_rag.ingestion.metadata`

**Purpose**: Extract structured metadata from Notion pages for RAG embedding.

**Key Functions**:

- `get_literature_semantic_metadata(parent_page_id: str, item: ZoteroItem) -> Dict[str, Any]`
  - Extracts metadata from literature Notion pages
  - Includes diseases, model systems, matrices, signatures
  - Returns structured metadata dictionary

- `get_email_semantic_metadata(email_page: Dict[str, Any]) -> Dict[str, Any]`
  - Extracts metadata from email Notion pages
  - Includes related datasets, experiments, signatures
  - Returns structured metadata dictionary

- `get_experiment_semantic_metadata(exp_page: Dict[str, Any]) -> Dict[str, Any]`
  - Extracts metadata from experiment ELN pages
  - Includes disease, matrix, model systems, readout signatures
  - Returns structured metadata dictionary

- `get_dataset_semantic_metadata(dataset_page: Dict[str, Any]) -> Dict[str, Any]`
  - Extracts metadata from dataset pages
  - Includes omics type, programs, experiments, signature matches
  - Returns structured metadata dictionary

**Example**:
```python
from amprenta_rag.ingestion.metadata import get_dataset_semantic_metadata
from amprenta_rag.ingestion.metadata.helpers import _fetch_notion_page

page = _fetch_notion_page("dataset-page-id-123")
metadata = get_dataset_semantic_metadata(page)
print(metadata["omics_type"])
print(metadata["related_programs"])
```

---

## Client Modules

### Notion Client

#### `amprenta_rag.clients.notion_client`

**Purpose**: Low-level Notion API interactions.

**Key Functions**:

- `notion_headers() -> Dict[str, str]`
  - Returns Notion API headers with authentication

- `get_page_text(page_id: str) -> str`
  - Fetches full text content from a Notion page
  - Returns plain text representation

---

### OpenAI Client

#### `amprenta_rag.clients.openai_client`

**Purpose**: OpenAI API interactions for embeddings and chat.

**Key Functions**:

- `get_openai_client() -> OpenAI`
  - Returns configured OpenAI client instance

- `get_default_models() -> Tuple[str, str]`
  - Returns (chat_model, embedding_model) tuple
  - Defaults: gpt-4o-mini, text-embedding-3-small

---

## Configuration

### `amprenta_rag.config`

**Purpose**: Centralized configuration management.

**Key Functions**:

- `get_config() -> PipelineConfig`
  - Returns global configuration object
  - Loads from environment variables (.env file)

**Configuration Options**:

- `NOTION_API_KEY`: Notion integration token
- `NOTION_SIGNATURE_DB_ID`: Lipid Signatures database ID
- `NOTION_DATASET_DB_ID`: Experimental Data Assets database ID
- `PINECONE_API_KEY`: Pinecone API key
- `PINECONE_INDEX_NAME`: Pinecone index name
- `OPENAI_API_KEY`: OpenAI API key
- `SIGNATURE_OVERLAP_THRESHOLD`: Default overlap threshold (0.3)
- `ENABLE_SIGNATURE_SCORING`: Enable signature scoring (True)
- `ENABLE_LIPID_MAPPING`: Enable lipid mapping (True)

---

## Error Handling

All ingestion and query functions use consistent error handling:

- **Non-blocking**: Errors are logged but don't stop execution
- **Idempotent**: Safe to retry operations
- **Graceful degradation**: Missing data is handled gracefully

**Logging Prefixes**:

- `[INGEST][LIPIDOMICS]`: Lipidomics ingestion
- `[INGEST][METABOLOMICS]`: Metabolomics ingestion
- `[INGEST][PROTEOMICS]`: Proteomics ingestion
- `[INGEST][TRANSCRIPTOMICS]`: Transcriptomics ingestion
- `[INGEST][SIGNATURE-MATCH]`: Signature matching
- `[INGEST][FEATURE]`: Feature linking
- `[RAG][CROSS-OMICS]`: Cross-omics reasoning
- `[RAG][SIGNATURE-SCORE]`: Signature similarity queries

---

## See Also

- [Architecture Overview](ARCHITECTURE.md)
- [Usage Examples](USAGE_EXAMPLES.md)
- [Configuration Guide](CONFIGURATION.md)

