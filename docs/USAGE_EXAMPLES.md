# Usage Examples

**Last Updated**: 2025-12-03

This document provides practical examples for common use cases.

---

## Table of Contents

1. [Ingesting Omics Data](#ingesting-omics-data)
2. [Ingesting Signatures](#ingesting-signatures)
3. [Querying the RAG System](#querying-the-rag-system)
4. [Cross-Omics Reasoning](#cross-omics-reasoning)
5. [Signature Matching](#signature-matching)
6. [Feature Linking](#feature-linking)

---

## Ingesting Omics Data

### Lipidomics Dataset

```python
from amprenta_rag.ingestion.lipidomics_ingestion import ingest_lipidomics_file

# Ingest a lipidomics CSV file
result = ingest_lipidomics_file(
    file_path="data/lipidomics_study.csv",
    dataset_name="ST004396",
    program_ids=["program-abc-123"],
    experiment_ids=["experiment-xyz-789"]
)

print(f"Dataset page ID: {result['dataset_page_id']}")
print(f"Features linked: {result['features_linked']}")
print(f"Signature matches: {len(result['signature_matches'])}")
```

### Metabolomics Dataset

```python
from amprenta_rag.ingestion.metabolomics_ingestion import ingest_metabolomics_file

result = ingest_metabolomics_file(
    file_path="data/metabolomics_study.tsv",
    dataset_name="MTBLS123",
    program_ids=["program-abc-123"]
)
```

### Proteomics Dataset

```python
from amprenta_rag.ingestion.proteomics_ingestion import ingest_proteomics_file

result = ingest_proteomics_file(
    file_path="data/proteomics_study.csv",
    dataset_name="PXD012345",
    experiment_ids=["experiment-xyz-789"]
)
```

### Transcriptomics Dataset

```python
from amprenta_rag.ingestion.transcriptomics_ingestion import ingest_transcriptomics_file

result = ingest_transcriptomics_file(
    file_path="data/dge_results.tsv",
    dataset_name="GSE12345",
    program_ids=["program-abc-123"]
)
```

---

## Ingesting Signatures

### From TSV File

```python
from amprenta_rag.ingestion.signature_ingestion import ingest_signature_from_file

result = ingest_signature_from_file(
    signature_file="data/signatures/ALS-CSF-Core-6Ceramides.tsv",
    signature_type="Literature-derived",
    data_ownership="Public"
)

print(f"Signature page ID: {result['signature_page_id']}")
print(f"Components created: {result['components_created']}")
print(f"Species created: {result['species_created']}")
```

### Signature File Format

**TSV/CSV with columns**:
- `feature_type` (optional, auto-detected): gene, protein, metabolite, lipid
- `feature_name`: Name of the feature
- `direction`: ↑ (up) or ↓ (down)
- `weight`: Numeric weight (default: 1.0)

**Example**:
```tsv
feature_type	feature_name	direction	weight
lipid	Cer(d18:1/16:0)	↑	1.0
lipid	Cer(d18:1/18:0)	↑	0.8
lipid	SM(d18:1/16:0)	↓	0.5
gene	TP53	↑	1.0
protein	P04637	↓	0.7
```

### Bulk Ingestion

```python
from scripts.bulk_ingest_signatures import bulk_ingest_signatures

# Ingest all signatures in SIGNATURES_DIR
results = bulk_ingest_signatures()

for result in results:
    print(f"{result['signature_name']}: {result['status']}")
```

---

## Querying the RAG System

### Basic RAG Query

```python
from amprenta_rag.query.rag_engine import query_rag

results = query_rag(
    query="What are the key ceramide signatures in ALS?",
    top_k=10,
    source_type="Signature"
)

for result in results:
    print(f"Score: {result['score']}")
    print(f"Page ID: {result['metadata']['source_page_id']}")
    print(f"Text: {result['text'][:200]}...")
    print()
```

### Filter by Source Type

```python
# Query only datasets
results = query_rag(
    query="lipidomics studies in CSF",
    top_k=20,
    source_type="Dataset"
)

# Query only experiments
results = query_rag(
    query="ALS experiments with ceramide measurements",
    top_k=15,
    source_type="Experiment"
)
```

### Signature Similarity Query

```python
from amprenta_rag.query.rag_engine import signature_similarity_query

# Find signatures matching a dataset
matches = signature_similarity_query(
    dataset_page_id="dataset-abc-123",
    top_k=10
)

for match in matches:
    print(f"{match['signature_name']}:")
    print(f"  Score: {match['score']}")
    print(f"  Overlap: {match['overlap_fraction']}")
    print(f"  Matched: {match['matched_components']}")
    print(f"  Missing: {match['missing_components']}")
```

---

## Cross-Omics Reasoning

### Program Summary

```python
from amprenta_rag.query.cross_omics import cross_omics_program_summary

summary = cross_omics_program_summary(
    program_page_id="program-abc-123",
    top_k_per_omics=20
)

print(summary)
```

**Output**: High-level summary of multi-omics evidence for the program, including:
- Per-omics findings (lipidomics, metabolomics, proteomics, transcriptomics)
- Cross-omics convergence
- Cross-omics divergence
- Disease/model system context
- Open questions

### Signature Summary

```python
from amprenta_rag.query.cross_omics import cross_omics_signature_summary

summary = cross_omics_signature_summary(
    signature_page_id="signature-xyz-789",
    top_k_datasets=20,
    top_k_chunks=100
)

print(summary)
```

**Output**: Analysis of how the signature appears across datasets, including:
- Datasets matching the signature
- Feature-level patterns
- Disease/model system associations
- Cross-omics implications

### Feature Summary

```python
from amprenta_rag.query.cross_omics import cross_omics_feature_summary

summary = cross_omics_feature_summary(
    feature_name="Cer(d18:1/16:0)",
    feature_type="lipid",
    top_k_datasets=20,
    top_k_chunks=100
)

print(summary)
```

**Output**: Comprehensive analysis of a specific feature across all omics types.

### Dataset Summary

```python
from amprenta_rag.query.cross_omics import cross_omics_dataset_summary

summary = cross_omics_dataset_summary(
    dataset_page_id="dataset-abc-123",
    top_k_chunks=100
)

print(summary)
```

**Output**: Multi-omics context for a dataset, including:
- Related experiments and programs
- Matching signatures
- Feature-level insights
- Cross-omics connections

---

## Signature Matching

### Find Matching Signatures for a Dataset

```python
from amprenta_rag.ingestion.signature_matching import find_matching_signatures_for_dataset

matches = find_matching_signatures_for_dataset(
    dataset_page_id="dataset-abc-123",
    overlap_threshold=0.3
)

for match in matches:
    print(f"{match.signature_name}:")
    print(f"  Score: {match.score}")
    print(f"  Overlap: {match.overlap_fraction}")
    print(f"  Matched components: {match.matched_components}")
    print(f"  Missing components: {match.missing_components}")
```

### Update Dataset with Matches

```python
from amprenta_rag.ingestion.signature_matching import (
    find_matching_signatures_for_dataset,
    update_dataset_with_signature_matches
)

# Find matches
matches = find_matching_signatures_for_dataset(
    dataset_page_id="dataset-abc-123",
    overlap_threshold=0.3
)

# Update Notion page
update_dataset_with_signature_matches(
    dataset_page_id="dataset-abc-123",
    matches=matches
)
```

### Batch Signature Scoring

```python
from amprenta_rag.ingestion.batch_signature_scoring import score_datasets_against_signatures_batch

dataset_ids = [
    "dataset-abc-123",
    "dataset-xyz-789",
    "dataset-def-456"
]

results = score_datasets_against_signatures_batch(
    dataset_page_ids=dataset_ids,
    overlap_threshold=0.3,
    preload_cache=True
)

for dataset_id, matches in results.items():
    print(f"{dataset_id}: {len(matches)} matches")
```

---

## Feature Linking

### Link a Single Feature

```python
from amprenta_rag.ingestion.features import link_feature

# Link a metabolite to a dataset
link_feature(
    feature_type="metabolite",
    feature_name="Glutamate",
    dataset_page_id="dataset-abc-123"
)

# Link a gene to a dataset
link_feature(
    feature_type="gene",
    feature_name="TP53",
    dataset_page_id="dataset-abc-123"
)
```

### Extract and Link Features from Text

```python
from amprenta_rag.ingestion.features import (
    extract_features_from_text,
    link_features_to_notion_items
)

# Extract features from text
text = "Increased levels of glutamate, citrate, and alanine were observed."
features = extract_features_from_text(text)

# Link all features to a dataset
link_features_to_notion_items(
    feature_names=features,
    item_page_id="dataset-abc-123",
    item_type="dataset"
)
```

---

## CLI Examples

### Ingest Dataset

```bash
python scripts/ingest_lipidomics.py \
    --file data/lipidomics_study.csv \
    --dataset-name ST004396 \
    --program-id program-abc-123
```

### Ingest Signature

```bash
python scripts/ingest_signature.py \
    --signature-file data/signatures/ALS-CSF-Core-6Ceramides.tsv \
    --signature-type "Literature-derived"
```

### RAG Query

```bash
python scripts/rag_query.py \
    --query "ceramide signature in ALS" \
    --source-type Signature \
    --top-k 10
```

### Cross-Omics Summary

```bash
python scripts/rag_query.py \
    --cross-omics-program program-abc-123

python scripts/rag_query.py \
    --cross-omics-signature signature-xyz-789

python scripts/rag_query.py \
    --cross-omics-feature "Cer(d18:1/16:0)" lipid

python scripts/rag_query.py \
    --cross-omics-dataset dataset-abc-123
```

### Signature Similarity

```bash
python scripts/rag_query.py \
    --signature-score dataset-abc-123 \
    --top-k 10
```

---

## Error Handling

All functions are designed to be non-blocking and idempotent:

```python
try:
    result = ingest_lipidomics_file(...)
except Exception as e:
    # Errors are logged but don't stop execution
    logger.warning(f"Ingestion failed: {e}")
    # Continue with other operations
```

---

## Performance Tips

1. **Use Feature Cache**: Enable caching for batch operations
   ```python
   from amprenta_rag.ingestion.dataset_feature_cache import get_feature_cache
   cache = get_feature_cache(ttl_seconds=3600)
   ```

2. **Batch Operations**: Use batch functions for multiple datasets
   ```python
   score_datasets_against_signatures_batch(
       dataset_page_ids=[...],
       preload_cache=True
   )
   ```

3. **Parallel Processing**: Process multiple files in parallel (if needed)
   ```python
   from concurrent.futures import ThreadPoolExecutor
   with ThreadPoolExecutor(max_workers=4) as executor:
       results = executor.map(ingest_lipidomics_file, file_list)
   ```

---

## See Also

- [API Reference](API_REFERENCE.md)
- [Architecture Overview](ARCHITECTURE.md)
- [Configuration Guide](CONFIGURATION.md)

