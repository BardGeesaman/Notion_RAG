# RAG Signature Similarity Query Implementation

## Overview

This document describes the implementation of the signature similarity query interface in the RAG layer, enabling queries like "Which datasets best match the ALS ceramide signature?" via a new CLI interface.

## Implementation Summary

### Files Modified

1. **`amprenta_rag/query/rag_engine.py`**
   - Added `signature_similarity_query()` function
   - Implements dataset → signature scoring and ranking
   - Reuses existing signature matching logic from `signature_matching.py`

2. **`amprenta_rag/query/rag_query_engine.py`**
   - Exported `signature_similarity_query` for backward compatibility

3. **`scripts/rag_query.py`**
   - Added `--signature-score` CLI argument
   - Made `--query` optional when using `--signature-score`
   - Added formatted output for signature similarity results

### Functions Added

#### `signature_similarity_query(dataset_page_id: str, top_k: int = 10) -> List[Dict[str, Any]]`

**Purpose**: Compute and return top matching signatures for a given dataset.

**Process**:
1. Fetches dataset page from Notion
2. Extracts page content and mwTab data (with API fallback)
3. Extracts species from mwTab metabolite data
4. Fetches all signatures from Notion
5. Scores each signature against the dataset
6. Sorts by score (descending), then by overlap_fraction
7. Returns top_k results with full metadata

**Returns**: List of dictionaries containing:
- `signature_page_id`: Notion page ID
- `signature_name`: Signature name
- `signature_short_id`: Short ID (if available)
- `disease`: Disease context (if available)
- `matrix`: Matrix type (if available)
- `score`: Total signature score
- `overlap_fraction`: Fraction of components matched
- `matched_components`: List of matched species
- `missing_components`: List of missing species
- `conflicting_components`: List of conflicting species

## CLI Usage

### Basic Usage

```bash
python scripts/rag_query.py \
  --signature-score <dataset_page_id> \
  --top-k 5
```

### Examples

**ST004217 (Ceramide Accumulation Study)**:
```bash
python scripts/rag_query.py \
  --signature-score 2beadf61-42ab-81b7-b681-cf9fefb348d0 \
  --top-k 5
```

**ST004223 (Sulfatide/Sphingolipid Study)**:
```bash
python scripts/rag_query.py \
  --signature-score 2beadf61-42ab-81a3-b07b-d60f5ee254a4 \
  --top-k 5
```

**Synthetic Test Dataset**:
```bash
python scripts/rag_query.py \
  --signature-score 2beadf61-42ab-81d0-b92a-d5ef76a3186c \
  --top-k 5
```

## Test Results

### ST004217 (Perfect Match)
- **Score**: 0.650
- **Overlap**: 1.00 (6/6 components)
- **Top Matches**:
  1. ALS-CSF-Core-6Ceramides (score: 0.650, overlap: 1.00)
  2. test_signature_verification (score: 0.650, overlap: 1.00)

### ST004223 (Partial Match)
- **Score**: 0.543
- **Overlap**: 0.33 (2/6 components)
- **Top Matches**:
  1. ALS-CSF-Core-6Ceramides (score: 0.543, overlap: 0.33)
  2. test_signature_verification (score: 0.533, overlap: 0.33)

### Synthetic Dataset
- **Score**: 0.650
- **Overlap**: 1.00 (6/6 components)
- **Top Matches**:
  1. ALS-CSF-Core-6Ceramides (score: 0.650, overlap: 1.00)
  2. test_signature_verification (score: 0.650, overlap: 1.00)

### Error Handling
- **Invalid Page ID**: Gracefully handles errors, returns empty list
- **No mwTab Data**: Logs warning, returns empty list
- **No Species Extracted**: Logs warning, returns empty list

## Output Format

```
🔎 Signature similarity for dataset <page_id>

1. <Signature Name>
   - Score: X.XXX
   - Overlap: X.XX (N/M components)
   - Matched: <species list>
   - Missing: <species list>
   - Conflicting: <species list>
   - Short ID: <short_id> (if available)
   - Disease: <disease> (if available)
   - Matrix: <matrix> (if available)

2. <Signature Name>
   ...
```

## Logging

### Info Logs
- `[RAG][SIGNATURE-SCORE] Computing signature similarity for dataset <page_id> (top_k=<n>)`
- `[RAG][SIGNATURE-SCORE] Extracted N species from dataset <page_id>`
- `[RAG][SIGNATURE-SCORE] Found N matching signatures (returning top N)`

### Debug Logs
- `[RAG][SIGNATURE-SCORE] Match: <name> (score=X.XXX, overlap=X.XX, matched=N, missing=M)`

### Warning Logs
- `[RAG][SIGNATURE-SCORE] Could not extract content from dataset page <page_id>`
- `[RAG][SIGNATURE-SCORE] No mwTab data found for dataset <page_id>`
- `[RAG][SIGNATURE-SCORE] No species extracted from dataset <page_id>`
- `[RAG][SIGNATURE-SCORE] No signatures found in Notion`
- `[RAG][SIGNATURE-SCORE] Error processing signature <id>: <error>`

### Error Logs
- `[RAG][SIGNATURE-SCORE] ERROR: Failed to compute signature similarity for dataset <page_id>: <error>`

## Implementation Details

### Species Extraction
- Reuses `extract_features_from_mwtab()` from `feature_extraction.py`
- Extracts from `MS_METABOLITE_DATA`, `GC_METABOLITE_DATA`, `LC_METABOLITE_DATA`, `METABOLITE_DATA` sections
- Maps raw lipid names to canonical species using `map_raw_lipid_to_canonical_species()`

### Signature Scoring
- Reuses `score_signature_against_dataset()` from `signature_matching.py`
- Uses same scoring algorithm as ingestion pipeline
- Calculates overlap fraction: `matched_count / total_components`

### Sorting
- Primary: By `score` (descending)
- Secondary: By `overlap_fraction` (descending)
- Returns top_k results

### Error Handling
- All errors are logged but non-fatal
- Returns empty list on failure
- Graceful degradation for missing data

## Limitations & Future Enhancements

### Current Limitations
1. **No Caching**: Signatures are fetched from Notion on every query
2. **No Dataset Scope Limiting**: Scores all signatures against dataset
3. **No Inverse Query**: "Signature → Dataset" direction not yet implemented
4. **No Score Threshold**: Returns all matches, not just above threshold

### Recommended Improvements
1. **Add Caching**: Cache signature definitions to reduce Notion API calls
2. **Add Score Threshold**: Filter results by minimum score/overlap
3. **Implement Inverse Query**: `--signature-score-from-signature <signature_page_id>`
4. **Add Dataset Filtering**: Limit to specific dataset subsets
5. **Add Performance Metrics**: Track query execution time
6. **Add Batch Queries**: Support multiple datasets in one query

## Integration Points

### Reused Components
- `signature_matching.py`: Signature loading and scoring logic
- `feature_extraction.py`: Species extraction from mwTab
- `mwtab_extraction.py`: mwTab data extraction and API fallback
- `dataset_notion_utils.py`: Dataset page fetching
- `notion_pages.py`: Page content extraction

### New Dependencies
- None - all functionality uses existing modules

## Verification

✅ **Functionality**: All core features working  
✅ **Error Handling**: Graceful degradation verified  
✅ **Logging**: Comprehensive logging in place  
✅ **CLI Interface**: User-friendly output format  
✅ **Test Coverage**: All test datasets verified  

## Summary

The signature similarity query interface is fully implemented and tested. Users can now query which signatures best match a given dataset using a simple CLI command. The implementation reuses existing signature matching logic, ensuring consistency with the ingestion pipeline.

