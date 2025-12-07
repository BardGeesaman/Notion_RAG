# Internal Lipidomics Ingestion Pipeline Implementation

## Overview

This document describes the implementation of the internal Amprenta lipidomics ingestion pipeline, enabling ingestion of internal CSV/TSV lipidomics files with species normalization, signature scoring, Notion integration, and RAG embedding.

## Implementation Summary

### Files Created

1. **`amprenta_rag/ingestion/lipidomics_ingestion.py`** (600+ lines)
   - Main ingestion module with all core functionality
   - Species normalization with vendor format support
   - Notion page creation/update
   - Signature scoring integration
   - RAG embedding

2. **`scripts/ingest_lipidomics.py`** (80+ lines)
   - CLI entry point for lipidomics ingestion
   - Argument parsing and validation
   - User-friendly output

### Key Functions

#### `ingest_lipidomics_file()`
Main orchestration function that:
- Extracts species from CSV/TSV files
- Normalizes lipid names to canonical format
- Creates/updates Notion pages
- Scores against signatures
- Embeds into Pinecone

#### `normalize_lipid_species()`
Extended species normalization supporting:
- Vendor formats: `CER 16:0`, `SM 24:1;O2`, `hex_cer_24_0`
- Canonical formats: `Cer(d18:1/16:0)`, `SM(d18:1/24:1)`
- Edge cases: adducts, separators, modifications

#### `extract_species_from_file()`
Flexible CSV/TSV parsing with:
- Automatic delimiter detection
- Column name detection (species, lipid, Name, etc.)
- Species extraction and normalization

#### `create_lipidomics_dataset_page()`
Creates new Experimental Data Asset pages with:
- Data Origin: "Internal – Amprenta"
- Dataset Source Type: "Processed table"
- Summary with file info and species count

#### `embed_lipidomics_dataset()`
RAG embedding with:
- Text representation of dataset
- Species list
- Signature match summary
- Pinecone upsert with proper metadata

## Test Results

### Test File 1: Canonical Species

**File**: `test_data/test_canonical.csv`
```csv
species,intensity
Cer(d18:1/16:0),1234
SM(d18:1/24:1),5678
HexCer(d18:1/22:0),2345
```

**Results**:
- ✅ **Page Created**: `2beadf61-42ab-8167-bce5-d8cb15ced746`
- ✅ **Species Extracted**: 3 unique species
- ✅ **Normalization**: All already canonical (no changes needed)
- ✅ **Signature Matches**: 2 signatures (ALS-CSF-Core-6Ceramides, test_signature_verification)
- ✅ **Signature Match Score**: 0.650
- ✅ **Embedding**: 1 vector created and upserted to Pinecone

**Logs**:
```
[INGEST][LIPIDOMICS] Starting ingestion of lipidomics file: test_data/test_canonical.csv
[INGEST][LIPIDOMICS] Using column 'species' for lipid identity
[INGEST][LIPIDOMICS] Extracted 3 unique species from 3 rows
[INGEST][LIPIDOMICS] Created new dataset page 2beadf61-42ab-8167-bce5-d8cb15ced746
[INGEST][SIGNATURE-MATCH] Found 2 matching signature(s)
[INGEST][SIGNATURE-MATCH] Writing Signature Match Score = 0.650
[INGEST][LIPIDOMICS] Embedded dataset to Pinecone (1 vectors)
```

### Test File 2: Vendor-Style Names

**File**: `test_data/test_vendor.csv`
```csv
Lipid,Abundance
CER 16:0,1234
SM 24:1;O2,5678
hex_cer_24_0,6789
```

**Results**:
- ✅ **Page Created**: `2beadf61-42ab-816e-a097-cac73c6d3ed7`
- ✅ **Species Extracted**: 3 unique species
- ✅ **Normalization**: All vendor formats successfully normalized:
  - `CER 16:0` → `Cer(d18:1/16:0)`
  - `SM 24:1;O2` → `SM(d18:1/24:1)`
  - `hex_cer_24_0` → `HexCer(d18:1/24:0)`
- ✅ **Signature Matches**: 2 signatures (same as test 1)
- ✅ **Signature Match Score**: 0.650
- ✅ **Embedding**: 1 vector created and upserted

**Logs**:
```
[INGEST][LIPIDOMICS] Normalized raw lipid 'CER 16:0' -> 'Cer(d18:1/16:0)'
[INGEST][LIPIDOMICS] Normalized raw lipid 'SM 24:1;O2' -> 'SM(d18:1/24:1)'
[INGEST][LIPIDOMICS] Normalized raw lipid 'hex_cer_24_0' -> 'HexCer(d18:1/24:0)'
[INGEST][LIPIDOMICS] Extracted 3 unique species from 3 rows
```

### Species Normalization Tests

All normalization test cases pass:

| Input | Output | Status |
|-------|--------|--------|
| `CER 16:0` | `Cer(d18:1/16:0)` | ✅ |
| `SM 24:1;O2` | `SM(d18:1/24:1)` | ✅ |
| `hex_cer_24_0` | `HexCer(d18:1/24:0)` | ✅ |
| `Cer(d18:1/16:0)+H` | `Cer(d18:1/16:0)` | ✅ |
| `SM(d18:1_24:1)` | `SM(d18:1/24:1)` | ✅ |

## CLI Usage

### Basic Usage

**Create new page**:
```bash
python scripts/ingest_lipidomics.py \
  --file test_data/test_canonical.csv \
  --create-page
```

**Update existing page**:
```bash
python scripts/ingest_lipidomics.py \
  --file test_data/test_vendor.csv \
  --dataset-page-id 2beadf61-42ab-8167-bce5-d8cb15ced746
```

**With Program/Experiment links** (when schema is confirmed):
```bash
python scripts/ingest_lipidomics.py \
  --file test_data/test_canonical.csv \
  --create-page \
  --program-id <program_page_id> \
  --experiment-id <experiment_page_id>
```

## Features Implemented

### ✅ 1. CSV/TSV Support
- Automatic delimiter detection (comma or tab)
- Flexible column name detection
- Handles missing values gracefully

### ✅ 2. Species Normalization
- Vendor format support (CER 16:0, SM 24:1, etc.)
- Underscore format support (hex_cer_24_0)
- Adduct removal (+H, -H2O, etc.)
- Modification handling (;O2, etc.)
- Separator normalization (_ → /)

### ✅ 3. Notion Integration
- Creates new Experimental Data Asset pages
- Updates existing pages
- Sets Data Origin: "Internal – Amprenta"
- Sets Dataset Source Type: "Processed table"
- Adds file info to Summary

### ✅ 4. Signature Scoring
- Reuses existing signature matching logic
- Scores all signatures against dataset
- Writes Signature Match Score
- Updates Related Signature(s) relation
- Appends match summary to Summary field

### ✅ 5. RAG Embedding
- Creates text representation of dataset
- Includes species list and signature matches
- Chunks and embeds using existing pipeline
- Upserts to Pinecone with proper metadata
- Updates Embedding IDs on Notion page

### ✅ 6. Error Handling
- Non-fatal errors where possible
- Comprehensive logging
- Graceful degradation
- Clear error messages

## Notion Page Properties

### Created Pages Include:
- **Experiment Name**: "Internal Lipidomics — <filename>"
- **Data Origin**: "Internal – Amprenta"
- **Dataset Source Type**: "Processed table"
- **Summary**: File path, row count, species count
- **Signature Match Score**: Highest match score (if matches found)
- **Related Signature(s)**: Relations to matching signatures
- **Embedding IDs**: List of Pinecone vector IDs
- **Last Embedded**: Timestamp

## Logging

### Log Prefixes
- `[INGEST][LIPIDOMICS]` - Internal lipidomics ingestion
- `[INGEST][SIGNATURE-MATCH]` - Signature matching (reused)
- `[INGEST][DATASET]` - Dataset operations (reused)

### Key Log Messages
- File ingestion start/completion
- Column detection
- Species normalization (success/warning)
- Species extraction summary
- Page creation/update
- Signature matching results
- Embedding completion

## Limitations & Future Enhancements

### Current Limitations
1. **File Attachment**: File path noted in Summary; full Notion file attachment not yet implemented
2. **Program/Experiment Relations**: Schema-dependent; placeholder implementation ready
3. **Quantitative Data**: Intensity/abundance columns detected but not yet used in scoring
4. **Batch Processing**: Single file at a time

### Recommended Enhancements
1. **Full File Attachment**: Implement Notion file upload API
2. **Program/Experiment Relations**: Complete when schema is confirmed
3. **Quantitative Scoring**: Use intensity/abundance in signature scoring
4. **Batch Processing**: Support multiple files in one run
5. **mzTab Support**: Add mzTab format support (Phase 2)
6. **Validation**: Add file format validation and schema checking

## Integration Points

### Reused Components
- `signature_matching.py`: Signature scoring and matching
- `dataset_notion_utils.py`: Dataset page operations
- `text_embedding_utils.py`: Text chunking and embedding
- `pinecone_utils.py`: Metadata sanitization
- `notion_pages.py`: Page content extraction

### New Dependencies
- `pandas`: CSV/TSV parsing (already in environment)

## Verification Checklist

✅ **File Parsing**: CSV/TSV files parsed correctly  
✅ **Column Detection**: Lipid identity column detected automatically  
✅ **Species Normalization**: Vendor formats normalized correctly  
✅ **Page Creation**: New Notion pages created with correct properties  
✅ **Page Update**: Existing pages updated correctly  
✅ **Signature Scoring**: Signatures matched and scored correctly  
✅ **Notion Writeback**: Signature Match Score and relations written  
✅ **RAG Embedding**: Datasets embedded into Pinecone  
✅ **Error Handling**: Graceful error handling verified  
✅ **Logging**: Comprehensive logging in place  

## Summary

The internal lipidomics ingestion pipeline is fully implemented and tested. It successfully:
- Parses CSV/TSV files with flexible column detection
- Normalizes vendor-style lipid names to canonical format
- Creates/updates Notion Experimental Data Asset pages
- Scores datasets against existing signatures
- Embeds datasets into Pinecone for RAG queries

All core requirements are met and the system is ready for production use with internal lipidomics datasets.

