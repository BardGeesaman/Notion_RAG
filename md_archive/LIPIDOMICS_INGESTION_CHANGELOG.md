# Internal Lipidomics Ingestion Pipeline - Changelog

## Version: 1.0.0
**Date**: December 3, 2025  
**Status**: ✅ Implemented and Tested

## What Changed

### New Feature
- Added complete internal lipidomics ingestion pipeline for CSV/TSV files
- Supports vendor-style lipid name normalization
- Automatic signature scoring and Notion integration
- RAG embedding for internal datasets

### Files Created
- `amprenta_rag/ingestion/lipidomics_ingestion.py` - Main ingestion module
- `scripts/ingest_lipidomics.py` - CLI entry point

### Key Functions Added

#### `ingest_lipidomics_file()`
Main orchestration function that:
- Extracts species from CSV/TSV files
- Normalizes lipid names to canonical format
- Creates/updates Notion Experimental Data Asset pages
- Scores against existing signatures
- Embeds into Pinecone for RAG

#### `normalize_lipid_species()`
Extended species normalization supporting:
- Vendor formats: `CER 16:0`, `SM 24:1;O2`, `hex_cer_24_0`
- Canonical formats: `Cer(d18:1/16:0)`, `SM(d18:1/24:1)`
- Edge cases: adducts (+H, -H2O), modifications (;O2), separators (_ → /)

#### `extract_species_from_file()`
Flexible CSV/TSV parsing with:
- Automatic delimiter detection (comma or tab)
- Column name detection (species, lipid, Name, Molecule, etc.)
- Species extraction and normalization

#### `create_lipidomics_dataset_page()`
Creates new Experimental Data Asset pages with:
- Data Origin: "Internal – Amprenta"
- Dataset Source Type: "Processed table"
- Summary with file info and species count

#### `embed_lipidomics_dataset()`
RAG embedding with:
- Text representation of dataset
- Species list (truncated if >50)
- Signature match summary
- Pinecone upsert with proper metadata

## Technical Details

### Supported Input Formats
- **CSV**: Comma-separated values
- **TSV**: Tab-separated values
- **Delimiter**: Auto-detected

### Column Detection
Automatically detects lipid identity column from:
- `species`, `lipid`, `Lipid`, `Name`, `Molecule`, `compound`, `metabolite`
- Case-insensitive matching

### Species Normalization
Handles multiple vendor formats:
- Simple: `CER 16:0` → `Cer(d18:1/16:0)`
- With backbone: `Cer d18:1/16:0` → `Cer(d18:1/16:0)`
- Underscore: `hex_cer_24_0` → `HexCer(d18:1/24:0)`
- With adducts: `Cer(d18:1/16:0)+H` → `Cer(d18:1/16:0)`
- With modifications: `SM 24:1;O2` → `SM(d18:1/24:1)`
- Separator fix: `SM(d18:1_24:1)` → `SM(d18:1/24:1)`

### Notion Integration
- Creates pages in Experimental Data Assets database
- Sets Data Origin: "Internal – Amprenta"
- Sets Dataset Source Type: "Processed table"
- Updates Summary with file info
- Writes Signature Match Score (if matches found)
- Creates Related Signature(s) relations

### Signature Scoring
- Reuses existing `find_matching_signatures_for_dataset()` logic
- Uses same overlap threshold as MW datasets (default: 0.3)
- Writes highest match score to Notion
- Appends match summary to Summary field

### RAG Embedding
- Creates text representation: dataset name, origin, species list, signature matches
- Chunks using existing `chunk_text()` utility
- Embeds using existing `embed_texts()` utility
- Upserts to Pinecone with proper metadata
- Updates Embedding IDs on Notion page

## Testing

### Test Results
✅ **Canonical File**: 
- Page ID: `2beadf61-42ab-8167-bce5-d8cb15ced746`
- 3 species extracted
- 2 signature matches (score: 0.650)
- 1 vector embedded

✅ **Vendor File**: 
- Page ID: `2beadf61-42ab-816e-a097-cac73c6d3ed7`
- 3 species extracted and normalized
- 2 signature matches (score: 0.650)
- 1 vector embedded

✅ **Normalization Tests**: All 5 test cases pass

### Test Coverage
- [x] CSV file parsing
- [x] TSV file parsing
- [x] Column detection
- [x] Species normalization (vendor formats)
- [x] Notion page creation
- [x] Notion page update
- [x] Signature scoring
- [x] Signature Match Score writeback
- [x] RAG embedding
- [x] Error handling

## Breaking Changes
**None** - This is a purely additive feature.

## Migration Notes
**No migration required** - The feature works independently of existing datasets.

## Known Issues
1. **File Attachment**: File path noted in Summary; full Notion file attachment not yet implemented (requires external storage)
2. **Program/Experiment Relations**: Placeholder implementation ready; requires schema confirmation

## Next Steps
- Implement full Notion file attachment API
- Complete Program/Experiment relation updates when schema confirmed
- Add quantitative data support (intensity/abundance in scoring)
- Add batch processing for multiple files
- Add mzTab format support (Phase 2)

