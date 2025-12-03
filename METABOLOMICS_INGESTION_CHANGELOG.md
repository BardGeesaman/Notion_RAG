# Internal Metabolomics Ingestion Pipeline - Changelog

## Version: 1.0.0
**Date**: December 3, 2025  
**Status**: ✅ Implemented and Tested

## What Changed

### New Feature
- Added complete internal polar metabolomics ingestion pipeline for CSV/TSV files
- Supports various metabolite name formats with normalization
- Automatic Notion integration and RAG embedding
- Program/Experiment linking support

### Files Created
- `amprenta_rag/ingestion/metabolomics_ingestion.py` - Main ingestion module
- `scripts/ingest_metabolomics.py` - CLI entry point

### Key Functions Added

#### `ingest_metabolomics_file()`
Main orchestration function that:
- Extracts metabolites from CSV/TSV files
- Normalizes metabolite names to canonical format
- Creates/updates Notion Experimental Data Asset pages
- Links to Programs and Experiments
- Embeds into Pinecone for RAG

#### `normalize_metabolite_name()`
Metabolite name normalization supporting:
- Adduct removal: `[M+H]`, `[M-H]`, `[M+Na]`, etc.
- Annotation removal: `(pos)`, `(neg)`, `(+)`, `(-)`
- Underscore to space conversion
- Synonym mapping (20+ common metabolites)
- Case normalization (Title Case)

#### `extract_metabolite_set_from_file()`
Flexible CSV/TSV parsing with:
- Automatic delimiter detection (comma or tab)
- Column name detection (metabolite, compound, Name, Molecule, etc.)
- Metabolite extraction and normalization

#### `create_metabolomics_dataset_page()`
Creates new Experimental Data Asset pages with:
- Data Origin: "Internal – Amprenta"
- Dataset Source Type: "Processed table"
- Optional: Omics Type = Metabolomics (gracefully handles missing property)
- Summary with file info and metabolite count

#### `embed_metabolomics_dataset()`
RAG embedding with:
- Text representation of dataset
- Metabolite list (truncated if >100)
- Program/Experiment link counts
- Pinecone upsert with proper metadata (omics_type = "Metabolomics")
- Updates Embedding IDs and Last Embedded on Notion page

## Technical Details

### Supported Input Formats
- **CSV**: Comma-separated values
- **TSV**: Tab-separated values
- **Delimiter**: Auto-detected

### Column Detection
Automatically detects metabolite identity column from:
- `metabolite`, `Metabolite`, `compound`, `Compound`, `name`, `Name`, `molecule`, `Molecule`
- Case-insensitive matching

### Metabolite Normalization
Handles multiple formats:
- Adducts: `Serine [M+H]+` → `Serine`
- Annotations: `Glutamine (pos)` → `Glutamine`
- Case variations: `glutamine` → `Glutamine`, `GLUTAMIC ACID` → `Glutamate`
- Synonyms: `L-Glutamic acid` → `Glutamate`, `L-Glutamine` → `Glutamine`
- Underscores: `L_Serine` → `L Serine` → `Serine`

### Notion Integration
- Creates pages in Experimental Data Assets database
- Sets Data Origin: "Internal – Amprenta"
- Sets Dataset Source Type: "Processed table"
- Attempts to set Omics Type: "Metabolomics" (gracefully handles missing property)
- Updates Summary with file info
- Writes Embedding IDs and Last Embedded

### RAG Embedding
- Creates text representation: dataset name, origin, metabolite list, Program/Experiment links
- Chunks using existing `chunk_text()` utility
- Embeds using existing `embed_texts()` utility
- Upserts to Pinecone with proper metadata:
  - `source_type = "Dataset"`
  - `omics_type = "Metabolomics"`
- Updates Embedding IDs on Notion page

## Testing

### Test Results
✅ **Canonical File**: 
- Page ID: `2beadf61-42ab-81bf-8544-fd2c918be46f`
- 3 metabolites extracted and normalized
- 1 vector embedded

✅ **Mixed Format File**: 
- Page ID: `2beadf61-42ab-8171-9600-dd0aaa1eb7a2`
- 3 metabolites extracted and normalized
- 1 vector embedded

✅ **Normalization Tests**: All 6 test cases pass

### Test Coverage
- [x] CSV file parsing
- [x] TSV file parsing
- [x] Column detection
- [x] Metabolite normalization (various formats)
- [x] Notion page creation
- [x] Notion page update
- [x] RAG embedding
- [x] Error handling
- [x] Omics Type property handling (graceful degradation)

## Breaking Changes
**None** - This is a purely additive feature.

## Migration Notes
**No migration required** - The feature works independently of existing datasets.

## Known Issues
1. **File Attachment**: File path noted in Summary; full Notion file attachment not yet implemented (requires external storage)
2. **Program/Experiment Relations**: Placeholder implementation ready; requires schema confirmation
3. **Omics Type Property**: Property may not exist in all databases; handled gracefully with retry logic

## Next Steps
- Implement full Notion file attachment API
- Complete Program/Experiment relation updates when schema confirmed
- Add quantitative data support (intensity/abundance in future analysis)
- Add batch processing for multiple files
- Expand synonym mapping with external metabolite databases
- Add validation: file format and schema checking

