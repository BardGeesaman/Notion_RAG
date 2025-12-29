---
name: Publication Data Extraction
overview: Extract structured experiment data from scientific publication PDFs and supplementary files (Excel/CSV), linking extracted data to existing Literature and Dataset records.
todos:
  - id: pub-batch-1
    content: PDF Methods & Experiment Extraction (publication_extractor.py)
    status: in_progress
  - id: pub-batch-2
    content: Supplementary File Processing (supplementary_parser.py)
    status: pending
    dependencies:
      - pub-batch-1
  - id: pub-batch-3
    content: Database Schema & API endpoints
    status: pending
    dependencies:
      - pub-batch-2
  - id: pub-batch-4
    content: Dashboard Integration & Tests
    status: pending
    dependencies:
      - pub-batch-3
---

# Publication & Supplementary Data Extraction

## Overview

Build on existing extraction infrastructure to extract structured experiment data from scientific papers and their supplementary materials, then link this data to repository datasets.

## Current State

Existing infrastructure to leverage:
- [text_extraction.py](amprenta_rag/ingestion/text_extraction.py) - PDF text extraction (pypdf)
- [parsers.py](amprenta_rag/extraction/parsers.py) - DOCX/PPTX/Excel/CSV parsers
- [structured_extractor.py](amprenta_rag/extraction/structured_extractor.py) - LLM entity extraction
- [ai_extraction.py](scripts/dashboard/pages/ai_extraction.py) - AI Extraction dashboard
- Literature model with paper metadata

## Architecture

```
┌─────────────────────────────────────────────────────────────┐
│                    Input Sources                             │
├─────────────────────────────────────────────────────────────┤
│  Publication PDF  │  Supplementary Excel  │  Supplementary  │
│  (Methods, Figs)  │  (Data Tables)        │  CSV (Raw Data) │
└─────────────────────────────────────────────────────────────┘
                              │
                              ▼
┌─────────────────────────────────────────────────────────────┐
│              Extraction Pipeline (New)                       │
├─────────────────────────────────────────────────────────────┤
│  PublicationExtractor    │  SupplementaryParser             │
│  - Section detection     │  - Table detection               │
│  - Methods extraction    │  - Column type inference         │
│  - Figure/table parsing  │  - Data normalization            │
│  - LLM structured output │  - Schema matching               │
└─────────────────────────────────────────────────────────────┘
                              │
                              ▼
┌─────────────────────────────────────────────────────────────┐
│              Database (Extend)                               │
├─────────────────────────────────────────────────────────────┤
│  PublicationExtraction   │  SupplementaryFile               │
│  - literature_id (FK)    │  - literature_id (FK)            │
│  - experiment_type       │  - file_type                     │
│  - conditions            │  - extracted_tables              │
│  - measurements          │  - linked_dataset_id             │
└─────────────────────────────────────────────────────────────┘
```

## Batch Organization

### Batch 1: PDF Methods & Experiment Extraction

**Files to create:**
- `amprenta_rag/extraction/publication_extractor.py`

**Features:**
- Section detection (Abstract, Methods, Results, etc.)
- Experimental conditions extraction (cell lines, concentrations, timepoints)
- Assay type detection (qPCR, Western blot, RNA-seq, etc.)
- LLM-based structured extraction with Pydantic schemas
- Link to Literature record

**Pydantic schemas:**
```python
class PublicationExtraction(BaseModel):
    experiment_type: str  # RNA-seq, Western blot, ELISA, etc.
    cell_line: Optional[str]
    treatment: Optional[str]
    concentration: Optional[str]
    timepoint: Optional[str]
    replicate_count: Optional[int]
    measured_entities: List[str]  # genes, proteins, metabolites
```

**LLM Cost Optimization:**
- Extract Methods section only (not full PDF)
- Max 8000 tokens per extraction call
- Cache extraction results by literature_id + content_hash

### Batch 2: Supplementary File Processing

**Files to create:**
- `amprenta_rag/extraction/supplementary_parser.py`

**Features:**
- Auto-detect table schema (gene list, expression matrix, compound activity)
- Column type inference (gene symbol, p-value, fold change, IC50)
- Extract and normalize data tables
- Store as JSON with schema metadata

**Schema detection:**
- Gene expression: columns like "Gene", "log2FC", "padj"
- Compound activity: "Compound", "IC50", "Target"
- Proteomics: "Protein", "Abundance", "Sample"

### Batch 3: Database Schema & API

**Files to modify/create:**
- `amprenta_rag/models/content.py` - Add PublicationExtraction, SupplementaryFile models
- `alembic/versions/xxx_publication_extraction.py` - Migration
- `amprenta_rag/api/routers/papers.py` - Add extraction endpoints

**New endpoints:**
- POST /papers/{id}/extract - Extract experiments from PDF
- POST /papers/{id}/supplementary - Upload and parse supplementary file
- GET /papers/{id}/experiments - List extracted experiments
- POST /papers/{id}/link-dataset - Link supplementary to Dataset

### Batch 4: Dashboard Integration & Tests

**Files to modify/create:**
- `scripts/dashboard/pages/paper_search.py` - Add extraction UI
- `amprenta_rag/tests/extraction/test_publication_extractor.py`
- `amprenta_rag/tests/extraction/test_supplementary_parser.py`
- `amprenta_rag/tests/api/test_publication_extraction_api.py`

**UI features:**
- "Extract Experiments" button on paper detail
- Supplementary file upload
- View extracted data table
- Link to existing datasets

## Success Criteria

- Extract structured experiments from publication PDFs with >80% accuracy
- Parse supplementary Excel/CSV with auto-schema detection
- Link extracted data to Literature records
- 20+ tests covering extraction and API
- Zero @pytest.mark.skip (No Bandaids policy)
