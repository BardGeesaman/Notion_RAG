# Master Context Document for New Chat Sessions

**Purpose**: Upload this file to new Cursor chat sessions to preserve all context and continue work seamlessly.

**Usage**: 
1. Open a new chat in Cursor
2. Upload this file or paste its contents
3. Say: "I'm continuing work from a previous session. Please read this context and confirm you understand the system state."

**Last Updated**: 2025-12-03

**Update Command**: Run `python scripts/update_context.py` or say "update context" to Cursor

---

## 📋 Quick Reference: What This Document Contains

1. ✅ High-level strategic context
2. ✅ Current system architecture
3. ✅ Key configuration and database IDs
4. ✅ Recent work completed
5. ✅ Current system state
6. ✅ Important principles and patterns
7. ✅ Roadmap and next steps
8. ✅ File locations and structure

---

# 1. PROJECT OVERVIEW

## Project Name
**Amprenta Multi-Omics Platform** - RAG-enabled multi-omics knowledge system

## Core Components
- **Notion** - Canonical knowledge + schema (source of truth)
- **Pinecone** - Vector index for RAG
- **OpenAI** - Embedding + reasoning
- **Multi-omics ingestion** - Lipidomics, Metabolomics, Proteomics, Transcriptomics
- **Multi-omics signatures** - Gene, Protein, Metabolite, Lipid features
- **Feature-level knowledge graph** - Individual feature pages (genes, proteins, metabolites, lipids)
- **Program / Experiment / Dataset linking**
- **RAG Engine** - Chunk storage + retrieval + cross-omics reasoning

## Key Principle
**Notion is the canonical source of truth for the knowledge graph.**

---

# 2. HIGH-LEVEL STRATEGIC CONTEXT

## Authoritative Reference
**File**: `context/HIGH_LEVEL_STRATEGY_CONTEXT.md`

### Core Principles
1. ✅ **Idempotency** - All operations safe to re-run
2. ✅ **Non-blocking** - Errors don't stop ingestion
3. ✅ **Schema resilience** - Graceful property handling
4. ✅ **Comprehensive logging** - Clear prefixes (`[INGEST]`, `[RAG]`, etc.)
5. ✅ **Notion as canonical source** - All truth comes from Notion

### System Expectations
- All ingestion pipelines are idempotent
- Error handling is graceful and non-fatal
- Logging is consistent and informative
- Notion schema changes are coordinated (use agent instructions)

---

# 3. CURRENT SYSTEM ARCHITECTURE

## Ingestion Pipelines (All Implemented)

### Multi-Omics Ingestion
1. **Lipidomics**: `amprenta_rag/ingestion/lipidomics_ingestion.py`
   - Supports CSV/TSV
   - Species normalization (vendor formats → canonical)
   - Feature linking to Lipid Species DB
   - Signature scoring integration

2. **Metabolomics**: `amprenta_rag/ingestion/metabolomics_ingestion.py`
   - Supports CSV/TSV
   - Metabolite name normalization
   - Feature linking to Metabolite Features DB
   - Signature scoring integration

3. **Proteomics**: `amprenta_rag/ingestion/proteomics_ingestion.py`
   - Supports CSV/TSV
   - Protein identifier normalization
   - Feature linking to Protein Features DB
   - Signature scoring integration

4. **Transcriptomics**: `amprenta_rag/ingestion/transcriptomics_ingestion.py`
   - Supports CSV/TSV DGE tables
   - Gene identifier normalization
   - Feature linking to Gene Features DB
   - DGE summary in RAG chunks
   - Signature scoring integration

### Shared Utilities
- `amprenta_rag/ingestion/omics_ingestion_utils.py` - Common functions
- `amprenta_rag/ingestion/feature_extraction.py` - Feature linking

### Signature System (Multi-Omics)
- **Loader**: `amprenta_rag/signatures/signature_loader.py`
- **Feature Type Inference**: `amprenta_rag/signatures/feature_type_inference.py`
- **Notion CRUD**: `amprenta_rag/ingestion/signature_notion_crud.py`
- **Linking**: `amprenta_rag/ingestion/signature_linking.py`
- **Embedding**: `amprenta_rag/ingestion/signature_embedding.py`
- **Matching/Scoring**: `amprenta_rag/ingestion/signature_matching.py`
- **Multi-Omics Scoring**: `amprenta_rag/ingestion/multi_omics_scoring.py`

### RAG Engine
- **Query Engine**: `amprenta_rag/query/rag_engine.py`
- **Cross-Omics Reasoning**: `amprenta_rag/query/cross_omics_reasoning.py`
- **Query CLI**: `scripts/rag_query.py`

### Dataset Ingestion
- **MW Datasets**: `amprenta_rag/ingestion/dataset_ingestion.py`
- **Experiments**: `amprenta_rag/ingestion/experiments_ingestion.py`

---

# 4. KEY CONFIGURATION

## Configuration File
**Location**: `amprenta_rag/config.py`

## Essential Environment Variables (.env file)

### Notion Database IDs
```
NOTION_LITERATURE_DB_ID=
NOTION_EXPERIMENTAL_DATA_DB_ID=
NOTION_SIGNATURE_DB_ID=
NOTION_SIGNATURE_COMPONENT_DB_ID=
NOTION_LIPID_SPECIES_DB_ID=
NOTION_METABOLITE_FEATURES_DB_ID=
NOTION_PROTEIN_FEATURES_DB_ID=
NOTION_GENE_FEATURES_DB_ID=
NOTION_PROGRAMS_DB_ID=bde04fc4ed6640a1b392445d7e1a08ae
NOTION_EXPERIMENTS_DB_ID=
NOTION_EMAIL_DB_ID=
```

### API Keys
```
NOTION_API_KEY=
PINECONE_API_KEY=
OPENAI_API_KEY=
```

### Directories
```
SIGNATURES_DIR=
```

### Pipeline Configuration
```
SIGNATURE_OVERLAP_THRESHOLD=0.3
ENABLE_SIGNATURE_SCORING=True
ENABLE_LIPID_MAPPING=True
```

---

# 5. NOTION DATABASE SCHEMAS

## Experimental Data Assets
- **Key Properties**: Name, Data Origin, Dataset Source Type, Summary
- **Relations**: Related Programs, Related Experiments, Related Signature(s)
- **Metrics**: Signature Match Score (number), Embedding IDs, Last Embedded
- **File Attachment**: Data File property

## Lipid Signatures
- **Key Properties**: Name, Short ID, Status, Signature Type, Data Ownership
- **Multi-Omics**: Modalities (multi-select: Gene, Protein, Metabolite, Lipid)
- **Relations**: Components, Related Programs

## Signature Components
- **Key Properties**: Component Name, Feature Type (select), Direction, Weight
- **Relations**: Signature, Feature (gene/protein/metabolite/lipid)

## Feature Databases
- **Lipid Species**: Name, Canonical Name, Lipid Class
- **Metabolite Features**: Name, Normalized Name
- **Protein Features**: Name, UniProt ID (optional)
- **Gene Features**: Name, Gene Symbol

All feature databases have relation properties linking back to datasets.

---

# 6. RECENT WORK COMPLETED

## ✅ Fully Implemented and Working

1. **Multi-Omics Signature System**
   - Feature type inference (gene, protein, metabolite, lipid)
   - Multi-omics signature ingestion
   - Component linking to feature databases
   - Signature embedding with modalities

2. **Multi-Omics Scoring Engine**
   - Dataset feature extraction by omics type
   - Multi-omics signature scoring
   - Signature match writeback to datasets

3. **Cross-Omics RAG Reasoning**
   - `cross_omics_program_summary()`
   - `cross_omics_signature_summary()`
   - `cross_omics_feature_summary()`
   - `cross_omics_dataset_summary()`

4. **Feature Linking (All Omics Types)**
   - Automatic feature page creation/linking
   - Dynamic relation property handling
   - Idempotent operations

5. **Programs and Experiments DB Support**
   - Programs database configured and working
   - Experiments database support
   - Cross-omics reasoning integrated

6. **Unified Strategic Roadmap**
   - Integrated ChatGPT's roadmap
   - 30+ enhancements organized into 5 tiers
   - Implementation timelines

## 🔧 Recently Fixed Issues

1. **Programs DB Multiple Data Sources** - Notion API limitation identified, workaround implemented
2. **Metabolite Features DB** - Configuration updated and verified
3. **Feature Relation Properties** - Dynamic property detection implemented
4. **Circular Imports** - Refactored shared utilities

---

# 7. CURRENT SYSTEM STATE

## Implementation Status

### ✅ Complete (Production Ready)
- Multi-omics ingestion pipelines (all 4 types)
- Feature linking (all 4 types)
- Multi-omics signature ingestion
- Multi-omics signature scoring
- Cross-omics RAG reasoning
- All database configurations

### 🚧 In Progress
- None currently

### 📋 Planned (See Roadmap)
- Feature caching (Tier 1)
- Batch ingestion (Tier 1)
- Signature discovery (Tier 2)
- Pathway analysis (Tier 2)

## Known Limitations

1. **Dataset Feature Extraction** - Currently queries Notion on-demand (performance bottleneck)
2. **Pathway Analysis** - Not yet implemented (requires Notion schema changes)
3. **Public Repository Ingestion** - Stretch goal

---

# 8. KEY FILE LOCATIONS

## Core Modules
```
amprenta_rag/
├── config.py                          # Configuration (DB IDs, API keys)
├── ingestion/
│   ├── lipidomics_ingestion.py       # Lipidomics pipeline
│   ├── metabolomics_ingestion.py     # Metabolomics pipeline
│   ├── proteomics_ingestion.py       # Proteomics pipeline
│   ├── transcriptomics_ingestion.py  # Transcriptomics pipeline
│   ├── omics_ingestion_utils.py      # Shared utilities
│   ├── feature_extraction.py         # Feature linking
│   ├── signature_matching.py         # Signature scoring
│   ├── multi_omics_scoring.py        # Multi-omics scoring
│   └── ...
├── signatures/
│   ├── signature_loader.py           # Signature loading
│   └── feature_type_inference.py     # Feature type detection
└── query/
    ├── rag_engine.py                 # RAG orchestration
    └── cross_omics_reasoning.py      # Cross-omics summaries
```

## Scripts
```
scripts/
├── ingest_lipidomics.py
├── ingest_metabolomics.py
├── ingest_proteomics.py
├── ingest_transcriptomics.py
├── ingest_signature.py
├── ingest_dataset.py
├── ingest_experiment.py
├── rag_query.py
└── list_cross_omics_test_ids.py
```

## Documentation
```
context/
├── HIGH_LEVEL_STRATEGY_CONTEXT.md        # Authoritative strategy
├── UNIFIED_STRATEGIC_ROADMAP.md          # Complete roadmap
├── STRATEGIC_FEATURE_ENHANCEMENTS.md     # Strategic analysis
├── FEATURE_ENHANCEMENT_TOP_3_DETAILED.md # Top 3 detailed specs
├── ROADMAP_INTEGRATION_SUMMARY.md        # Integration reference
├── MASTER_CONTEXT_FOR_NEW_CHAT.md        # This file
├── NEW_CHAT_QUICK_START.md               # Quick start guide
└── UPDATE_CONTEXT_GUIDE.md               # How to update context files
```

---

# 9. IMPORTANT PATTERNS & CONVENTIONS

## Logging Prefixes
- `[INGEST][LIPIDOMICS]` - Lipidomics ingestion
- `[INGEST][METABOLOMICS]` - Metabolomics ingestion
- `[INGEST][PROTEOMICS]` - Proteomics ingestion
- `[INGEST][TRANSCRIPTOMICS]` - Transcriptomics ingestion
- `[INGEST][SIGNATURES]` - Signature ingestion
- `[INGEST][SIGNATURE-MATCH]` - Signature matching/scoring
- `[INGEST][MULTI-OMICS-SCORE]` - Multi-omics scoring
- `[RAG][CROSS-OMICS]` - Cross-omics reasoning
- `[SIGNATURE][FEATURE-TYPE]` - Feature type inference

## Error Handling
- All errors are non-blocking
- Graceful degradation (log warnings, continue)
- Idempotent operations (safe to re-run)

## Notion API Patterns
- Always use `timeout=30` for requests
- Handle missing properties gracefully
- Use dynamic property detection when possible

## Feature Normalization
- Lipidomics: `normalize_lipid_species()` in `signature_matching.py`
- Metabolomics: `normalize_metabolite_name()` in `metabolomics_ingestion.py`
- Proteomics: `normalize_protein_identifier()` in `proteomics_ingestion.py`
- Transcriptomics: `normalize_gene_identifier()` in `transcriptomics_ingestion.py`

---

# 10. ROADMAP STATUS

## Current Roadmap Document
**File**: `context/UNIFIED_STRATEGIC_ROADMAP.md`

## Priority Tiers

### Tier 1: Immediate High Value (Next 2-4 Weeks)
1. Enhanced Dataset Feature Extraction & Caching (Performance)
2. Batch Ingestion Framework (Operations)
3. Enhanced Cross-Omics Reasoning (Better summaries)

### Tier 2: Strategic Capabilities (Next 4-8 Weeks)
1. Automated Signature Discovery
2. Evidence Report Engine
3. Program-Level Signature Maps
4. Cross-Omics Pathway Analysis
5. Dataset Comparison & Clustering

### Tier 3: Quality & Operations (Ongoing)
1. Quality Control Extraction
2. Signature Validation & Quality Metrics
3. Cross-Feature Mapping
4. Performance Optimizations

See `UNIFIED_STRATEGIC_ROADMAP.md` for complete details.

---

# 11. TESTING & VERIFICATION

## Test Scripts
- `scripts/list_cross_omics_test_ids.py` - List available test IDs
- `scripts/rag_query.py` - RAG query interface

## Verification Commands
```bash
# List available test IDs
python scripts/list_cross_omics_test_ids.py --type all

# Test cross-omics summaries
python scripts/rag_query.py --cross-omics-program <program_id>
python scripts/rag_query.py --cross-omics-signature <signature_id>
python scripts/rag_query.py --cross-omics-feature "gene:TP53"

# Test signature scoring
python scripts/rag_query.py --signature-score <dataset_page_id>
```

---

# 12. COMMON WORKFLOWS

## Ingest a New Omics Dataset
```bash
# Lipidomics
python scripts/ingest_lipidomics.py --file <path> --create-page

# Metabolomics
python scripts/ingest_metabolomics.py --file <path> --create-page

# Proteomics
python scripts/ingest_proteomics.py --file <path> --create-page

# Transcriptomics
python scripts/ingest_transcriptomics.py --file <path> --create-page
```

## Ingest a Multi-Omics Signature
```bash
python scripts/ingest_signature.py --signature-file <path>
```

## Query RAG
```bash
python scripts/rag_query.py --query "your question"
```

---

# 13. TROUBLESHOOTING GUIDE

## Common Issues

### Missing Database IDs
- Check `.env` file
- Verify database IDs are correct
- Use `scripts/list_cross_omics_test_ids.py` to verify access

### Notion API Errors
- Check API key is valid
- Verify database permissions
- Check for rate limiting (429 errors)
- Use dynamic property detection for schema changes

### Feature Linking Failures
- Check feature database IDs are configured
- Verify relation properties exist on feature pages
- Check logs for specific error messages

### Signature Scoring Not Working
- Verify `ENABLE_SIGNATURE_SCORING=True` in config
- Check signature components are linked to features
- Verify dataset features are extracted correctly

---

# 14. NOTION AGENT INSTRUCTIONS

## When Schema Changes Are Needed

Document the requirements clearly:
1. Database name
2. Property names and types
3. Relation targets
4. Property options (for selects)

## Current Pending Instructions
- **Pathways Database** (for Pathway Analysis) - See `FEATURE_ENHANCEMENT_TOP_3_DETAILED.md`
- **QC Properties** (for Quality Control) - See `UNIFIED_STRATEGIC_ROADMAP.md`

---

# 15. KEY PRINCIPLES FOR CONTINUING WORK

## When Starting New Work

1. ✅ **Check existing implementations first** - Avoid duplication
2. ✅ **Maintain idempotency** - All operations safe to re-run
3. ✅ **Use consistent logging** - Follow existing prefixes
4. ✅ **Handle errors gracefully** - Non-blocking, log warnings
5. ✅ **Preserve Notion schema** - Only change with agent instructions
6. ✅ **Test with real data** - Verify with actual datasets
7. ✅ **Document changes** - Update relevant .md files

## When Adding New Features

1. Check the unified roadmap first
2. Follow existing patterns
3. Integrate with Notion as source of truth
4. Add appropriate logging
5. Include error handling
6. Test idempotency

## When Fixing Issues

1. Check logs for error patterns
2. Verify configuration (DB IDs, API keys)
3. Test with minimal examples
4. Preserve existing functionality
5. Document the fix

---

# 16. QUICK START FOR NEW CHAT

## What to Say in New Chat

```
I'm continuing work on the Amprenta Multi-Omics Platform. I've uploaded 
the master context document. Please:

1. Read and confirm you understand the system architecture
2. Confirm you understand the current implementation status
3. Tell me what you see as the next logical steps based on the roadmap

The system uses Notion as canonical source of truth, all ingestion is 
idempotent, and errors should be non-blocking. All work should align 
with HIGH_LEVEL_STRATEGY_CONTEXT.md principles.
```

## What to Expect

The AI should:
- ✅ Acknowledge understanding of the system
- ✅ Reference the strategic context
- ✅ Understand current state
- ✅ Propose next steps aligned with roadmap
- ✅ Maintain architectural patterns

---

# 17. UPDATING THIS DOCUMENT

## Automatic Updates

Run the update script:
```bash
python scripts/update_context.py
```

Or say to Cursor: **"update context"**

This automatically:
- Archives extraneous .md files
- Updates timestamp
- Verifies context files

## Manual Updates Needed

Update this master context document when:
- ✅ Major features are completed
- ✅ Configuration changes are made
- ✅ Architecture decisions are finalized
- ✅ New databases are added
- ✅ Significant issues are resolved
- ✅ Roadmap priorities change

## How to Update

1. Run `python scripts/update_context.py` first
2. Then manually edit this file:
   - Add new sections as needed
   - Update "Recent Work Completed"
   - Update "Current System State"
   - Update roadmap status
   - Add new file locations
   - Document new patterns

See `context/UPDATE_CONTEXT_GUIDE.md` for detailed instructions.

---

# 18. CONTACT & CONTINUITY

## If Context Is Lost

1. Upload this document
2. Reference `HIGH_LEVEL_STRATEGY_CONTEXT.md`
3. Check `UNIFIED_STRATEGIC_ROADMAP.md` for priorities
4. Review recent `.md` files for latest work

## Key Reference Files (Read These First)

1. `context/HIGH_LEVEL_STRATEGY_CONTEXT.md` - Strategic principles
2. `context/UNIFIED_STRATEGIC_ROADMAP.md` - Complete roadmap
3. `context/MASTER_CONTEXT_FOR_NEW_CHAT.md` - This file

---

# ✅ END OF MASTER CONTEXT

**Remember**: This system is production-ready for multi-omics ingestion, 
signature scoring, and cross-omics reasoning. All enhancements should 
build on this solid foundation while maintaining idempotency, graceful 
error handling, and Notion as the canonical source of truth.

**Next Steps**: See `UNIFIED_STRATEGIC_ROADMAP.md` Tier 1 items for 
immediate high-value enhancements.

