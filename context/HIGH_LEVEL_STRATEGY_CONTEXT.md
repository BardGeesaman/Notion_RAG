# Amprenta Multi-Omics Platform — Full Technical Context for Cursor

**STATUS**: ✅ AUTHORITATIVE CONTEXT - USE AS REFERENCE FOR ALL FUTURE WORK

This file provides the complete architecture and system expectations for the Amprenta RAG-enabled multi-omics platform.  

Cursor should use this as the authoritative context for all future work.

---

# 1. PROJECT OVERVIEW

---

Amprenta has built a unified AI-native multi-omics knowledge system integrating:

- Notion (canonical knowledge + schema)
- Pinecone (vector index for RAG)
- Cursor (code execution + ingestion pipelines)
- OpenAI (embedding + reasoning)
- Multi-omics ingestion (lipidomics, metabolomics, proteomics, transcriptomics)
- Signature ingestion (now multi-omics)
- Feature-level knowledge graph (genes, proteins, metabolites, lipids)
- Program / Experiment / Dataset linking
- RAG Engine (chunk storage + retrieval)

Cursor is responsible for:

- Ingesting data from all omics pipelines
- Building normalized feature lists
- Creating/Updating Notion pages
- Linking features to datasets (and vice versa where schema supports)
- Extracting multi-omics signatures
- Scoring datasets against signatures
- Generating RAG chunks and embedding them
- Maintaining idempotency and non-blocking error handling

**Notion is the canonical source of truth for the knowledge graph.**

---

# 2. INGESTION PIPELINES (FULLY IMPLEMENTED)

---

Cursor manages ingestion for four omics types, each with its own ingestion module:

- `lipidomics_ingestion.py`
- `metabolomics_ingestion.py`
- `proteomics_ingestion.py`
- `transcriptomics_ingestion.py`

All ingestion pipelines follow the same pattern:

1. **Parse CSV/TSV**  
   Auto-detect delimiter. Flexible column detection.

2. **Normalize features**  
   - Lipids → canonical species (Cer(d18:1/16:0), SM(d18:1/24:1), etc.)
   - Metabolites → normalized metabolite names (Glutamate, Citrate)
   - Proteins → UniProt or canonical gene (P04637, TP53)
   - Genes → HGNC or Ensembl (TP53, ENSG00000141510)

3. **Create or update dataset page in Notion**  
   Database: **Experimental Data Assets**

4. **Link dataset → Experiments/Programs** (if provided)

5. **Feature linking (CRITICAL)**  
   Each normalized feature is linked to its feature DB:
   - Lipids → Lipid Species  
   - Metabolites → Metabolite Features  
   - Proteins → Protein Features  
   - Genes → Gene Features  

   A feature page is created if missing.  
   Relations are idempotent.

6. **Multi-omics signature scoring**  
   Every dataset is compared to all known signatures:
   - Score = weighted, direction-aware matching
   - Update dataset Notion fields:
     - `Related Signature(s)`
     - `Signature Match Score`
     - Signature overlap summary appended to Summary

7. **RAG embedding**  
   - Build text representation (features + summary + signature matches)
   - Chunk into ~2000 token blocks
   - Embed via OpenAI
   - Upsert vectors to Pinecone (batch upserts)
   - Write back:
     - `Embedding IDs`
     - `Last Embedded`

8. **Logging**  
   All ingestion uses consistent prefixes:
   - `[INGEST][LIPIDOMICS]`
   - `[INGEST][PROTEOMICS]`
   - `[INGEST][FEATURE]`
   - `[INGEST][SIGNATURE-MATCH]`
   - `[RAG]`

9. **Error handling (REQUIRED)**  
   - Never block ingestion on Notion failures for non-critical fields
   - Missing schema elements → warnings only
   - All ingestion runs must be idempotent

---

# 3. NOTION SCHEMA (REQUIRED STRUCTURE)

---

Cursor MUST expect these Notion databases with these key properties:

----------------------------------------------------

## 3.1 Experimental Data Assets (datasets)

Properties:

- Title: `Experiment Name`
- `Omics Type` (select: Lipidomics, Metabolomics, Proteomics, Transcriptomics, Multi-omics)
- `Data Origin` (select)
- `Dataset Source Type` (select)
- `Summary` (text)
- `Results` (text)
- `Conclusions` (text)
- `Embedding IDs` (text)
- `Last Embedded` (date)
- `Related Signature(s)` (relation → Signatures)
- `Signature Match Score` (number)
- `Lipid Species` (relation → Lipid Species)
- `Related Programs` (relation)
- `Experiments` (relation)
- Hidden RAG fields (`Chunks`, `Full Text`)

----------------------------------------------------

## 3.2 Feature Databases (CRITICAL)

Cursor must use these 4 DBs:

### Gene Features

- Name (title)
- External IDs (text)
- Synonyms (text)
- Species (multi-select)
- Transcriptomics Datasets (relation → Experimental Data Assets)
- Proteomics Datasets (optional)
- Experiments (relation)
- Signatures (relation)

### Protein Features

- Name (title)
- UniProt ID (text)
- Transcriptomics Datasets (optional)
- Proteomics Datasets (relation)
- Experiments (relation)
- Signatures (relation)

### Metabolite Features

- Name (title)
- Class (select)
- HMDB, KEGG IDs
- Metabolomics Datasets (relation)
- Experiments (relation)
- Signatures (relation)

### Lipid Species

- Name (title)
- Class (select)
- Synonyms (text)
- Experiments (relation)
- **Experimental Data Assets** (relation → datasets)
- Lipid Signature Components (relation)
- Notes (text)

----------------------------------------------------

## 3.3 Multi-Omics Signatures DB

Properties:

- `Title`
- `Short ID`
- `Signature Type`
- `Modalities` (multi-select: Gene, Protein, Metabolite, Lipid)
- `Components` (relation to signature components)
- `Embedding IDs`
- `Last Embedded`
- `Related Datasets`
- `Related Experiments`

----------------------------------------------------

## 3.4 Signature Components

Properties:

- `Component Name` (title)
- `Feature Type` (select: Gene, Protein, Metabolite, Lipid)
- `Raw Name`
- `Direction` (Up, Down, Neutral)
- `Weight` (number)
- `Feature` (relation → feature DB)
- `Signature` (relation → signature)

----------------------------------------------------

## 3.5 RAG Engine

Chunk DB with:

- `Chunk ID` (title)
- `Chunk Text`
- `Order`
- `Parent Type` (select: Literature, Dataset, Experiment, Note, Email)
- Relations:
  - Parent Item → Literature
  - Parent Email/Note Item → Email/Notes
  - Parent Experimental Data → Experimental Data Assets
  - Parent Research Note → Research Notes

---

# 4. MULTI-OMICS SIGNATURE SYSTEM (FULL BEHAVIOR)

---

Cursor must support ingestion of signature files (TSV/CSV) with:

Columns:

- feature_name (required)
- feature_type (optional; infer if missing)
- direction (optional)
- weight (optional)

Signature ingestion steps:

1. Load TSV/CSV  
2. Infer feature types (lipid/protein/metabolite/gene)  
3. Normalize each feature name  
4. Find/create feature pages  
5. Create signature component pages  
6. Link component → feature  
7. Link component → signature  
8. Update signature modalities  
9. Embed signature into Pinecone  
10. Write embedding metadata  
11. Score all datasets against signature  
12. Update dataset pages accordingly

All steps must be idempotent.

---

# 5. MULTI-OMICS SCORING ENGINE

---

Dataset scoring requires:

- Extract dataset feature sets by omics type
  - lipids → lipid species list
  - metabolites → metabolite list
  - proteins → protein list
  - genes → gene list

- For each signature:
  - Match features by type
  - Weighted + direction-aware scoring
  - Compute overlap fraction
  - Write back highest score to dataset

- Update:
  - `Related Signature(s)`
  - `Signature Match Score`
  - Overlap summary appended to `Summary`

---

# 6. RAG ENGINE BEHAVIOR

---

Cursor must:

- Build structured text representations for all entities
- Chunk text into ~2000 token segments
- Embed via OpenAI
- Upsert vectors into Pinecone using:
  - `source_type`
  - `dataset_page_id`
  - `experiment_page_id`
  - `signature_page_id`
  - `omics_type`

Cursor must not depend on specific RAG Engine schema fields.

---

# 7. CODING RULES

---

### Idempotency (MANDATORY)

- Never create duplicate features, signatures, components, datasets, or chunks
- Always check existing pages first

### Non-blocking (MANDATORY)

- Missing Notion properties → warn, skip
- Missing DB IDs → warn, skip
- API failures → warn and retry only where appropriate

### Clean logging (MANDATORY)

Use prefixes:

- `[INGEST][LIPIDOMICS]`
- `[INGEST][METABOLOMICS]`
- `[INGEST][PROTEOMICS]`
- `[INGEST][TRANSCRIPTOMICS]`
- `[INGEST][FEATURE]`
- `[INGEST][SIGNATURE-MATCH]`
- `[INGEST][SIGNATURE]`
- `[RAG]`

### File Structure Requirements

All ingestion modules should follow existing patterns in:

- `amprenta_rag/ingestion/*`
- `amprenta_rag/signatures/*`
- `amprenta_rag/query/*`

---

# 8. WHAT CURSOR MUST DO GOING FORWARD

---

Cursor should:

- Maintain and extend ingestion pipelines
- Extend multi-omics scoring
- Add cross-omics RAG reasoning functions
- Maintain schema compatibility with Notion
- Build new CLI tools following patterns
- Keep all code idempotent
- Keep ingestion resilient to schema drift
- Add new features modularly

Cursor should NOT:

- Break existing pipelines
- Change Notion schema without instruction
- Assume availability of any ingestion-controlled field
- Produce embeddings without metadata

---

# 9. CURRENT SYSTEM STATUS

---

## ✅ Fully Implemented & Operational

- **Multi-omics ingestion pipelines** (all 4 types)
- **Feature linking** (all 4 types)
- **Multi-omics signature ingestion**
- **Multi-omics signature scoring** (framework complete, dataset feature extraction integrated)
- **RAG embedding** (all sources)
- **Cross-omics RAG reasoning** (all 4 summary types)
- **Signature similarity queries**
- **Experiment ELN ingestion**

## 📊 Database Status

- ✅ All feature databases configured
- ✅ Signatures database working
- ✅ Experimental Data Assets working
- ✅ Programs database working (new single-source DB)
- ✅ Experiments database working

## 🎯 Key Principles

- **Idempotency**: All operations are safe to re-run
- **Non-blocking**: Errors don't stop ingestion
- **Schema resilience**: Graceful handling of missing properties
- **Comprehensive logging**: Clear prefixes for all operations

---

# End of Master Context Packet for Cursor

**Last Updated**: 2025-12-03  
**Status**: ✅ Authoritative - Use for all future work
