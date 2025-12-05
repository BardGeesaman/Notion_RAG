# Chemistry & HTS Integration Architecture

**Status**: 📋 Roadmap Item - Tier 2.6

**Priority**: 🔥🔥🔥 HIGH STRATEGIC VALUE

This document defines the architecture and implementation plan for integrating chemistry and High-Throughput Screening (HTS) data into the Amprenta multi-omics platform.

---

## 0. GENERAL PRINCIPLES

1. **SQLite is the system of record for chemistry and screening.**
   - All compound definitions
   - HTS-level data (up to 1M molecules)
   - Biochemical follow-up data

2. **Notion is the knowledge lens, not the raw chemistry store.**
   - Only "promoted" compounds (those reaching cell-based/animal assays) appear as **Compound Features** in Notion.
   - Notion stores Programs, Experiments, Datasets, Signatures, Features, and curated compound entities.

3. **Pinecone is a semantic index, not a DB.**
   - Only receive summarised text (campaign summaries, biological assay summaries, etc.).

4. **HTS & biochemical data should NOT be fully mirrored into Notion.**
   - Only summaries and selected hits.

---

## 1. SQLite CHEMISTRY + SCREENING LAYER

### 1.1 Module Structure

Create new modules:
- `amprenta_rag/chemistry/compound_registry.py` - Compound registration and normalization
- `amprenta_rag/chemistry/screening_db.py` - SQLite database operations
- `amprenta_rag/chemistry/normalization.py` - RDKit-based SMILES normalization (optional)

### 1.2 Database Schema

**Location**: `amprenta_chemistry.sqlite` (local file DB)

#### `compounds` Table
- `compound_id` (TEXT, primary key) - Internal AMPR ID, e.g., "AMPR-000001"
- `smiles` (TEXT)
- `inchi_key` (TEXT)
- `salt_stripped_smiles` (TEXT)
- `source` (TEXT) - vendor/library
- `created_at` (TIMESTAMP)
- `updated_at` (TIMESTAMP)

#### `libraries` Table
- `library_id` (TEXT, primary key)
- `name` (TEXT)
- `vendor` (TEXT)
- `description` (TEXT)

#### `hts_campaigns` Table
- `campaign_id` (TEXT, primary key)
- `program_id` (TEXT) - Maps to Notion Program page ID
- `assay_name` (TEXT)
- `target` (TEXT)
- `library_id` (TEXT)
- `screen_date` (TIMESTAMP)
- `num_molecules` (INTEGER)
- `hit_cutoff` (REAL or TEXT)
- `num_hits` (INTEGER)
- `qc_metrics_json` (TEXT) - JSON for Z', CV, etc.

#### `hts_results` Table
- `campaign_id` (TEXT)
- `compound_id` (TEXT)
- `plate_id` (TEXT)
- `well_id` (TEXT)
- `raw_signal` (REAL)
- `normalized_value` (REAL)
- `zscore` (REAL)
- `is_hit` (INTEGER/BOOLEAN)
- Composite primary key: `(campaign_id, compound_id, plate_id, well_id)`

#### `biochemical_results` Table
- `assay_id` (TEXT, primary key)
- `campaign_id` (TEXT, nullable)
- `compound_id` (TEXT)
- `ic50` (REAL, nullable)
- `ec50` (REAL, nullable)
- `curve_fit_quality` (REAL/TEXT)
- `replicate_count` (INTEGER)
- `qc_flags` (TEXT)
- `stage` (TEXT) - e.g., "primary confirmation", "orthogonal", "selectivity", etc.

#### `compound_program` Table
- `compound_id` (TEXT)
- `program_id` (TEXT)
- `role` (TEXT) - e.g., HIT, LEAD, TOOL, CANDIDATE
- Composite primary key: `(compound_id, program_id)`

**Indexes**: Create indexes on `compound_id`, `campaign_id`, and `program_id` for performance.

---

## 2. RDKit-BASED NORMALIZATION (IF AVAILABLE)

In `compound_registry.py`, add functions to:

1. **Standardize SMILES**:
   - Strip salts
   - Normalize tautomers where possible (optional)
   - Canonicalize

2. **Generate `inchi_key`**

3. **Assign or reuse an internal `compound_id`**:
   - Pattern: `AMPR-000001`, `AMPR-000002`, etc.
   - Maintain a counter in the DB.

**Fallback**: If RDKit is not available, document the assumption and implement a minimal version (SMILES + user-provided IDs).

---

## 3. NOTION SCHEMA INTEGRATION

### 3.1 Compound Features DB (Notion)

**Title**: `🧪 Compound Features`

**Properties**:
- `Name` (title) - Internal ID, e.g., AMPR-000123 or AMPR-001 (series)
- `SMILES` (text)
- `InChIKey` (text)
- `Series` (select)
- `Stage` (select: HTS Hit, Biochemical Hit, Cell-Based, Lead, Candidate, Tool)
- `Programs` (relation → Pipeline/Programs)
- `Experiments` (relation → Experiments)
- `Datasets` (relation → Experimental Data Assets)
- `Signatures` (relation → Signatures)
- `Notes` (text)

**Note**: Only compounds that pass biochemical triage and/or have in vitro/in vivo biology are "promoted" into this DB.

### 3.2 HTS Campaigns DB (Notion)

**Title**: `🧪 HTS Campaigns`

**Properties**:
- `Campaign Name` (title)
- `Program` (relation → Pipeline/Programs)
- `Assay Name` (text)
- `Target` (text)
- `Library` (text or relation)
- `Date` (date)
- `Molecules Screened` (number)
- `Hit Cutoff` (text/number)
- `Hits Found` (number)
- `QC Summary` (text)
- `Data Files` (files or links)
- `Notes` (text)

### 3.3 Biochemical Hits DB (Optional, Notion)

**Title**: `🧪 Biochemical Hits`

**Properties**:
- `Compound` (relation → Compound Features or text ID initially)
- `Program` (relation → Pipeline/Programs)
- `Assay Name` (text)
- `IC50` / `EC50` (number)
- `Curve Quality` (number/text)
- `Decision` (select: Promote, Drop, Investigate)
- `Notes` (text)

**Implementation Note**: Cursor should generate a Notion schema message (for the user to paste to the Notion agent) that creates these databases, similar to how other schemas have been handled.

---

## 4. SCREENING INGESTION MODULE

### 4.1 Module Structure

Create:
- `amprenta_rag/ingestion/screening_ingestion.py`
- `scripts/ingest_screening.py`

### 4.2 Responsibilities

- Ingest HTS campaign summary files (metadata, QC)
- Ingest hit lists and dose-response data into SQLite
- Optionally create/refresh HTS Campaigns DB entries in Notion with summary only (no full data matrix)
- Maintain mapping from vendor IDs → internal `compound_id`

### 4.3 CLI Use Cases

#### A. Register/Update Campaign
```bash
python scripts/ingest_screening.py \
  --campaign-metadata-file path/to/hts_campaign_meta.csv \
  --program-id <notion_program_id>
```

#### B. Ingest Hit List
```bash
python scripts/ingest_screening.py \
  --hit-list-file path/to/hits.csv \
  --campaign-id <campaign_id>
```

#### C. Ingest Dose-Response Data
```bash
python scripts/ingest_screening.py \
  --dose-response-file path/to/dr_data.csv \
  --campaign-id <campaign_id>
```

#### D. Promote Compounds to Notion
```bash
python scripts/ingest_screening.py \
  --promote-compounds \
  --program-id <notion_program_id> \
  --stage Biochemical Hit
```

---

## 5. COMPOUND-BIOLOGY LINKING

### 5.1 Integration Points

- Link compounds to Programs/Experiments/Datasets
- Link compounds to Signatures (if compound affects signature features)
- Cross-reference compound features with omics features
- RAG embedding of compound summaries (campaign summaries, biological assay summaries)

### 5.2 Module Structure

- `amprenta_rag/ingestion/compound_linking.py`
- Extend existing feature linking to support compounds
- Extend RAG embedding to include compound context

---

## 6. RAG INTEGRATION FOR CHEMISTRY

### 6.1 Features

- Embed HTS campaign summaries into Pinecone
- Embed biochemical hit summaries
- Query support: "Which compounds match this signature?", "What HTS campaigns tested this target?"
- Cross-omics + chemistry reasoning

### 6.2 Module Structure

- `amprenta_rag/query/chemistry_query.py`
- Extend `cross_omics_reasoning.py` to include compound context

---

## 7. IMPLEMENTATION PHASES

### Phase 1: SQLite Chemistry & Screening Layer (Weeks 13-14)
- Create database schema
- Implement compound registry
- RDKit normalization (if available)
- Database initialization and migration utilities

**Estimated Effort**: 8-10 days

### Phase 2: Screening Ingestion Pipeline (Week 14-15)
- Campaign metadata ingestion
- Hit list ingestion
- Dose-response data ingestion
- Notion integration (summary only)

**Estimated Effort**: 4-5 days

### Phase 3: Compound-Biology Linking (Week 15)
- Link compounds to Programs/Experiments/Datasets
- Link compounds to Signatures
- Cross-reference with omics features

**Estimated Effort**: 3-4 days

### Phase 4: RAG Integration for Chemistry (Week 16)
- Campaign summary embedding
- Biochemical hit summary embedding
- Chemistry query interface
- Cross-omics + chemistry reasoning

**Estimated Effort**: 3-4 days

---

## 8. DEPENDENCIES

- SQLite (Python `sqlite3` module)
- RDKit (optional, for SMILES normalization)
- Notion schema updates (Compound Features, HTS Campaigns, Biochemical Hits DBs)
- Existing multi-omics infrastructure (Programs, Experiments, Datasets, Signatures)

---

## 9. NOTION AGENT INSTRUCTIONS

When ready to implement, Cursor should generate a detailed Notion schema message describing:
1. Compound Features database structure
2. HTS Campaigns database structure
3. Biochemical Hits database structure (optional)
4. Relation properties to existing databases (Programs, Experiments, Datasets, Signatures)

This message should be formatted for the user to paste directly to the Notion agent.

---

**Last Updated**: 2025-12-04
**Status**: Ready for implementation

