# Notion Database Setup Guide

This guide provides step-by-step instructions for creating the required Notion databases for the RAG system.

---

## 🧪 Chemistry & HTS Databases

### 1. Compound Features Database

**Purpose**: Store promoted compounds from HTS campaigns

**Steps**:
1. Create a new database in Notion
2. Name it: "🧪 Compound Features" (or "Compound Features")
3. Add the following properties:

| Property Name | Type | Required | Notes |
|--------------|------|----------|-------|
| Name | Title | ✅ Yes | Compound ID |
| SMILES | Rich Text | ✅ Yes | SMILES string |
| Canonical SMILES | Rich Text | ❌ No | Normalized SMILES |
| InChI Key | Rich Text | ❌ No | InChI key identifier |
| Molecular Formula | Rich Text | ❌ No | Chemical formula |
| Molecular Weight | Number | ❌ No | Molecular weight (Da) |
| LogP | Number | ❌ No | Lipophilicity |
| HBD Count | Number | ❌ No | Hydrogen bond donors |
| HBA Count | Number | ❌ No | Hydrogen bond acceptors |
| Rotatable Bonds | Number | ❌ No | Rotatable bond count |
| Related Programs | Relation | ❌ No | Link to Programs database |

4. Copy the database ID (from URL: `notion.so/YOUR_WORKSPACE/DATABASE_ID?v=...`)
5. Add to `.env`: `NOTION_COMPOUND_FEATURES_DB_ID=<database_id>`

---

### 2. HTS Campaigns Database

**Purpose**: Store HTS screening campaign summaries

**Steps**:
1. Create a new database in Notion
2. Name it: "🧪 HTS Campaigns" (or "HTS Campaigns")
3. Add the following properties:

| Property Name | Type | Required | Notes |
|--------------|------|----------|-------|
| Campaign Name | Title | ✅ Yes | Campaign name |
| Campaign ID | Rich Text | ✅ Yes | Unique campaign identifier |
| Description | Rich Text | ❌ No | Campaign description |
| Assay Type | Select | ❌ No | e.g., "Primary", "Secondary" |
| Target | Rich Text | ❌ No | Target protein/gene |
| Total Wells | Number | ❌ No | Total number of wells screened |
| Hit Count | Number | ❌ No | Number of hits identified |
| Run Date | Date | ❌ No | Date campaign was run |
| Related Programs | Relation | ❌ No | Link to Programs database |

4. Copy the database ID
5. Add to `.env`: `NOTION_HTS_CAMPAIGNS_DB_ID=<database_id>`

---

### 3. Biochemical Hits Database

**Purpose**: Store detailed biochemical assay results

**Steps**:
1. Create a new database in Notion
2. Name it: "🧪 Biochemical Hits" (or "Biochemical Hits")
3. Add the following properties:

| Property Name | Type | Required | Notes |
|--------------|------|----------|-------|
| Assay Name | Title | ✅ Yes | Name of the assay |
| Result ID | Rich Text | ✅ Yes | Unique result identifier |
| Compound | Relation | ✅ Yes | Link to Compound Features database |
| Target | Rich Text | ❌ No | Target protein/gene |
| IC50 | Number | ❌ No | IC50 value |
| EC50 | Number | ❌ No | EC50 value |
| Ki | Number | ❌ No | Ki value |
| Kd | Number | ❌ No | Kd value |
| Activity Type | Select | ❌ No | e.g., "IC50", "EC50", "Ki", "Kd" |
| Units | Rich Text | ❌ No | e.g., "nM", "μM" |
| Run Date | Date | ❌ No | Date assay was run |
| Related Programs | Relation | ❌ No | Link to Programs database |

4. Copy the database ID
5. Add to `.env`: `NOTION_BIOCHEMICAL_HITS_DB_ID=<database_id>`

---

## 🧬 Pathway Analysis Database (Optional)

### Pathways Database

**Purpose**: Store biological pathway information

**Steps**:
1. Create a new database in Notion
2. Name it: "🧬 Pathways" (or "Pathways")
3. Add the following properties:

| Property Name | Type | Required | Notes |
|--------------|------|----------|-------|
| Pathway Name | Title | ✅ Yes | Pathway name |
| Pathway ID | Rich Text | ✅ Yes | KEGG/Reactome ID |
| Source | Select | ✅ Yes | "KEGG" or "Reactome" |
| Description | Rich Text | ❌ No | Pathway description |
| Related Features | Relation | ❌ No | Link to feature databases |
| Related Datasets | Relation | ❌ No | Link to Experimental Data Assets |
| Related Signatures | Relation | ❌ No | Link to Signatures database |

4. Copy the database ID
5. Add to `.env`: `NOTION_PATHWAYS_DB_ID=<database_id>` (if implementing)

---

## ✅ Verification

After creating the databases, run:

```bash
python scripts/verify_notion_setup.py
```

This will verify that all required databases are configured and accessible.

---

## 📝 Notes

- Database IDs are the 32-character hex string in the Notion URL (without dashes)
- Relations can be set up after creating the databases
- Some properties are optional but recommended for full functionality
- The system will gracefully handle missing optional properties

