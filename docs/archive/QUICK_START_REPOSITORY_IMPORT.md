# Quick Start: Import All Omics Data from Repositories

## 🚀 Fastest Way to Get Started

### Option 1: Simple Command (All Omics Types)

```bash
# Discover and import all omics types for a disease
python scripts/import_all_omics_repositories.py \
    --disease "ALS" \
    --import \
    --ingest
```

### Option 2: Specific Omics Type

```bash
# Import metabolomics studies
python scripts/import_all_omics_repositories.py \
    --omics-type metabolomics \
    --disease "Alzheimer" \
    --import \
    --ingest

# Import transcriptomics studies
python scripts/import_all_omics_repositories.py \
    --omics-type transcriptomics \
    --keywords "Parkinson" \
    --import \
    --ingest
```

### Option 3: Two-Step Process (Recommended for Large Imports)

```bash
# Step 1: Discover studies
python scripts/import_all_omics_repositories.py \
    --disease "ALS" \
    --discover-only \
    --output discovered_studies.json

# Step 2: Review and import
python scripts/import_all_omics_repositories.py \
    --input discovered_studies.json \
    --import \
    --ingest
```

---

## 📋 Available Repositories

- **GEO** - Transcriptomics (gene expression)
- **PRIDE** - Proteomics
- **MetaboLights** - Metabolomics
- **MW** - Metabolomics Workbench (all)
- **MW_LIPIDOMICS** - Lipidomics only
- **MW_METABOLOMICS** - Metabolomics only

List all repositories:
```bash
python scripts/import_all_omics_repositories.py --list-repositories
```

---

## 🎯 Common Use Cases

### Import All Omics for a Disease

```bash
python scripts/import_all_omics_repositories.py \
    --disease "ALS" \
    --max-results 30 \
    --import \
    --ingest
```

This will:
- Search all repositories for "ALS"
- Import up to 30 studies per repository
- Create datasets in Postgres
- Ingest to Pinecone for RAG queries

### Import Specific Omics Type

```bash
# Lipidomics only
python scripts/import_all_omics_repositories.py \
    --omics-type lipidomics \
    --disease "Alzheimer" \
    --import \
    --ingest

# Proteomics only
python scripts/import_all_omics_repositories.py \
    --omics-type proteomics \
    --keywords "Parkinson" \
    --import \
    --ingest
```

### Search by Keywords

```bash
python scripts/import_all_omics_repositories.py \
    --keywords "ceramide" "sphingolipid" \
    --repository MW_LIPIDOMICS \
    --import \
    --ingest
```

### Filter by Sample Type

```bash
python scripts/import_all_omics_repositories.py \
    --disease "ALS" \
    --sample-type "CSF" \
    --omics-type metabolomics \
    --import \
    --ingest
```

---

## 🔧 Command Options

### Discovery Options
- `--disease DISEASE` - Filter by disease
- `--omics-type TYPE` - Filter by omics type (transcriptomics, proteomics, metabolomics, lipidomics)
- `--repository REPO` - Search specific repository
- `--keywords KEYWORD ...` - Search keywords
- `--organism ORGANISM` - Filter by organism
- `--sample-type TYPE` - Filter by sample type
- `--max-results N` - Maximum results per repository (default: 50)

### Import Options
- `--import` - Import discovered studies
- `--ingest` - Ingest to Pinecone (enables RAG queries)
- `--create-notion` - Also create Notion pages (optional)
- `--discover-only` - Only discover, don't import
- `--input FILE` - Load studies from JSON file
- `--output FILE` - Save discovery results to JSON file
- `--stop-on-error` - Stop on first error

---

## 📊 What Happens During Import

1. **Repository Query** - Fetches study metadata from repository API
2. **Postgres Creation** - Creates Dataset record in Postgres
3. **Feature Extraction** - Extracts features (genes, proteins, metabolites, lipids)
4. **Feature Linking** - Links features to dataset in Postgres
5. **Pinecone Embedding** - Embeds study content for RAG queries

All data is stored in Postgres. Notion is optional.

---

## 💡 Tips

1. **Start Small** - Use `--max-results 5` to test first
2. **Review Before Import** - Use `--discover-only` to see what's available
3. **Filter Smart** - Use disease/organism/sample-type filters to narrow results
4. **Monitor Progress** - Watch logs for import status

---

## 📚 Full Documentation

See `docs/REPOSITORY_IMPORT_GUIDE.md` for complete documentation.

