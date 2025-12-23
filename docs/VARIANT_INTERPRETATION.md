# Variant Interpretation

## Overview

This pipeline supports:

1. Ingesting a VEP `--tab` TSV as a `VariantSet`
2. Creating `Variant` rows (locus + gene + consequence/impact + gnomAD AF)
3. Downloading/reading ClinVar `variant_summary.txt.gz` and matching by (chr,pos,ref,alt)
4. Creating `VariantAnnotation` rows for ClinVar matches
5. Computing `GeneBurden` per gene for the set

## VEP TSV Format Requirements

The ingestion service expects a VEP `--tab`-style TSV with headers that include (at minimum):

- `Location` (e.g. `1:12345` or `1:12345-12345`)
- `Allele` (ALT allele)
- `Consequence`
- `IMPACT`

Optional columns that will be parsed when present:

- `SIFT` (e.g. `deleterious(0.02)` → stored as `deleterious`)
- `PolyPhen` (e.g. `benign(0.01)` → stored as `benign`)
- `gnomAD_AF`
- `CADD_PHRED`
- `SYMBOL` (gene symbol)
- `Uploaded_variation` / `Existing_variation` / `HGVSg` / `HGVSc` / `HGVSp`

Notes:

- The MVP ref/alt inference uses `Uploaded_variation` when possible. If not available, it falls back to `ref="N"` and `alt=Allele`.
- The Variant unique key is `(variant_set_id, chromosome, position, ref_allele, alt_allele)`.

## ClinVar Data Download

ClinVar is downloaded from:

- `variant_summary.txt.gz` (default):
  - `https://ftp.ncbi.nlm.nih.gov/pub/clinvar/tab_delimited/variant_summary.txt.gz`

The download URL can be overridden with:

- `CLINVAR_VARIANT_SUMMARY_URL`

The ingestion service caches the decompressed file at:

- `data/clinvar/variant_summary.txt`

## Burden Score Formula

Per gene, we compute:

- `n_variants`: count of variants in the set for that gene
- `n_pathogenic`: ClinVar `clinical_significance` contains `"pathogenic"`
- `n_vus`: contains `"uncertain"` or `"vus"`
- `n_benign`: contains `"benign"`

Burden score:

```text
burden_score = (n_pathogenic * 10) + n_vus
```

## API Usage Examples

### Ingest a VEP TSV into a VariantSet

```bash
curl -X POST "http://localhost:8000/api/variants/sets" \
  -H "Content-Type: application/json" \
  -d '{
    "name": "My Variant Set",
    "vep_tsv_path": "/absolute/path/to/vep_output.tsv"
  }'
```

### List variant sets

```bash
curl -X GET "http://localhost:8000/api/variants/sets"
```

### Query variants with filters

```bash
curl -X GET "http://localhost:8000/api/variants/sets/<set_id>/variants?gene=TP53&impact=HIGH&significance=pathogenic&limit=100"
```

### List gene burdens

```bash
curl -X GET "http://localhost:8000/api/variants/sets/<set_id>/burdens?min_score=5&limit=50"
```

### Run pathway enrichment on high-burden genes

```bash
curl -X POST "http://localhost:8000/api/variants/sets/<set_id>/pathway-enrichment" \
  -H "Content-Type: application/json" \
  -d '{"min_burden_score": 5}'
```


