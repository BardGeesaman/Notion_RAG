# MetaboLights MAF (Metabolite Assignment File) Extraction

## Overview

MetaboLights uses the **ISA-Tab (Investigation-Study-Assay)** standard. The metabolite identification and quantification data is stored in **Metabolite Assignment File (MAF)** files with naming pattern `m_MTBLS[ID]_...v2_maf.tsv`.

## MAF File Structure

### Fixed Metadata Columns (Left Side)

The first set of columns describes the chemical identity of the metabolite:

- `database_identifier`: External ID (e.g., `CHEBI:17234`, `HMDB0000123`)
- `chemical_formula`: Chemical formula (e.g., `C6H12O6`)
- `smiles`: SMILES string structure representation
- `inchi`: InChI string
- `metabolite_identification`: Common name (e.g., "Glucose") ⭐ **Primary column for extraction**
- `mass_to_charge`: Observed m/z value
- `retention_time`: Time the metabolite eluted
- `reliability`: Confidence score (1-4, 1 is highest/confirmed standard)

### Sample/Abundance Columns (Right Side)

Following the metadata columns, the file lists every **Sample Name** used in the study:
- Headers must match the `Sample Name` entries in the `s_*.txt` file
- Cells contain quantification values (peak area, intensity, concentration)

## Implementation

### File Discovery Priority

1. **MAF files** (`m_*_maf.tsv`) - Highest priority, contains metabolite data
2. **Other m_ files** (`m_*.tsv`) - Fallback for metabolite data files
3. **Assay files** (`a_*.txt`) - Last resort (contains sample metadata, not metabolite data)

### Column Detection Logic

The extraction automatically identifies the metabolite identification column with this priority:

1. **`metabolite_identification`** - Primary column for metabolite names
2. **`database_identifier`** - Fallback for database IDs (CHEBI, HMDB, etc.)
3. **`chemical_formula`** or other compound identifiers - Last resort

The extraction stops at sample columns (when "Sample Name" or "sample_" columns are detected) to avoid extracting sample IDs as metabolites.

### Extraction Process

1. **Find MAF file**: Checks investigation file for `m_*_maf.tsv` files
2. **Stream download**: Uses efficient streaming approach (same as GEO/PRIDE)
3. **Parse header**: Identifies metabolite identification column
4. **Extract metabolites**: Reads metadata columns only (stops at sample columns)
5. **Normalize**: Uses `normalize_metabolite_name()` to canonicalize names
6. **Link to dataset**: Links extracted metabolites to Postgres dataset

## Usage

```python
from amprenta_rag.ingestion.repository_feature_extraction import extract_metabolights_metabolites_from_isa_tab

# Extract metabolites from a MetaboLights study
metabolite_set = extract_metabolights_metabolites_from_isa_tab(study_id="MTBLS1")

# Results are normalized metabolite names ready for linking
print(f"Found {len(metabolite_set)} unique metabolites")
```

## Testing

```bash
# Test with default study (MTBLS1)
python scripts/test_metabolights_feature_extraction.py

# Test with specific study
python scripts/test_metabolights_feature_extraction.py MTBLS200
```

## Notes

- MAF files are TSV (tab-separated) format
- Sample columns are automatically skipped
- Empty values, "NA", and sample names are filtered out
- Metabolite names are normalized using the standard normalization function

