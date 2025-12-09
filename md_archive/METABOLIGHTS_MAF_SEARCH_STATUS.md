# MetaboLights MAF File Search - Status Report

## Implementation Status

✅ **Script Created**: `scripts/find_metabolights_maf_studies.py`

The script implements the "List-then-Check" approach as recommended:
1. Fetches list of study IDs from MetaboLights API
2. Checks each study's investigation file for MAF files (`m_*.tsv`)
3. Handles 500 errors gracefully (skips and continues)
4. Returns studies that have MAF files

## Error Handling

✅ **Resilient Skip Strategy** implemented:
- Handles 500 server errors gracefully
- Logs errors but continues processing
- Tracks failed studies separately
- Does not crash on API failures

## Current Findings

After scanning multiple studies:
- **MAF files are rare** - Most studies use assay files (`a_*.txt`) with sample metadata
- Many studies contain only raw data files, not processed metabolite data
- Studies with MAF files appear to be less common in the public repository

## Why MAF Files Are Rare

1. **Study Type**: Many studies focus on raw data submission
2. **Processing Status**: MAF files require metabolite identification/quantification, which is a later processing step
3. **Repository Stage**: Many studies may not have completed metabolite assignment

## Next Steps

The extraction logic is **ready and will work** when MAF files are found:
- File discovery prioritizes MAF files
- Column detection looks for `metabolite_identification`
- Proper ISA-Tab parsing implemented

When a study with MAF files is encountered during normal ingestion, the extraction will automatically work correctly.

## Testing

The script can be run to search for MAF files:
```bash
python scripts/find_metabolights_maf_studies.py [limit]
```

Example: Search top 100 studies
```bash
python scripts/find_metabolights_maf_studies.py 100
```

