# Internal Lipidomics Ingestion - Usage Guide

## Quick Start

### Basic Ingestion

**Create a new dataset page:**
```bash
python scripts/ingest_lipidomics.py \
  --file path/to/your/lipidomics.csv \
  --create-page
```

**Update an existing dataset page:**
```bash
python scripts/ingest_lipidomics.py \
  --file path/to/your/lipidomics.csv \
  --dataset-page-id <notion_page_id>
```

## File Format Requirements

### Supported Formats
- **CSV**: Comma-separated values (`.csv`)
- **TSV**: Tab-separated values (`.tsv` or `.txt` with tabs)

### Required Columns
- **Lipid Identity Column**: Must contain one of these names (case-insensitive):
  - `species`
  - `lipid`
  - `Lipid`
  - `Name`
  - `Molecule`
  - `compound`
  - `metabolite`

### Optional Columns
- `intensity`, `abundance`, `value` - Quantitative data (detected but not yet used in scoring)
- `condition`, `group`, `sample` - Sample metadata
- `fold_change`, `log2_fc`, `pvalue` - Statistical data

### Example Files

**Minimal format:**
```csv
species,intensity
Cer(d18:1/16:0),1234
SM(d18:1/24:1),5678
HexCer(d18:1/22:0),2345
```

**Vendor format:**
```csv
Lipid,Abundance
CER 16:0,1234
SM 24:1;O2,5678
hex_cer_24_0,6789
```

**With metadata:**
```csv
species,intensity,condition,fold_change
Cer(d18:1/16:0),1234,Control,1.0
Cer(d18:1/16:0),2345,Treatment,1.9
SM(d18:1/24:1),5678,Control,1.0
```

## Species Name Formats

### Supported Formats

The ingestion pipeline automatically normalizes various vendor formats to canonical species names:

| Input Format | Output (Canonical) | Notes |
|-------------|-------------------|-------|
| `CER 16:0` | `Cer(d18:1/16:0)` | Simple vendor format |
| `Cer 16:0` | `Cer(d18:1/16:0)` | Case variations |
| `SM 24:1;O2` | `SM(d18:1/24:1)` | Modifications removed |
| `hex_cer_24_0` | `HexCer(d18:1/24:0)` | Underscore format |
| `Cer(d18:1/16:0)+H` | `Cer(d18:1/16:0)` | Adducts removed |
| `SM(d18:1_24:1)` | `SM(d18:1/24:1)` | Separator normalized |
| `Cer d18:1/16:0` | `Cer(d18:1/16:0)` | With backbone |

### Already Canonical
If your file already uses canonical format (`Class(d18:1/chain)`), no normalization is needed:
- `Cer(d18:1/16:0)` ✅
- `SM(d18:1/24:1)` ✅
- `HexCer(d18:1/22:0)` ✅

## CLI Options

### Required Arguments
- `--file <path>`: Path to the lipidomics file (CSV/TSV)

### Page Management
- `--create-page`: Create a new Experimental Data Asset page
- `--dataset-page-id <id>`: Use existing page ID (alternative to `--create-page`)

**Note**: Either `--create-page` or `--dataset-page-id` must be provided.

### Linking (Optional)
- `--program-id <id>`: Link to Program page (can be specified multiple times)
- `--experiment-id <id>`: Link to Experiment page (can be specified multiple times)

### Example Commands

**Create new page:**
```bash
python scripts/ingest_lipidomics.py \
  --file data/my_lipidomics.csv \
  --create-page
```

**Update existing page:**
```bash
python scripts/ingest_lipidomics.py \
  --file data/updated_lipidomics.csv \
  --dataset-page-id 2beadf61-42ab-8167-bce5-d8cb15ced746
```

**With program/experiment links:**
```bash
python scripts/ingest_lipidomics.py \
  --file data/my_lipidomics.csv \
  --create-page \
  --program-id <program_page_id> \
  --experiment-id <experiment_page_id>
```

## Output

### Console Output
```
✅ Ingestion complete!
   Dataset page ID: 2beadf61-42ab-8167-bce5-d8cb15ced746
   File: test_data/test_canonical.csv
   Linked to 2 program(s)
   Linked to 1 experiment(s)
```

### Logs
The ingestion process logs detailed information:
- File parsing status
- Column detection
- Species normalization (success/warnings)
- Species extraction summary
- Page creation/update
- Signature matching results
- Embedding completion

### Notion Page
After ingestion, the Notion page will contain:
- **Experiment Name**: "Internal Lipidomics — <filename>"
- **Data Origin**: "Internal – Amprenta"
- **Dataset Source Type**: "Processed table"
- **Summary**: File path, row count, species count, signature matches
- **Signature Match Score**: Highest match score (if matches found)
- **Related Signature(s)**: Relations to matching signatures
- **Embedding IDs**: List of Pinecone vector IDs
- **Last Embedded**: Timestamp

## Troubleshooting

### Common Issues

**"Could not find lipid identity column"**
- Ensure your file has a column named: `species`, `lipid`, `Lipid`, `Name`, `Molecule`, `compound`, or `metabolite`
- Check for typos or extra spaces in column names

**"Could not normalize raw lipid"**
- Some vendor formats may not be supported yet
- Check logs for the specific lipid name
- Consider pre-processing your file to canonical format

**"No signature matches found"**
- This is normal if your dataset doesn't overlap with existing signatures
- The dataset will still be ingested and embedded
- Signature Match Score will remain empty

**"Error creating dataset page"**
- Check that `NOTION_EXP_DATA_DB_ID` is configured in `.env`
- Verify Notion API key is valid
- Check Notion API rate limits

### Error Messages

**File not found:**
```
FileNotFoundError: Lipidomics file not found: <path>
```
→ Check file path is correct

**No species extracted:**
```
ValueError: No species extracted from file: <path>
```
→ Check file format and column names

**Invalid arguments:**
```
error: Either --dataset-page-id must be provided or --create-page must be specified
```
→ Provide either `--create-page` or `--dataset-page-id`

## Best Practices

### File Preparation
1. **Use canonical format when possible**: Reduces normalization errors
2. **Include quantitative data**: Even if not used yet, it's preserved for future use
3. **Clean data**: Remove empty rows, handle missing values
4. **Consistent naming**: Use consistent lipid naming within a file

### Workflow
1. **Test with small file first**: Verify format and normalization
2. **Check logs**: Review normalization warnings
3. **Verify Notion page**: Confirm all data is correct
4. **Check signature matches**: Review match scores and components

### Data Quality
- **Normalize before ingestion**: Pre-process vendor formats if possible
- **Validate species names**: Ensure they match expected lipid classes
- **Check for duplicates**: Remove duplicate species entries if needed

## Integration with Existing Pipeline

### Consistency
The internal lipidomics pipeline uses the same:
- Signature scoring algorithm as MW datasets
- Notion database (Experimental Data Assets)
- RAG embedding pipeline
- Metadata structure

### Differences
- **Data Origin**: "Internal – Amprenta" vs "External – Published"
- **Source**: File path vs MW Study ID
- **Format**: CSV/TSV vs mwTab JSON

### Signature Matching
Both internal and external datasets:
- Use same overlap threshold (default: 0.3)
- Score against same signature database
- Write same Signature Match Score property
- Create same Related Signature(s) relations

## Future Enhancements

### Planned Features
- **Full file attachment**: Upload files to Notion (requires external storage)
- **Quantitative scoring**: Use intensity/abundance in signature scoring
- **Batch processing**: Process multiple files in one run
- **mzTab support**: Add mzTab format support
- **Validation**: File format and schema validation
- **Program/Experiment relations**: Complete implementation when schema confirmed

### Requesting Features
If you need additional features or encounter issues:
1. Check existing documentation
2. Review logs for error details
3. Test with minimal example file
4. Report specific use case and expected behavior

