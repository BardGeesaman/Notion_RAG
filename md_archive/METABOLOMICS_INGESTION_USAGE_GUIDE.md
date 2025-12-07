# Internal Metabolomics Ingestion - Usage Guide

## Quick Start

### Basic Ingestion

**Create a new dataset page:**
```bash
python scripts/ingest_metabolomics.py \
  --file path/to/your/metabolomics.csv \
  --create-page
```

**Update an existing dataset page:**
```bash
python scripts/ingest_metabolomics.py \
  --file path/to/your/metabolomics.csv \
  --dataset-page-id <notion_page_id>
```

## File Format Requirements

### Supported Formats
- **CSV**: Comma-separated values (`.csv`)
- **TSV**: Tab-separated values (`.tsv` or `.txt` with tabs)

### Required Columns
- **Metabolite Identity Column**: Must contain one of these names (case-insensitive):
  - `metabolite`
  - `Metabolite`
  - `compound`
  - `Compound`
  - `name`
  - `Name`
  - `molecule`
  - `Molecule`

### Optional Columns
- `intensity`, `abundance`, `area` - Quantitative data (detected but not yet used)
- `sample`, `group`, `condition` - Sample metadata
- `fold_change`, `log2_fc`, `pvalue` - Statistical data

### Example Files

**Minimal format:**
```csv
metabolite,intensity
L-Glutamine,1234
L-Glutamic acid,5678
L-Serine,2345
```

**Mixed formatting:**
```csv
Name,Abundance
glutamine,1234
GLUTAMIC ACID,5678
Serine [M+H]+,2345
```

**With metadata:**
```csv
metabolite,intensity,condition,fold_change
Glutamine,1234,Control,1.0
Glutamine,2345,Treatment,1.9
Serine,5678,Control,1.0
```

## Metabolite Name Formats

### Supported Formats

The ingestion pipeline automatically normalizes various formats to canonical metabolite names:

| Input Format | Output (Canonical) | Notes |
|-------------|-------------------|-------|
| `L-Glutamine` | `Glutamine` | L-prefix removed via synonym |
| `L-Glutamic acid` | `Glutamate` | Synonym mapping |
| `glutamine` | `Glutamine` | Case normalization |
| `GLUTAMIC ACID` | `Glutamate` | Case + synonym |
| `Serine [M+H]+` | `Serine` | Adduct removed |
| `Glutamine (pos)` | `Glutamine` | Annotation removed |
| `L_Serine` | `Serine` | Underscore → space |

### Normalization Features
- **Adduct Removal**: `[M+H]`, `[M-H]`, `[M+Na]`, etc.
- **Annotation Removal**: `(pos)`, `(neg)`, `(+)`, `(-)`
- **Underscore Handling**: `_` → space
- **Synonym Mapping**: 20+ common metabolites (L-Glutamic acid → Glutamate, etc.)
- **Case Normalization**: Title Case for common metabolites

## CLI Options

### Required Arguments
- `--file <path>`: Path to the metabolomics file (CSV/TSV)

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
python scripts/ingest_metabolomics.py \
  --file data/my_metabolomics.csv \
  --create-page
```

**Update existing page:**
```bash
python scripts/ingest_metabolomics.py \
  --file data/updated_metabolomics.csv \
  --dataset-page-id 2beadf61-42ab-81bf-8544-fd2c918be46f
```

**With program/experiment links:**
```bash
python scripts/ingest_metabolomics.py \
  --file data/my_metabolomics.csv \
  --create-page \
  --program-id <program_page_id> \
  --experiment-id <experiment_page_id>
```

## Output

### Console Output
```
✅ Ingestion complete!
   Dataset page ID: 2beadf61-42ab-81bf-8544-fd2c918be46f
   File: test_data/test_metabolomics_canonical.csv
   Linked to 2 program(s)
   Linked to 1 experiment(s)
```

### Logs
The ingestion process logs detailed information:
- File parsing status
- Column detection
- Metabolite normalization (success)
- Metabolite extraction summary
- Page creation/update
- Embedding completion

### Notion Page
After ingestion, the Notion page will contain:
- **Experiment Name**: "Internal Metabolomics — <filename>"
- **Data Origin**: "Internal – Amprenta"
- **Dataset Source Type**: "Processed table"
- **Omics Type**: "Metabolomics" (if property exists)
- **Summary**: File path, row count, metabolite count
- **Embedding IDs**: List of Pinecone vector IDs
- **Last Embedded**: Timestamp

## Troubleshooting

### Common Issues

**"Could not find metabolite identity column"**
- Ensure your file has a column named: `metabolite`, `Metabolite`, `compound`, `Compound`, `name`, `Name`, `molecule`, or `Molecule`
- Check for typos or extra spaces in column names

**"Omics Type is not a property that exists"**
- This is normal - the property may not exist in your database
- The system automatically retries without this property
- Ingestion continues successfully

**"No metabolites extracted"**
- Check file format and column names
- Verify file is not empty
- Check for encoding issues (should be UTF-8)

**"Error creating dataset page"**
- Check that `NOTION_EXP_DATA_DB_ID` is configured in `.env`
- Verify Notion API key is valid
- Check Notion API rate limits

### Error Messages

**File not found:**
```
FileNotFoundError: Metabolomics file not found: <path>
```
→ Check file path is correct

**No metabolites extracted:**
```
ValueError: No metabolites extracted from file: <path>
```
→ Check file format and column names

**Invalid arguments:**
```
error: Either --dataset-page-id must be provided or --create-page must be specified
```
→ Provide either `--create-page` or `--dataset-page-id`

## Best Practices

### File Preparation
1. **Use consistent naming**: Use consistent metabolite naming within a file
2. **Include quantitative data**: Even if not used yet, it's preserved for future use
3. **Clean data**: Remove empty rows, handle missing values
4. **Standard formats**: Prefer standard metabolite names when possible

### Workflow
1. **Test with small file first**: Verify format and normalization
2. **Check logs**: Review normalization results
3. **Verify Notion page**: Confirm all data is correct
4. **Check embeddings**: Verify RAG can retrieve content

### Data Quality
- **Normalize before ingestion**: Pre-process metabolite names if possible
- **Validate metabolite names**: Ensure they match expected formats
- **Check for duplicates**: Remove duplicate metabolite entries if needed

## Integration with Existing Pipeline

### Consistency
The internal metabolomics pipeline uses the same:
- Notion database (Experimental Data Assets)
- RAG embedding pipeline
- Metadata structure (as lipidomics)

### Differences
- **Data Origin**: "Internal – Amprenta" (same as lipidomics)
- **Source**: File path (same as lipidomics)
- **Format**: CSV/TSV (same as lipidomics)
- **Omics Type**: "Metabolomics" (vs no omics type for lipidomics)
- **No Signature Scoring**: Metabolomics datasets don't use lipid signature matching

### RAG Embedding
Both internal metabolomics and lipidomics datasets:
- Use same embedding pipeline
- Create similar text representations
- Upsert to Pinecone with proper metadata
- Update Embedding IDs and Last Embedded

## Future Enhancements

### Planned Features
- **Full file attachment**: Upload files to Notion (requires external storage)
- **Quantitative analysis**: Use intensity/abundance in future analysis
- **Batch processing**: Process multiple files in one run
- **Enhanced normalization**: Expand synonym mapping with external databases
- **Validation**: File format and schema validation
- **Program/Experiment relations**: Complete implementation when schema confirmed

### Requesting Features
If you need additional features or encounter issues:
1. Check existing documentation
2. Review logs for error details
3. Test with minimal example file
4. Report specific use case and expected behavior

