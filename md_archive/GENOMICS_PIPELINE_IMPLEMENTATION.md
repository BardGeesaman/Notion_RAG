# Genomics Pipeline Protocol Implementation

## Summary

Successfully implemented the Genomics Pipeline Protocol from Gemini, extending the basic ENA repository access with complete FASTQ download and quantification capabilities.

## Implementation Status

### ✅ Completed

1. **Genomics Pipeline Module**
   - Created: `amprenta_rag/ingestion/genomics/pipeline.py`
   - Functions: Search, download, quantification, gene count extraction

2. **FTP to HTTP Conversion**
   - Automatic conversion for `requests` library compatibility
   - Pipeline Protocol compliant

3. **Download Functionality**
   - Streaming downloads for large files
   - Subset download option (for testing - first N lines)
   - User confirmation required (Master Protocol)

4. **Quantification Support**
   - Salmon integration (via subprocess)
   - Kallisto integration (via subprocess)
   - Gene count extraction from quantification output

5. **Example Script**
   - `scripts/genomics_pipeline_example.py`
   - Complete workflow demonstration

### Integration

- Uses ENA repository for reliable searching
- Leverages existing repository infrastructure
- Maintains Master Protocol compliance

## Pipeline Functions

### 1. `get_ena_fastqs(keyword, limit, filters)`

**Purpose**: Find FASTQ download links for a study/keyword.

**Features**:
- Uses ENA repository for reliable searching
- Filters runs with no FASTQ files
- Returns HTTP-converted URLs (Pipeline Protocol)

**Example**:
```python
runs = get_ena_fastqs("Homo sapiens", limit=3)
# Returns: [{"Run": "DRR054724", "Sample": "...", "URL": "http://...", ...}]
```

### 2. `download_fastq(run_info, output_dir, confirm, subset, max_lines)`

**Purpose**: Download FASTQ file from ENA.

**Features**:
- Converts FTP to HTTP automatically
- Streaming download for large files
- Requires user confirmation (`confirm=True`)
- Supports subset download (first N lines) for testing

**Example**:
```python
fastq_path = download_fastq(
    run_info=runs[0],
    confirm=True,  # Required
    subset=True,  # Test with first 1000 lines
    max_lines=1000,
)
```

### 3. `quantify_with_salmon(fastq_path, index_path, output_dir, library_type)`

**Purpose**: Run Salmon quantification on FASTQ file.

**Requirements**:
- Salmon must be installed: `conda install -c bioconda salmon`
- Pre-built transcriptome index required

**Output**: `quant.sf` file with transcript counts (TPM)

**Example**:
```python
quant_file = quantify_with_salmon(
    fastq_path=Path("./fastq/DRR054724.fastq.gz"),
    index_path=Path("./human_index"),
    output_dir=Path("./quants/DRR054724"),
)
```

### 4. `quantify_with_kallisto(fastq_path, index_path, output_dir)`

**Purpose**: Alternative quantification using Kallisto.

**Requirements**:
- Kallisto must be installed: `conda install -c bioconda kallisto`
- Pre-built transcriptome index required

**Output**: `abundance.tsv` file with transcript counts

### 5. `extract_gene_counts_from_salmon(quant_file)`

**Purpose**: Extract gene counts from Salmon quantification output.

**Returns**: Dictionary mapping transcript/gene IDs to TPM values

**Example**:
```python
gene_counts = extract_gene_counts_from_salmon(quant_file)
# Returns: {"ENST00000123456": 1250.5, "ENST00000789012": 890.3, ...}
```

## Protocol Compliance

### Master Protocol Rules

| Rule | Status | Implementation |
|------|--------|----------------|
| **Rate Limiting** | ✅ | 1 second between requests |
| **User-Agent** | ✅ | Included in all requests |
| **Error Handling** | ✅ | 404/500 handled gracefully |
| **No Auto-Download** | ✅ | Requires `confirm=True` |

### Pipeline Protocol Features

| Feature | Status | Notes |
|---------|--------|-------|
| **FTP → HTTP Conversion** | ✅ | Automatic conversion |
| **Streaming Download** | ✅ | Chunked processing |
| **Subset Download** | ✅ | For testing (first N lines) |
| **Salmon Integration** | ✅ | Via subprocess |
| **Kallisto Integration** | ✅ | Via subprocess |
| **Gene Count Extraction** | ✅ | From quant.sf |

## Workflow

```
1. Search ENA → Find FASTQ files
   ↓
2. Download FASTQ → Convert FTP to HTTP, stream download
   ↓
3. Quantify → Run Salmon/Kallisto (requires index)
   ↓
4. Extract Gene Counts → Parse quant.sf file
   ↓
5. Link to Postgres → Store gene counts as features
```

## Files Created

- `amprenta_rag/ingestion/genomics/__init__.py` - Module initialization
- `amprenta_rag/ingestion/genomics/pipeline.py` - Pipeline implementation
- `scripts/genomics_pipeline_example.py` - Example usage
- `docs/GENOMICS_PIPELINE_PROTOCOL.md` - Protocol documentation
- `docs/GENOMICS_PIPELINE_IMPLEMENTATION.md` - This document

## Requirements

### For Quantification

**Salmon**:
```bash
conda install -c bioconda salmon
```

**Kallisto**:
```bash
conda install -c bioconda kallisto
```

### Index Building

Users must build transcriptome indices separately:

**Salmon Index**:
```bash
salmon index -t transcripts.fa -i salmon_index
```

**Kallisto Index**:
```bash
kallisto index -i kallisto_index transcripts.fa
```

## Test Results

✅ **Pipeline Functions**: All working correctly
✅ **Search Integration**: Uses ENA repository successfully
✅ **HTTP Conversion**: FTP links converted correctly
✅ **Code Structure**: Clean, modular, well-documented

## Notes

- **Large Files**: FASTQ files can be terabytes. Always require explicit confirmation.
- **Tool Installation**: Users must install Salmon/Kallisto separately.
- **Index Requirement**: Quantification requires pre-built transcriptome indices.
- **Processing Time**: Quantification can take hours for large datasets.

## Next Steps

The pipeline is implemented and ready to use. Users need to:
1. Install Salmon or Kallisto
2. Build transcriptome indices
3. Run the pipeline with their FASTQ files

The repository layer (link generation) works immediately. The quantification pipeline requires external tools.

