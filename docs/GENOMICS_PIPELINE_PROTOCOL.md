# Genomics Pipeline Protocol (Extended)

## Source
This protocol is from Gemini, acting as an expert Bioinformatics Data Engineer. It extends the basic ENA repository access with a complete pipeline for processing raw sequencing data into gene count matrices.

## Overview

This protocol defines the complete genomics pipeline:
1. **Search & Download**: Find and download FASTQ files from ENA
2. **Feature Extraction**: Quantify using Salmon/Kallisto
3. **Gene Count Matrix**: Generate quantified gene expression data

---

## A. The Search & Download Protocol

### Repository
- **Primary Repository:** ENA (European Nucleotide Archive)
- **API Endpoint:** `https://www.ebi.ac.uk/ena/portal/api/search`
- **Do NOT use:** SRA Toolkit (requires complex configuration)

### Critical Fields
When searching ENA, you must request:
- `fastq_ftp` - FTP download links
- `run_accession` - Run identifier  
- `sample_alias` - Sample identifier

### Download Method

1. **FTP to HTTP Conversion**: 
   - ENA API returns FTP links (e.g., `ftp.sra.ebi.ac.uk/...`)
   - Convert to HTTP: Prepend `http://` and use `requests.get(stream=True)`
   - Example: `ftp.sra.ebi.ac.uk/vol1/fastq/...` → `http://ftp.sra.ebi.ac.uk/vol1/fastq/...`

2. **File Size Warning**:
   - FASTQ files are HUGE (gigabytes to terabytes)
   - **Only download if user explicitly confirms**
   - **Alternative**: Download a subset (first 1000 lines) for testing

3. **Streaming Download**:
   - Use `requests.get(url, stream=True)` for large files
   - Process in chunks to avoid memory issues

---

## B. The Feature Extraction Protocol (Quantification)

### Problem Statement
ENA provides **raw sequencing data** (FASTQ files). You cannot "extract" a gene count matrix directly - you must **generate** it through quantification.

### Tools

**Primary:** **Salmon** (preferred)
- Fast pseudo-alignment tool
- Much faster than traditional aligners (STAR/Bowtie)
- Runs on standard laptops
- Optimized for transcript quantification

**Alternative:** **Kallisto**
- Another pseudo-alignment tool
- Similar performance to Salmon

### Why Pseudo-Alignment?

Traditional aligners (STAR, Bowtie, BWA):
- Require full alignment to reference genome
- Very slow (hours to days)
- Require powerful servers

Pseudo-alignment tools (Salmon, Kallisto):
- Fast approximate alignment
- Minutes to hours
- Run on standard laptops
- Perfect for transcript quantification

### Agent Restrictions

- **Cannot install binary tools**: You (the agent) cannot install Salmon/Kallisto
- **Use subprocess**: Write Python scripts that call tools via `subprocess`
- **Assumes pre-installed**: Tools must be installed in user's environment

### Quantification Process

1. **Requires Transcriptome Index**:
   - Pre-built index of reference transcripts
   - Built separately (not part of this pipeline)
   - Example: Human transcriptome index

2. **Salmon Command**:
   ```bash
   salmon quant -i index_path -l A -r fastq_file -o output_dir
   ```

3. **Output**:
   - `quant.sf` file with transcript counts
   - TPM (Transcripts Per Million) values
   - Can be aggregated to gene level

---

## Implementation Structure

### Repository Layer (Complete ✅)
- `amprenta_rag/ingestion/repositories/ena.py`
- Provides: Search, metadata, FASTQ link generation
- **Status**: Fully implemented and tested

### Pipeline Layer (New ✅)
- `amprenta_rag/ingestion/genomics/pipeline.py`
- Provides: Download, quantification, feature extraction
- **Status**: Implemented, requires Salmon/Kallisto installation

---

## Example Implementation

### Complete Pipeline Script

```python
from amprenta_rag.ingestion.genomics.pipeline import (
    get_ena_fastqs,
    download_fastq,
    quantify_with_salmon,
    extract_gene_counts_from_salmon,
)

# Step 1: Search for FASTQ files
runs = get_ena_fastqs("breast cancer", limit=1)

# Step 2: Download FASTQ (with confirmation)
fastq_path = download_fastq(
    run_info=runs[0],
    confirm=True,  # User must confirm
    subset=False,  # Full download (or True for subset)
)

# Step 3: Quantify with Salmon
quant_file = quantify_with_salmon(
    fastq_path=fastq_path,
    index_path=Path("./human_index"),
    output_dir=Path("./quants"),
)

# Step 4: Extract gene counts
gene_counts = extract_gene_counts_from_salmon(quant_file)
```

---

## Functions Provided

### 1. `get_ena_fastqs(keyword, limit, filters)`
- Searches ENA for FASTQ files
- Returns list of runs with HTTP-converted download URLs
- Filters out runs without FASTQ files

### 2. `download_fastq(run_info, output_dir, confirm, subset, max_lines)`
- Downloads FASTQ file from ENA
- Converts FTP to HTTP automatically
- Requires user confirmation (`confirm=True`)
- Supports subset download for testing

### 3. `quantify_with_salmon(fastq_path, index_path, output_dir, library_type)`
- Runs Salmon quantification via subprocess
- Generates `quant.sf` output file
- Requires pre-built Salmon index

### 4. `quantify_with_kallisto(fastq_path, index_path, output_dir)`
- Alternative quantification using Kallisto
- Generates `abundance.tsv` output file
- Requires pre-built Kallisto index

### 5. `extract_gene_counts_from_salmon(quant_file)`
- Reads Salmon `quant.sf` file
- Extracts transcript/gene counts (TPM values)
- Returns dictionary of gene counts

---

## Installation Requirements

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

---

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
| **Subset Download** | ✅ | For testing |
| **Salmon Integration** | ✅ | Via subprocess |
| **Kallisto Integration** | ✅ | Via subprocess |
| **Gene Count Extraction** | ✅ | From quant.sf |

---

## Files Created

- `amprenta_rag/ingestion/genomics/__init__.py` - Genomics module init
- `amprenta_rag/ingestion/genomics/pipeline.py` - Pipeline implementation
- `scripts/genomics_pipeline_example.py` - Example usage script
- `docs/GENOMICS_PIPELINE_PROTOCOL.md` - This documentation

---

## Usage Example

See `scripts/genomics_pipeline_example.py` for a complete example demonstrating:
1. Searching ENA
2. Downloading FASTQ files
3. Quantification (commented, requires Salmon)
4. Gene count extraction

---

## Notes

- **Large Files**: FASTQ files can be terabytes. Always require explicit confirmation.
- **Tool Installation**: Users must install Salmon/Kallisto separately.
- **Index Requirement**: Quantification requires pre-built transcriptome indices.
- **Processing Time**: Quantification can take hours for large datasets.
- **Memory Requirements**: Large FASTQ files require significant disk space.

---

## Next Steps

The pipeline is implemented and ready to use. Users need to:
1. Install Salmon or Kallisto
2. Build transcriptome indices
3. Run the pipeline with their FASTQ files

The repository layer (link generation) works immediately. The quantification pipeline requires external tools.
