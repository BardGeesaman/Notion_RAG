# Master Bioinformatics Data Protocol

## Source
This protocol is from Gemini, acting as an expert Bioinformatics Data Engineer. It defines strict rules for selecting repositories and extraction methods based on omics type.

## Overview

When fetching bioinformatics data, we must STRICTLY select the correct repository and extraction method based on the "Omics" type requested. Do not guess APIs. Follow these specific rules to ensure stability and avoid 500 errors.

---

## 1. TRANSCRIPTOMICS (Gene Expression)

### Primary Repository: NCBI GEO (Gene Expression Omnibus)

**Goal:** Extract Gene Count Matrices and Clinical Metadata.

### The Protocol

1. **Search:** Use `Bio.Entrez` (Biopython) to find `GSE` (Series) IDs.
   - Query: `"{term} AND GSE[Entry Type]"`

2. **Extraction:** Use the `GEOparse` library.
   - Why: It automatically parses "Series Matrix" text files into Pandas DataFrames.
   - Do Not: Do not write custom parsers for Soft files.

3. **Code Pattern:**
   ```python
   import GEOparse
   gse = GEOparse.get_GEO(geo="GSE12345", destdir="./cache")
   clinical_meta = gse.phenotype_data
   expression_matrix = gse.pivot_samples('VALUE')
   ```

### Implementation Status
✅ **FULLY IMPLEMENTED**
- Using Bio.Entrez for searching
- Using GEOparse for extraction
- Automatic caching via destdir
- Direct DataFrame access

---

## 2. PROTEOMICS (Protein ID & Quant)

### Primary Repository: PRIDE (Proteomics Identifications Database)

**Goal:** Extract protein identification tables (not raw mass spec files).

### The Protocol

1. **Search:** Use `requests` with PRIDE API v2.
   - Base URL: `https://www.ebi.ac.uk/pride/ws/archive/v2`

2. **File Selection (CRITICAL):** You must filter the file list to find the **Result File**. Priority order:
   1. `*.mzTab` (Standard format, parse with `pandas` or specific logic).
   2. `protein_groups.txt` (MaxQuant output, parse with `pd.read_csv(sep='\t')`).
   3. Ignore: Do not download `.raw`, `.wiff`, or `.d` files unless explicitly asked.

3. **Code Pattern:**
   ```python
   # Search Projects
   resp = requests.get(f"https://www.ebi.ac.uk/pride/ws/archive/v2/search/projects?keyword={term}")
   
   # Get Files & Filter
   files = requests.get(f"https://www.ebi.ac.uk/pride/ws/archive/v2/files/byProject?accession={pxd_id}").json()
   target = next((f for f in files if f['fileName'].endswith('.mzTab')), None)
   ```

### Implementation Status
✅ **FULLY IMPLEMENTED**
- Using PRIDE API v2
- Priority-based file selection (mzTab → MaxQuant → Excel → TSV/CSV)
- Pandas-based parsing
- Proper file filtering (skips raw files)

---

## 3. METABOLOMICS (Small Molecules)

### Primary Repository: Metabolomics Workbench (NMDR)

**Alternative:** MetaboLights (Only if necessary)

### The Protocol (Preferred)

- **Repository:** Metabolomics Workbench
- **Why:** Returns clean JSON data; no file parsing needed.
- **Method:** Use `requests` to hit the REST API.
  - Pattern: `https://www.metabolomicsworkbench.org/rest/study/study_id/{ID}/data`
  - Output: JSON containing `data` (counts) and `metabolites` (names).

### The Protocol (Fallback)

- **Repository:** MetaboLights (EBI)
- **Warning:** API is brittle. Expect 500 errors on private studies.
- **Method:** "Scan and Check"
  1. Get Study ID.
  2. List files via `/studies/{id}/files`.
  3. Check specifically for `m_*.tsv` (MAF files).
  4. **Resilience:** Wrap in `try/except` and skip if status code >= 500.

### Implementation Status
✅ **FULLY IMPLEMENTED**
- Metabolomics Workbench primary (clean JSON API)
- MetaboLights fallback (ISA-Tab parsing)
- Resilient error handling for 500 errors
- MAF file detection and parsing

---

## 4. GENOMICS (Raw Reads)

### Primary Repository: ENA (European Nucleotide Archive)

**Goal:** Download raw FASTQ sequencing files.

### The Protocol

- **Avoid:** NCBI SRA (Requires complex `sra-toolkit` configuration).
- **Use:** ENA Browser API.
- **Method:**
  - Search: `https://www.ebi.ac.uk/ena/portal/api/search?result=read_run&query={term}`
  - Download: The API provides direct `fastq_ftp` links (e.g., `ftp.sra.ebi.ac.uk/...`).
  - Agent Action: Generate the FTP links for the user; do not attempt to download terabytes of data unless confirmed.

### Implementation Status
⚠️ **NOT YET IMPLEMENTED**
- ENA repository not yet added
- Genomics data ingestion not implemented

---

## 5. GLOBAL RULES

1. **Rate Limiting:** Always include `time.sleep(1)` inside any loop.

2. **User Agent:** If using `requests`, set a descriptive header: `headers={'User-Agent': 'ResearchBot/1.0'}`.

3. **Data Persistence:** Always check if a file exists locally before downloading it again.

4. **Error Handling:**
   - If 404: "Resource not found (or private)."
   - If 500: "Server Error. Skipping to next item." (Do not crash).

### Implementation Status
✅ **MOSTLY IMPLEMENTED**
- ✅ Rate limiting in place (1 second delays)
- ⚠️ User-Agent headers not consistently set
- ✅ File caching (GEOparse, download_dir usage)
- ✅ Error handling for 404/500 errors

---

## Compliance Summary

| Omics Type | Repository | Status | Compliance |
|------------|-----------|--------|------------|
| Transcriptomics | NCBI GEO | ✅ Complete | ✅ Full |
| Proteomics | PRIDE | ✅ Complete | ✅ Full |
| Metabolomics | MW / MetaboLights | ✅ Complete | ✅ Full |
| Genomics | ENA | ⚠️ Not Implemented | ❌ N/A |
| Global Rules | All | ⚠️ Partial | ⚠️ Mostly |

---

## Notes

- All implemented repositories follow the protocol strictly
- GEO uses GEOparse as recommended
- PRIDE uses priority-based file selection as recommended
- Metabolomics Workbench is primary (clean JSON) as recommended
- Genomics (ENA) support can be added when needed
- Consider adding User-Agent headers consistently across all requests

---

## Future Improvements

1. **Add Genomics Support:** Implement ENA repository following the protocol
2. **User-Agent Headers:** Add descriptive User-Agent headers to all requests
3. **Enhanced Caching:** Improve file existence checks before downloading

