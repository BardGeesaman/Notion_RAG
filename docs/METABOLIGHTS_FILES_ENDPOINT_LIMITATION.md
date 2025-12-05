# MetaboLights Files Endpoint Limitation

## Issue

The MetaboLights API endpoint `/studies/{study_id}/files` consistently returns **500 Internal Server Error**.

## Diagnostic Results

✅ **Working Endpoints:**
- `/studies/{study_id}` - Returns study metadata (Status: 200)
- `/studies` - Lists all studies (Status: 200)

❌ **Broken Endpoint:**
- `/studies/{study_id}/files` - Returns 500 Internal Server Error

## Workaround

Files can still be accessed via:

1. **HTTP/FTP URLs from Study Details:**
   - Study details include `studyHttpUrl` and `studyFtpUrl`
   - Example: `http://ftp.ebi.ac.uk/pub/databases/metabolights/studies/public/MTBLS1`

2. **ISA-Tab File Structure:**
   MetaboLights studies follow ISA-Tab format with standard file naming:
   - `i_Investigation.txt` - Investigation metadata
   - `s_[study_id].txt` - Sample metadata
   - `a_[study_id]_*.txt` - Assay files
   - `m_[study_id]_*.tsv` - Metabolite data files (for feature extraction)

## Impact

- `fetch_study_data_files()` will return an empty list
- Feature extraction can work directly with ISA-Tab files
- File discovery must be done via ISA-Tab structure or directory listing

## Status

This is a **MetaboLights API server-side issue**, not a problem with our implementation. The protocol-compliant code correctly handles this limitation and logs appropriate warnings.

