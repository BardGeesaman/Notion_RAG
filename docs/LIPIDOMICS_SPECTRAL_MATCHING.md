# Lipidomics Spectral Matching

## Overview

The lipidomics spectral matching workflow in Amprenta links an internal `Feature` (a detected lipid feature) to the best matching reference spectra from one or more ingested spectral libraries (MGF).

High-level flow:

1. **Ingest reference library**: Load an MGF file into the database as a `SpectralLibrary` with many `SpectralReference` entries (each contains precursor m/z and fragment peaks).
2. **Populate feature spectrum**: Ensure each lipid `Feature` includes an MS/MS spectrum in `Feature.external_ids["spectrum"]` (see convention below).
3. **Run matching**: For a given `feature_id`, the service:
   - Filters reference spectra by **precursor m/z tolerance (ppm)**.
   - Computes **cosine similarity** between query and reference fragment peaks.
   - Stores the top matches as `LipidAnnotation` rows.
4. **Review/accept**: Users can review matches and (optionally) accept the top hit to rename the `Feature`.

## Feature Spectrum Convention

For spectral matching to work, `Feature.external_ids` must contain:

```json
{
  "spectrum": {
    "precursor_mz": 760.585,
    "peaks": [[184.07, 100.0], [264.27, 45.2], [520.51, 30.1]]
  }
}
```

- `precursor_mz`: Precursor ion m/z in Daltons.
- `peaks`: List of `[m/z, intensity]` pairs. Intensity can be arbitrary scale (relative intensity is sufficient).

Notes:

- The matcher expects `Feature.external_ids` to be a dict, with a nested dict at key `"spectrum"`.
- Empty/invalid peak lists will cause matching to fail with a 400 error from the API.

## Example: Populating during lipidomics ingestion

### Populate a Feature with spectrum data

```python
from amprenta_rag.database.models import Feature

# feature: Feature (SQLAlchemy model instance)
precursor_mz = 760.585
peaks = [
    [184.07, 100.0],
    [264.27, 45.2],
    [520.51, 30.1],
]

ext = feature.external_ids or {}
ext.update(
    {
        "spectrum": {
            "precursor_mz": precursor_mz,
            "peaks": peaks,
        }
    }
)
feature.external_ids = ext
```

### Example mapping from MGF parsing output

If you are parsing an MGF (e.g., via `amprenta_rag.spectral.library_parser.parse_mgf()`), the parser yields dicts like:

```python
{
  "precursor_mz": 760.585,
  "name": "PC 34:1",
  "peaks": [(184.07, 100.0), (264.27, 45.2), ...]
}
```

To store this into `Feature.external_ids["spectrum"]`, convert peaks to JSON-friendly lists:

```python
parsed = spectrum_dict_from_mgf

feature.external_ids = feature.external_ids or {}
feature.external_ids["spectrum"] = {
    "precursor_mz": float(parsed["precursor_mz"]),
    "peaks": [[float(mz), float(inten)] for mz, inten in parsed["peaks"]],
}
```

## Confidence Thresholds

The spectral matching service assigns flags using:

- **Confident**: cosine \(>= 0.7\) AND `mz_error_ppm` \(<= 10\)
- **Ambiguous**: multiple confident hits within **0.05** cosine score of each other (based on top-2 confident hits)

These values are encoded in `amprenta_rag/spectral/matching_service.py` and stored on each `LipidAnnotation` row as:

- `is_confident`
- `is_ambiguous`

## API Usage Examples

### 1) Ingest a spectral library (MGF)

```bash
curl -X POST "http://localhost:8000/api/spectral/libraries/ingest" \
  -H "Content-Type: application/json" \
  -d '{
    "mgf_path": "/absolute/path/to/library.mgf",
    "name": "LipidBlast",
    "version": "v1"
  }'
```

### 2) Run matching for a Feature

```bash
curl -X POST "http://localhost:8000/api/spectral/match/<feature_id>"
```

### 3) Fetch annotations (match results) for a Feature

```bash
curl -X GET "http://localhost:8000/api/spectral/annotations/<feature_id>"
```


