# Automatic Linking Guide

**Last Updated**: 2025-12-07

This document describes the automatic linking system that intelligently connects datasets to programs and experiments based on metadata analysis and confidence scoring.

---

## Overview

The **Auto-Linking** system automatically infers program and experiment relationships for ingested datasets by analyzing metadata like disease context, keywords, filenames, matrix types, and model systems.

### Key Benefits

- **Reduced Manual Work**: Eliminates need to manually specify program/experiment IDs for every dataset
- **Intelligent Matching**: Uses confidence scoring to avoid ambiguous or low-quality matches
- **Safe by Default**: Only links when confidence threshold is met (default: 0.8 or 80%)
- **Non-Destructive**: Does not override explicit program/experiment IDs provided during ingestion
- **Audit Trail**: Logs confidence scores for all linking decisions

### When Auto-Linking Occurs

Auto-linking is invoked during:
- **Lipidomics ingestion**: `ingest_lipidomics_file()`
- **Metabolomics ingestion**: `ingest_metabolomics_file()`
- **Proteomics ingestion**: `ingest_proteomics_file()`
- **Transcriptomics ingestion**: `ingest_transcriptomics_file()`
- **Batch ingestion**: `batch_ingest_omics.py` (when `--program-id` not provided)

---

## How It Works

### Confidence Scoring Algorithm

Auto-linking uses a multi-factor scoring system to match datasets with programs/experiments.

#### Program Inference

**Scoring Factors** (`infer_program_from_metadata()`):

1. **Disease Overlap** (up to 0.7 points):
   - Base: 0.6 points for any disease match
   - Bonus: +0.1 per additional matching disease
   - Example: Dataset with ["ALS", "PD"] matching program with ["ALS", "PD", "AD"] → 0.8 points

2. **Keyword/Name Match** (0.3 points):
   - Matches keywords from:
     * Explicitly provided keywords
     * Filename tokens (split on `_` and `-`)
   - Example: Filename `ALS_CSF_study.csv` matches program name "ALS Research" → +0.3 points

**Total Score Range**: 0.0 to 1.0
- **0.8+ (default threshold)**: Auto-link to program
- **< 0.8**: Do not link (log warning with score)

**Ambiguity Handling**:
- If multiple programs have the same top score → no link (ambiguous)
- Ensures only confident, unambiguous matches are made

#### Experiment Inference

**Scoring Factors** (`infer_experiment_from_metadata()`):

1. **Disease Overlap** (up to 0.6 points):
   - Base: 0.5 points for any disease match
   - Bonus: +0.1 per additional matching disease

2. **Matrix Match** (0.25 points):
   - Matches sample matrix (CSF, Plasma, Serum, Tissue, etc.)
   - Example: Dataset with matrix=["CSF"] matches experiment with matrix=["CSF"] → +0.25 points

3. **Model System Match** (0.25 points):
   - Matches model systems (Mouse, Human, Cell Line, etc.)
   - Example: Dataset with model=["Mouse"] matches experiment with model=["Mouse"] → +0.25 points

**Total Score Range**: 0.0 to 1.1
- **0.8+ (default threshold)**: Auto-link to experiment
- **< 0.8**: Do not link

**Ambiguity Handling**:
- Same as program inference (no link if tied top scores)

---

## Configuration

Auto-linking is configured via environment variables in `.env`.

### Enable/Disable Auto-Linking

```bash
# Enable automatic linking (default: true)
AUTO_LINK_ENABLED=true

# Disable automatic linking (manual linking only)
AUTO_LINK_ENABLED=false
```

### Confidence Threshold

```bash
# Minimum confidence score for auto-linking (default: 0.8)
AUTO_LINK_MIN_CONFIDENCE=0.8

# More conservative (fewer auto-links, higher precision)
AUTO_LINK_MIN_CONFIDENCE=0.9

# More permissive (more auto-links, lower precision)
AUTO_LINK_MIN_CONFIDENCE=0.7
```

**Tuning Recommendations**:

| Threshold | Behavior | Use Case |
|-----------|----------|----------|
| 0.9 | Very conservative, only obvious matches | High-precision environments |
| 0.8 | Balanced (default) | Most use cases |
| 0.7 | More permissive, more auto-links | Well-structured metadata |
| 0.6 | Aggressive, many auto-links | Exploratory analysis |

---

## Usage Examples

### Example 1: Auto-Link During Single File Ingestion

```python
from amprenta_rag.ingestion.lipidomics_ingestion import ingest_lipidomics_file

# Ingest without explicit program_id
# Auto-linking will infer program based on metadata
page_id = ingest_lipidomics_file(
    file_path="data/ALS_CSF_lipidomics.csv",
    create_page=True
)

# Log output (if auto-link successful):
# [INGEST][AUTO-LINK] Linked dataset 2b9adf... to program abc123... (conf=0.85)
```

**What Happens**:
1. System extracts diseases from dataset (e.g., ["ALS"])
2. Searches Postgres for programs with matching diseases
3. Scores each program based on disease overlap + filename keywords
4. If top score ≥ 0.8 and unambiguous → links dataset to program

### Example 2: Auto-Link with Explicit Metadata

```python
from amprenta_rag.ingestion.metabolomics_ingestion import ingest_metabolomics_file

# Provide explicit metadata to improve auto-linking
page_id = ingest_metabolomics_file(
    file_path="data/metabolomics_study.csv",
    create_page=True,
    # These metadata fields are used for auto-linking
    diseases=["Alzheimer's Disease", "Dementia"],
    matrix=["CSF"],
    model_systems=["Human"]
)

# Auto-linking will use provided metadata for matching
```

**What Happens**:
1. System uses provided `diseases`, `matrix`, `model_systems` for scoring
2. Matches against both programs and experiments in Postgres
3. Links to program if confidence ≥ 0.8
4. Links to experiment if confidence ≥ 0.8 (independently scored)

### Example 3: Batch Ingestion with Auto-Linking

```bash
# Auto-linking enabled by default in batch mode
python scripts/batch_ingest_omics.py \
    --directory data/omics_files/ \
    --create-pages

# Each file auto-linked based on its metadata and filename
```

**Example Log Output**:
```text
[INGEST][AUTO-LINK] Linked dataset 2b9adf6142ab80... to program abc123... (conf=0.87)
[INGEST][AUTO-LINK] Linked dataset 2b9adf6142ab81... to experiment xyz789... (conf=0.92)
[INGEST][AUTO-LINK] No confident program match for dataset 2b9adf6142ab82... (best=0.65)
```

### Example 4: Override Auto-Linking with Explicit IDs

```python
# Explicit program_ids override auto-linking
page_id = ingest_lipidomics_file(
    file_path="data/dataset.csv",
    create_page=True,
    program_ids=["explicit-program-id-123"],  # Auto-linking skipped
    experiment_ids=["explicit-exp-id-456"]     # Auto-linking skipped
)

# No auto-linking occurs when explicit IDs are provided
```

### Example 5: Check Auto-Linking Confidence Programmatically

```python
from amprenta_rag.ingestion.auto_linking import (
    infer_program_from_metadata,
    infer_experiment_from_metadata
)

# Test program inference
program_id, confidence = infer_program_from_metadata(
    diseases=["ALS"],
    keywords=["CSF", "lipidomics"],
    filename="ALS_CSF_study.csv",
    min_confidence=0.8
)

if program_id:
    print(f"Would link to program {program_id} (confidence: {confidence:.2f})")
else:
    print(f"No confident match (best score: {confidence:.2f})")

# Test experiment inference
exp_id, confidence = infer_experiment_from_metadata(
    diseases=["ALS"],
    matrix=["CSF"],
    model_systems=["Human"],
    min_confidence=0.8
)

if exp_id:
    print(f"Would link to experiment {exp_id} (confidence: {confidence:.2f})")
else:
    print(f"No confident match (best score: {confidence:.2f})")
```

---

## Auto-Linking Scenarios

### Scenario 1: High Confidence Match

**Dataset**:
- Filename: `ALS_CSF_lipidomics_ST004396.csv`
- Diseases: ["ALS"]
- Matrix: ["CSF"]

**Program in Database**:
- Name: "ALS Research Program"
- Disease: ["ALS"]

**Scoring**:
- Disease overlap: 0.7 (base 0.6 + 0.1 for one match)
- Keyword match: 0.3 ("ALS" in filename matches "ALS" in program name)
- **Total: 1.0 → AUTO-LINK ✅**

### Scenario 2: Low Confidence (No Link)

**Dataset**:
- Filename: `study_123.csv`
- Diseases: ["Parkinson's Disease"]

**Program in Database**:
- Name: "Neurodegeneration Program"
- Disease: ["ALS", "AD", "HD"]

**Scoring**:
- Disease overlap: 0.0 (no matching diseases)
- Keyword match: 0.0 (no matching keywords)
- **Total: 0.0 → NO LINK ❌**

### Scenario 3: Ambiguous (No Link)

**Dataset**:
- Diseases: ["ALS"]

**Programs in Database**:
1. Program A: Disease ["ALS"] → score 0.7
2. Program B: Disease ["ALS"] → score 0.7

**Scoring**:
- Both programs score 0.7 (ambiguous tie)
- **Result: NO LINK ❌** (logs: "Ambiguous match")

### Scenario 4: Experiment Match

**Dataset**:
- Diseases: ["ALS"]
- Matrix: ["CSF"]
- Model Systems: ["Human"]

**Experiment in Database**:
- Disease: ["ALS"]
- Matrix: ["CSF"]
- Model Systems: ["Human"]

**Scoring**:
- Disease overlap: 0.6
- Matrix match: 0.25
- Model match: 0.25
- **Total: 1.1 → AUTO-LINK ✅**

---

## Best Practices

### 1. Use Descriptive Filenames

Good filenames improve auto-linking accuracy:

**Good**:
- `ALS_CSF_lipidomics_cohort1.csv`
- `PD_plasma_metabolomics_MTBLS123.csv`
- `AD_mouse_transcriptomics_GSE12345.csv`

**Bad** (poor for auto-linking):
- `data.csv`
- `study_v2.csv`
- `results_final.csv`

### 2. Populate Metadata in Postgres

Ensure programs and experiments have:
- **Diseases**: Primary field for matching
- **Program Names**: Include keywords for filename matching
- **Matrix** (experiments): CSF, Plasma, Serum, etc.
- **Model Systems** (experiments): Human, Mouse, etc.

### 3. Monitor Auto-Linking Logs

Check logs for auto-linking decisions:

```bash
# Search logs for auto-link events
grep "AUTO-LINK" logs/ingestion.log

# Look for low-confidence warnings
grep "No confident.*match" logs/ingestion.log
```

### 4. Tune Confidence Threshold

Start with default (0.8), adjust based on results:

- **Too many false positives** (wrong links) → increase to 0.9
- **Too few auto-links** (missing obvious matches) → decrease to 0.7

### 5. Explicitly Provide IDs When Known

Auto-linking is a convenience feature. When you know the correct program/experiment:

```bash
# Explicit IDs are always preferred
python scripts/ingest_lipidomics.py \
    --file data.csv \
    --program-id abc123 \
    --experiment-id xyz789
```

---

## Troubleshooting

### Issue: No Auto-Links Occurring

**Symptoms**: All datasets ingested without program/experiment links

**Solutions**:

1. **Check if auto-linking is enabled**:
   ```bash
   grep AUTO_LINK_ENABLED .env
   # Should be: AUTO_LINK_ENABLED=true
   ```

2. **Lower confidence threshold**:
   ```bash
   echo "AUTO_LINK_MIN_CONFIDENCE=0.7" >> .env
   ```

3. **Verify Postgres has programs/experiments**:
   ```python
   from amprenta_rag.database.base import get_db
   from amprenta_rag.database.models import Program
   
   db = next(get_db())
   programs = db.query(Program).all()
   print(f"Programs in database: {len(programs)}")
   for p in programs:
       print(f"  - {p.name}: {p.disease}")
   ```

### Issue: Wrong Programs Being Linked

**Symptoms**: Datasets linked to incorrect programs

**Solutions**:

1. **Increase confidence threshold**:
   ```bash
   echo "AUTO_LINK_MIN_CONFIDENCE=0.9" >> .env
   ```

2. **Check program metadata quality**:
   - Ensure programs have accurate disease fields
   - Avoid overlapping diseases across programs

3. **Review scoring logic**:
   ```python
   from amprenta_rag.ingestion.auto_linking import infer_program_from_metadata
   
   # Test with your dataset's metadata
   program_id, conf = infer_program_from_metadata(
       diseases=["YourDisease"],
       filename="your_file.csv",
       min_confidence=0.8
   )
   print(f"Matched: {program_id}, Score: {conf}")
   ```

### Issue: Ambiguous Matches

**Symptoms**: Log shows "Ambiguous match" or "tied top scores"

**Solutions**:

1. **Differentiate programs with unique keywords**:
   - Add distinguishing keywords to program names
   - Example: "ALS Program - Cohort A" vs "ALS Program - Cohort B"

2. **Add more metadata to datasets**:
   ```python
   # Provide more context for matching
   ingest_lipidomics_file(
       file_path="data.csv",
       diseases=["ALS"],
       keywords=["cohort_a", "early_stage"],  # Help disambiguate
       create_page=True
   )
   ```

### Issue: Auto-Linking Too Slow

**Symptoms**: Ingestion slower than expected

**Solutions**:

1. **Disable auto-linking for batch operations**:
   ```bash
   echo "AUTO_LINK_ENABLED=false" >> .env
   ```

2. **Link explicitly in batch**:
   ```bash
   # Faster: specify program_id once for all files
   python scripts/batch_ingest_omics.py \
       --directory data/ \
       --program-id abc123
   ```

---

## Implementation Details

### Source Code

Auto-linking logic is implemented in:
```
amprenta_rag/ingestion/auto_linking.py
```

**Key Functions**:
- `infer_program_from_metadata()`: Program inference with confidence scoring
- `infer_experiment_from_metadata()`: Experiment inference with confidence scoring

### Database Requirements

Auto-linking queries the **Postgres database** for programs and experiments:

**Programs Table** (`amprenta_rag.database.models.Program`):
- `id`: UUID
- `name`: Program name (used for keyword matching)
- `disease`: List of diseases (used for overlap scoring)

**Experiments Table** (`amprenta_rag.database.models.Experiment`):
- `id`: UUID
- `disease`: List of diseases
- `matrix`: List of sample matrices
- `model_systems`: List of model systems

### Integration Points

Auto-linking is called by:
1. `amprenta_rag/ingestion/omics_service.py` (central ingestion service)
2. All omics-specific ingestion modules (lipidomics, metabolomics, proteomics, transcriptomics)

### Thread Safety

Auto-linking uses database sessions safely:
- Each inference creates a new database session
- Sessions are closed in `finally` blocks
- Safe for parallel ingestion workers

---

## Advanced Configuration

### Custom Scoring Weights

To customize scoring weights, modify `auto_linking.py`:

```python
# In infer_program_from_metadata()

# Default weights:
# - Disease overlap: 0.6 base + 0.1 per extra match
# - Keyword match: 0.3

# Custom weights (example):
if overlap:
    score += 0.7 + 0.15 * len(overlap)  # Higher disease weight
if keywords_norm & name_tokens:
    score += 0.2  # Lower keyword weight
```

### Logging Auto-Link Decisions

Enable detailed auto-link logging:

```bash
# In .env
LOG_LEVEL=DEBUG

# Logs will show:
# - All candidate programs/experiments scored
# - Individual score components
# - Ambiguity details
```

### Dry-Run Mode (Testing)

Test auto-linking without making actual links:

```python
from amprenta_rag.ingestion.auto_linking import (
    infer_program_from_metadata,
    infer_experiment_from_metadata
)

# Test program inference
prog_id, prog_conf = infer_program_from_metadata(
    diseases=["ALS"],
    keywords=["CSF"],
    filename="ALS_study.csv"
)

# Test experiment inference
exp_id, exp_conf = infer_experiment_from_metadata(
    diseases=["ALS"],
    matrix=["CSF"],
    model_systems=["Human"]
)

print(f"Program: {prog_id} (confidence: {prog_conf:.2f})")
print(f"Experiment: {exp_id} (confidence: {exp_conf:.2f})")
```

---

## See Also

- [Configuration Guide](CONFIGURATION.md) - Configuration options
- [Usage Examples](USAGE_EXAMPLES.md) - Ingestion examples
- [Batch Ingestion Guide](USAGE_EXAMPLES.md#batch-ingestion) - Batch workflows
- [Database Models](../amprenta_rag/database/models.py) - Program/Experiment schemas

---

**Questions or Issues?**

See [Troubleshooting Guide](TROUBLESHOOTING.md) or file an issue on GitHub.

