# Metadata Editing Guide

**Last Updated**: December 7, 2025

This document describes how to view and edit metadata for datasets, programs, experiments, and signatures using the dashboard.

---

## Overview

The **Data Management** page provides a centralized interface for:

- **Editing metadata** - Update all editable fields for any entity
- **Linking entities** - Create relationships between datasets, programs, experiments, features, and signatures
- **Viewing raw metadata** - Inspect complete JSON metadata for debugging
- **Bulk operations** - Perform actions on multiple entities at once

**Access**: Dashboard → Data Management page

---

## Editable Fields by Entity Type

### Dataset Metadata

Datasets have the most comprehensive metadata schema:

| Field | Type | Description | Example |
|-------|------|-------------|---------|
| **Name** | Text | Dataset identifier/title | `ST004396_ALS_CSF_lipidomics` |
| **Description** | Long text | Detailed description | `Lipidomics analysis of CSF from 30 ALS patients vs 30 controls` |
| **Omics Type** | Select | Data type | `lipidomics`, `metabolomics`, `proteomics`, `transcriptomics`, `other` |
| **Disease** | Comma-separated | Disease contexts | `ALS, Neurodegeneration` |
| **Sample Type** | Comma-separated | Biological matrix | `CSF, Plasma, Serum` |
| **Organism** | Comma-separated | Species | `Human, Mouse` |
| **Methods** | Long text | Experimental methods | `Lipids extracted using Folch method...` |
| **Summary** | Long text | Brief summary | `We measured 245 lipid species in CSF samples` |
| **Results** | Long text | Key findings | `C16 and C18 ceramides elevated 3-fold in ALS` |
| **Conclusions** | Long text | Interpretations | `Ceramide dysregulation may contribute to pathogenesis` |
| **Data Origin** | Text | Source information | `Metabolomics Workbench`, `Internal` |
| **Dataset Source Type** | Text | Data type | `Public Repository`, `Internal Study`, `Collaborator` |

**Auto-populated fields** (not editable via dashboard):
- `id` - UUID assigned at creation
- `created_at` - Timestamp when dataset created
- `updated_at` - Timestamp of last modification
- `quality_score` - Computed by quality checks system
- `quality_status` - High/Medium/Low based on score
- `quality_issues` - List of quality problems

---

### Program Metadata

Programs represent research initiatives or projects:

| Field | Type | Description | Example |
|-------|------|-------------|---------|
| **Name** | Text | Program name | `ALS Research Program` |
| **Description** | Long text | Program overview | `Multi-year initiative to study ceramide metabolism in ALS` |
| **Disease** | Comma-separated | Disease focus | `ALS, Motor Neuron Disease` |
| **PI Name** | Text | Principal investigator | `Dr. Jane Smith` |
| **Institution** | Text | Organization | `University of California, San Francisco` |
| **Funding Source** | Text | Sponsor | `NIH, R01-NS123456` |
| **Start Date** | Date | Program start | `2023-01-01` |
| **End Date** | Date | Program end | `2027-12-31` |

---

### Experiment Metadata

Experiments are specific studies within programs:

| Field | Type | Description | Example |
|-------|------|-------------|---------|
| **Name** | Text | Experiment name | `ALS CSF Lipidomics - Cohort 1` |
| **Description** | Long text | Experiment details | `Cross-sectional lipidomics study of CSF from ALS patients` |
| **Disease** | Comma-separated | Disease context | `ALS` |
| **Sample Type** | Comma-separated | Matrix | `CSF` |
| **Organism** | Comma-separated | Species | `Human` |
| **Model System** | Comma-separated | Experimental system | `Human, Patient-derived` |
| **Methods** | Long text | Experimental methods | `CSF collected via lumbar puncture, stored at -80°C` |
| **Summary** | Long text | Brief overview | `30 ALS patients, 30 age-matched controls` |
| **Results** | Long text | Key findings | `Significant ceramide elevation in ALS group` |
| **Date** | Date | Experiment date | `2024-06-15` |

---

### Signature Metadata

Signatures are collections of features that characterize a biological state:

| Field | Type | Description | Example |
|-------|------|-------------|---------|
| **Name** | Text | Signature name | `ALS-CSF-Core-6Ceramides` |
| **Description** | Long text | Signature details | `Six ceramide species consistently elevated in ALS CSF` |
| **Version** | Text | Version identifier | `v1.0`, `v2.1` |
| **Signature Type** | Text | Classification | `Literature-derived`, `Consortium`, `Open Dataset` |
| **Modalities** | Comma-separated | Feature types | `Lipid, Metabolite, Protein` |
| **Confidence** | Number | Quality score | `0.92` (0-1 scale) |

---

## How to Edit Metadata

### Step-by-Step Workflow

#### 1. Access Data Management

```bash
# Start dashboard
cd scripts/dashboard
streamlit run app.py

# Navigate to Data Management page (sidebar)
```

#### 2. Select Edit Metadata Tab

Click the **"Edit Metadata"** tab (second tab)

#### 3. Choose Entity Type

Select from dropdown:
- Dataset
- Program
- Experiment
- Signature

#### 4. Select Specific Entity

Use the dropdown to choose which entity to edit

**Current Metadata Display**:
- Shows read-only current values
- Includes ID, created/updated timestamps
- Provides context before editing

#### 5. Edit Fields

**Text fields**:
- Click in field and type new value
- Clear field to remove value

**Text areas** (multi-line):
- Click in area and type
- Use Shift+Enter for new lines

**Comma-separated fields**:
- Enter values separated by commas
- Spaces around commas are trimmed automatically
- Example: `ALS, PD, AD` → `["ALS", "PD", "AD"]`

**Select fields**:
- Click dropdown and choose new value
- Example: Omics Type dropdown

#### 6. Save Changes

Click **"💾 Save Changes"** button

**Behavior**:
- Only changed fields are updated (efficient)
- "No changes detected" message if nothing modified
- Success message on successful save
- Page refreshes to show updated values

---

## Field Editing Examples

### Example 1: Add Disease Annotation

**Before**:
```
Disease: (empty)
```

**Edit**:
```
Disease: ALS, Amyotrophic Lateral Sclerosis, Motor Neuron Disease
```

**After save**:
- Dataset now has disease metadata
- Quality score improves (+5 points)
- Auto-linking can now match to ALS programs
- Filterable by disease in queries

---

### Example 2: Update Description

**Before**:
```
Description: (empty)
```

**Edit**:
```
Description: Lipidomics profiling of cerebrospinal fluid (CSF) from 30 patients 
with sporadic ALS and 30 age-matched healthy controls. Samples collected via 
lumbar puncture and analyzed using LC-MS/MS. A total of 245 lipid species 
quantified across major lipid classes including ceramides, sphingomyelins, 
phosphatidylcholines, and phosphatidylethanolamines.
```

**After save**:
- Improves quality score (+5 points)
- Enables semantic search of description
- Provides context for RAG queries
- Helps future users understand dataset

---

### Example 3: Correct Omics Type

**Before**:
```
Omics Type: other
```

**Edit**:
```
Omics Type: lipidomics (select from dropdown)
```

**After save**:
- Dataset correctly categorized
- Appears in lipidomics filters
- Uses appropriate pipelines for analysis
- Quality checks use correct thresholds

---

### Example 4: Add Methodology

**Before**:
```
Methods: (empty)
```

**Edit**:
```
Methods: Lipids were extracted from 200 µL CSF using the Folch method with 
chloroform:methanol (2:1, v/v). Internal standards (C17:0 ceramide, d7-SM) 
were added before extraction. Lipid extracts were analyzed by LC-MS/MS using 
an Agilent 6490 triple quadrupole mass spectrometer in positive ion mode with 
multiple reaction monitoring (MRM). Data were normalized to total lipid 
content and sample volume.
```

**After save**:
- Provides reproducibility information
- Helps assess data quality
- Enables comparison with other studies
- Important for meta-analysis

---

### Example 5: Link Multiple Diseases

**Before**:
```
Disease: ALS
```

**Edit**:
```
Disease: ALS, Motor Neuron Disease, Neurodegeneration, Sporadic ALS
```

**After save**:
- Richer metadata for filtering
- Matches broader search queries
- Better auto-linking to programs
- More findable in semantic search

---

## Raw Metadata Viewer

### Purpose

The raw metadata viewer shows the complete JSON representation of an entity's metadata, including fields not editable via the UI.

**Use Cases**:
- Debugging data issues
- Verifying imports from repositories
- Checking computed fields (quality scores, relationships)
- Understanding data structure for programmatic access

### How to Access

1. Navigate to Edit Metadata tab
2. Select entity type and specific entity
3. Scroll to bottom of page
4. Click **"View raw metadata"** expander

### What's Included

**All fields in database**:
- Editable fields (name, description, etc.)
- Auto-computed fields (quality scores, timestamps)
- Relationships (linked programs, experiments, features)
- Repository metadata (data_origin, accession IDs)
- Internal fields (UUIDs, database keys)

**Example Raw Metadata (Dataset)**:
```json
{
  "id": "2b9adf6142ab8026853ef58f725665a6",
  "name": "ST004396_ALS_CSF_lipidomics",
  "description": "Lipidomics analysis of CSF from ALS patients",
  "omics_type": "lipidomics",
  "disease": ["ALS", "Neurodegeneration"],
  "sample_type": ["CSF"],
  "organism": ["Human"],
  "methods": "Folch extraction, LC-MS/MS...",
  "summary": "245 lipid species quantified",
  "results": "C16/C18 ceramides elevated 3-fold",
  "conclusions": "Ceramide dysregulation in ALS",
  "data_origin": "Metabolomics Workbench",
  "dataset_source_type": "Public Repository",
  "quality_score": 95.0,
  "quality_status": "high",
  "quality_issues": [],
  "created_at": "2024-12-01T10:30:00Z",
  "updated_at": "2024-12-07T15:45:00Z",
  "programs": [
    {
      "id": "abc-123-def",
      "name": "ALS Research Program"
    }
  ],
  "experiments": [
    {
      "id": "xyz-789-uvw",
      "name": "ALS CSF Lipidomics - Cohort 1"
    }
  ],
  "features": [
    {"id": "feat-001", "name": "Cer(d18:1/16:0)", "type": "lipid"},
    {"id": "feat-002", "name": "Cer(d18:1/18:0)", "type": "lipid"}
  ]
}
```

---

## Study Grouping (Repository Imports)

### What is Study Grouping?

When importing datasets from public repositories (Metabolomics Workbench, GEO, PRIDE, MetaboLights), the system automatically creates **Experiment** entities to group related datasets from the same study.

**Why It Matters**:
- Multi-dataset studies stay organized
- Relationships preserved from source repository
- Easy to find all datasets from one study
- Maintains experimental context

### How It Works

#### Step 1: Repository Import

When you import a study (e.g., `ST004396` from Metabolomics Workbench):

```bash
python scripts/ingest_from_repository.py \
    --repository MW \
    --study-id ST004396 \
    --create-experiment
```

#### Step 2: Experiment Auto-Creation

System creates an Experiment with metadata from repository:

**Experiment Metadata**:
- **Name**: Study title from repository
- **Description**: Study abstract/summary
- **Disease**: Extracted from study metadata
- **Sample Type**: Extracted from methods
- **Organism**: From sample information
- **Methods**: Study protocols
- **Date**: Study submission date

**Example**:
```
Experiment created:
  Name: "Lipidomics analysis of CSF from ALS patients (ST004396)"
  Description: "Cross-sectional study comparing lipid profiles..."
  Disease: ["ALS"]
  Sample Type: ["CSF"]
  Organism: ["Human"]
```

#### Step 3: Dataset-Experiment Linking

All datasets from the study are automatically linked to the experiment:

**What Gets Linked**:
- Dataset 1: ALS patient group → Experiment
- Dataset 2: Control group → Experiment
- Dataset 3: Validation cohort → Experiment

**Result**: All three datasets now grouped under one experiment

### Multi-Dataset Study Example

**Study**: ST004396 (Metabolomics Workbench)

**Imported Datasets**:
1. `ST004396_ALS_patients` - 30 ALS CSF samples
2. `ST004396_controls` - 30 healthy control CSF samples
3. `ST004396_validation` - 15 independent ALS samples

**Auto-Created Experiment**:
```
Name: ST004396 - ALS CSF Lipidomics
Description: Multi-cohort lipidomics study of ALS cerebrospinal fluid
Disease: [ALS]
Linked Datasets: 3
  - ST004396_ALS_patients (n=30)
  - ST004396_controls (n=30)
  - ST004396_validation (n=15)
```

**Benefits**:
- Easily find all datasets from this study
- Understand dataset relationships (case vs control)
- Compare results across cohorts
- Maintain study context for interpretation

---

## Linking Entities

The Data Management page also provides tools to create relationships between entities.

### Available Link Types

1. **Dataset → Program**
   - Link datasets to research programs
   - Datasets can belong to multiple programs

2. **Dataset → Experiment**
   - Link datasets to experiments
   - Datasets can belong to multiple experiments

3. **Feature → Dataset**
   - Link features (genes, proteins, metabolites, lipids) to datasets
   - Automated during ingestion, manual linking for corrections

4. **Signature → Dataset**
   - Link signatures to datasets where they appear
   - Used for signature scoring and matching

### How to Link Entities

#### Access Link Entities Tab

Click **"Link Entities"** tab (first tab in Data Management)

#### Example: Link Dataset to Program

**Scenario**: You imported a dataset from a repository and want to link it to your research program.

**Steps**:
1. Select link type: "Dataset → Program"
2. Choose dataset from dropdown
3. Select one or more programs to link (multiselect)
4. Click "🔗 Link Dataset to Programs"

**Result**:
- Dataset now associated with program(s)
- Dataset appears in program's dataset list
- Can filter/query by program
- Cross-omics reasoning includes this dataset when analyzing program

---

## Best Practices

### 1. Fill Required Fields First

**Priority order**:
1. **Name** - Always required, must be descriptive
2. **Omics Type** - Essential for correct processing
3. **Description** - Improves findability and quality score
4. **Disease** - Enables filtering and auto-linking

**Example workflow**:
- Import dataset from repository
- Check if name is descriptive (edit if generic)
- Verify omics type is correct
- Add description if missing
- Add disease annotation

### 2. Use Consistent Terminology

**Disease names**:
- Use standard disease names: "ALS" not "Lou Gehrig's disease"
- Include synonyms for better search: "ALS, Amyotrophic Lateral Sclerosis"
- Use abbreviations consistently

**Sample types**:
- Standardize: "CSF" not "cerebrospinal fluid", "spinal fluid", "CSF fluid"
- Common terms: CSF, Plasma, Serum, Tissue, Cell Line

**Organism**:
- Use standard names: "Human", "Mouse", "Rat"
- Include strain if relevant: "Mouse (C57BL/6)"

### 3. Provide Context in Long Text Fields

**Methods**:
- Include extraction protocol
- List instruments used
- Specify normalization methods
- Note any deviations from standard protocols

**Results**:
- Summarize key findings
- Include statistics (p-values, fold-changes)
- Mention number of significant features

**Conclusions**:
- Interpret biological significance
- Connect to existing literature
- Note limitations or caveats

### 4. Keep Metadata Current

**When to update**:
- After re-ingestion or data correction
- When new analysis reveals issues
- After quality checks identify problems
- When linking to new programs/experiments

**What to update**:
- Quality issues resolved → update description
- New disease associations discovered → add to disease field
- Methods clarified → update methods section

### 5. Use Raw Metadata for Verification

**Check after**:
- Repository imports
- Bulk operations
- Auto-linking operations
- Data migrations

**Verify**:
- Relationships are correct (programs, experiments, features)
- Computed fields make sense (quality scores)
- Repository metadata preserved (accession IDs)

---

## Troubleshooting

### Issue: Changes Not Saving

**Problem**: Click "Save Changes" but metadata doesn't update

**Solutions**:
1. **Check for error message**: Red error text may indicate validation issue
2. **Verify field format**: Ensure comma-separated fields use commas, not semicolons
3. **Refresh page**: Browser cache may be stale
4. **Check database connection**: Dashboard logs may show connection errors

**Diagnostic**:
```bash
# Check dashboard logs for errors
# Terminal running streamlit run app.py will show error messages
```

---

### Issue: Can't Find Entity to Edit

**Problem**: Entity doesn't appear in dropdown

**Solutions**:
1. **Verify entity exists**: Check other dashboard pages (Datasets, Programs, etc.)
2. **Check entity type**: May be looking in wrong category (dataset vs experiment)
3. **Refresh data**: Click refresh button or reload page
4. **Check database**: Entity may not be ingested yet

---

### Issue: Quality Score Doesn't Update After Edits

**Problem**: Added metadata but quality score stays the same

**Explanation**: Quality scores are not automatically recomputed on metadata edit

**Solution**:
```bash
# Recompute quality scores for all datasets
python scripts/check_dataset_quality.py

# Refresh dashboard to see updated scores
```

---

### Issue: Can't Link to Program/Experiment

**Problem**: Link button doesn't work or shows error

**Solutions**:
1. **Verify entities exist**: Both source and target entities must exist
2. **Check for existing link**: May already be linked (check raw metadata)
3. **Refresh page**: Reload to get fresh entity lists
4. **Check permissions**: Database connection must have write access

---

## Programmatic Metadata Editing

For bulk edits or automation, use the Python API:

```python
from amprenta_rag.api.services.datasets import update_dataset
from amprenta_rag.api.schemas import DatasetUpdate
from amprenta_rag.database.base import get_db

db = next(get_db())

# Update dataset metadata
update_data = DatasetUpdate(
    description="New description text",
    disease=["ALS", "Neurodegeneration"],
    sample_type=["CSF"],
)

updated_dataset = update_dataset(
    db=db,
    dataset_id="your-dataset-uuid",
    dataset=update_data
)

db.commit()
print(f"Updated: {updated_dataset.name}")
```

---

## See Also

- [Quality Checks Guide](QUALITY_CHECKS.md) - Improve metadata quality
- [Auto-Linking Guide](AUTO_LINKING.md) - Automatic program/experiment linking
- [Usage Examples](USAGE_EXAMPLES.md) - Data ingestion workflows
- [API Reference](API_REFERENCE.md) - Programmatic metadata access

---

**Questions or Issues?**

See [Troubleshooting Guide](TROUBLESHOOTING.md) or file an issue on GitHub.

