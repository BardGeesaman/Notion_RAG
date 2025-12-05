# Repository Ingestion Test Plan

**Date**: 2025-12-04  
**Status**: Ready for Testing

## Test Overview

Test the complete end-to-end workflow for importing a study from a public repository into the system.

## Test Scenario

**Repository**: Metabolomics Workbench (MW)  
**Study ID**: ST004168 (Alzheimer metabolomics study)  
**Reason**: Small study, already discovered, metabolomics type

## Test Steps

### Step 1: Verify Discovery ✅ (Already Done)
```bash
python scripts/discover_omics_studies.py \
  --repository MW \
  --keywords "Alzheimer" \
  --max-results 1
```
**Expected**: Study ST004168 found

### Step 2: Test Dry Run ✅ (Already Done)
```bash
python scripts/harvest_repository_study.py \
  --study-id ST004168 \
  --repository MW \
  --dry-run
```
**Expected**: Shows metadata without creating anything

### Step 3: Test Harvest (Create Notion Page Only)
```bash
python scripts/harvest_repository_study.py \
  --study-id ST004168 \
  --repository MW \
  --create-notion
```
**Expected**: 
- Creates Notion Dataset page
- Populates metadata
- Returns page ID

### Step 4: Test Full Import (Harvest + Postgres + Ingestion)
```bash
python scripts/harvest_repository_study.py \
  --study-id ST004168 \
  --repository MW \
  --create-notion \
  --ingest
```
**Expected**:
- Creates Postgres dataset
- Creates Notion page
- Downloads data files (if available)
- Extracts features
- Links features to Postgres
- Embeds in Pinecone
- Runs signature matching

### Step 5: Verify in Dashboard
- Open dashboard
- Check Datasets page
- Verify study appears
- Check linked features
- Verify metadata

### Step 6: Verify in Postgres
```sql
SELECT * FROM datasets WHERE external_ids->>'mw_study_id' = 'ST004168';
SELECT COUNT(*) FROM dataset_feature WHERE dataset_id = (SELECT id FROM datasets WHERE external_ids->>'mw_study_id' = 'ST004168');
```

## What to Check

### Postgres
- [ ] Dataset created in `datasets` table
- [ ] Study ID stored in `external_ids.mw_study_id`
- [ ] Features created and linked
- [ ] Relationships exist

### Notion
- [ ] Dataset page created
- [ ] Metadata populated correctly
- [ ] Summary contains study ID

### Features
- [ ] Features extracted from data
- [ ] Features linked to dataset in Postgres
- [ ] Feature normalization working

### Embedding
- [ ] Chunks created in Pinecone
- [ ] Metadata correct
- [ ] Embedding IDs stored in Notion

## Test Script

Quick test command:
```bash
# Full test
python scripts/harvest_repository_study.py \
  --study-id ST004168 \
  --repository MW \
  --create-notion \
  --ingest
```

## Rollback Plan

If test fails or creates unwanted data:
1. Delete Notion page (if created)
2. Delete Postgres dataset (if created)
3. Delete Pinecone vectors (if created)

## Expected Results

**Success Criteria**:
- ✅ Postgres dataset created with correct metadata
- ✅ Notion page created and linked
- ✅ Features extracted and linked
- ✅ Data visible in dashboard
- ✅ All relationships working

**Potential Issues**:
- Data file download might fail (network, file size)
- Feature extraction might need adjustment
- Signature matching might need tuning

