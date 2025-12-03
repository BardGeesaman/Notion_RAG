# Signature Match Score Implementation

## Overview

This document describes the implementation of automatic "Signature Match Score" writing to the Experimental Data Assets database in Notion during dataset ingestion.

## Feature Description

The system now automatically writes the highest signature match score to a numeric property called "Signature Match Score" on dataset pages when signatures are matched during ingestion.

## Implementation Details

### Location
- **File**: `amprenta_rag/ingestion/signature_matching.py`
- **Function**: `update_dataset_with_signature_matches()`

### Changes Made

1. **Added Signature Match Score Writing**
   - Computes the highest match score among all matched signatures
   - Writes the score to the "Signature Match Score" numeric property
   - Only writes when matches exist (idempotent behavior)

2. **Error Handling**
   - Gracefully handles cases where the property doesn't exist
   - Logs warnings but doesn't break ingestion
   - Retries update without score property if initial write fails

3. **Logging**
   - Clear logging when score is written
   - Logging when no matches are found
   - Includes score value in success messages

### Code Changes

```python
# Write Signature Match Score (numeric property)
# Only write if we have matches (idempotent: don't clear if no matches)
if matches and highest_score > 0:
    properties["Signature Match Score"] = {
        "number": highest_score,
    }
    logger.info(
        "[INGEST][SIGNATURE-MATCH] Writing Signature Match Score = %.3f to dataset %s",
        highest_score,
        dataset_page_id,
    )
```

### Error Handling

The implementation includes robust error handling:

```python
if resp.status_code >= 300:
    # Check if error is due to missing property
    error_text = resp.text
    if "Signature Match Score" in error_text and "not a property" in error_text:
        logger.warning(
            "[INGEST][SIGNATURE-MATCH] Signature Match Score property not found on dataset %s. "
            "Property may not exist in database schema. Continuing without score update.",
            dataset_page_id,
        )
        # Retry without the Signature Match Score property
        # ... retry logic ...
```

## Behavior

### When Matches Are Found
- Computes highest score: `max(match.score for match in matches)`
- Writes score to "Signature Match Score" property
- Logs: `[INGEST][SIGNATURE-MATCH] Writing Signature Match Score = X.XXX to dataset <page_id>`
- Updates relations and summary as before

### When No Matches Are Found
- Does NOT modify "Signature Match Score" property (idempotent)
- Logs: `[INGEST][SIGNATURE-MATCH] No signature matches found for dataset <page_id> — Signature Match Score unchanged`
- Returns early without updating Notion

### Property Requirements

- **Property Name**: `Signature Match Score`
- **Property Type**: Number (numeric)
- **Database**: Experimental Data Assets
- **Behavior**: Optional - ingestion continues even if property doesn't exist

## Testing

### Test Datasets

1. **ST004217** (Ceramide Accumulation Study)
   - Page ID: `2beadf61-42ab-81b7-b681-cf9fefb348d0`
   - Expected Score: 0.650
   - Result: ✅ Successfully written

2. **ST004223** (Sulfatide/Sphingolipid Study)
   - Page ID: `2beadf61-42ab-81a3-b07b-d60f5ee254a4`
   - Expected Score: 0.543
   - Result: ✅ Successfully written

### Test Commands

```bash
# Test ST004217
python scripts/ingest_dataset.py \
  --dataset-page-id 2beadf61-42ab-81b7-b681-cf9fefb348d0 \
  --force

# Test ST004223
python scripts/ingest_dataset.py \
  --dataset-page-id 2beadf61-42ab-81a3-b07b-d60f5ee254a4 \
  --force
```

### Expected Log Output

**When matches are found:**
```
[INGEST][SIGNATURE-MATCH] Writing Signature Match Score = 0.650 to dataset 2beadf61-42ab-81b7-b681-cf9fefb348d0
[INGEST][SIGNATURE-MATCH] Updated dataset 2beadf61-42ab-81b7-b681-cf9fefb348d0 with 2 signature match(es) (Signature Match Score = 0.650)
```

**When no matches are found:**
```
[INGEST][SIGNATURE-MATCH] No signature matches found for dataset <page_id> — Signature Match Score unchanged
```

## Integration Points

### Dataset Ingestion Flow

1. Dataset is ingested via `ingest_dataset()` in `dataset_ingestion.py`
2. mwTab data is extracted
3. Species are extracted from mwTab
4. Signature matching is triggered: `find_matching_signatures_for_dataset()`
5. Matches are written back: `update_dataset_with_signature_matches()`
   - **NEW**: Signature Match Score is written here

### Existing Functionality Preserved

- ✅ Relations to matching signatures still created
- ✅ Summary field still updated with match details
- ✅ Other metadata fields (Methods, Results, etc.) unchanged
- ✅ All existing logging preserved

## Configuration

No new configuration required. The feature uses existing:
- `cfg.notion.base_url` - Notion API base URL
- `cfg.notion.exp_data_db_id` - Experimental Data Assets database ID
- Signature matching threshold (default: 0.3)

## Error Scenarios

### Property Doesn't Exist
- **Behavior**: Logs warning, retries update without score property
- **Impact**: Ingestion continues successfully
- **Log**: Warning message about missing property

### API Error
- **Behavior**: Logs warning with error details
- **Impact**: Ingestion continues (score write is non-critical)
- **Log**: Warning with full error response

### No Matches Found
- **Behavior**: Logs info message, returns early
- **Impact**: No Notion update (idempotent)
- **Log**: Info message about no matches

## Future Enhancements

Potential improvements:
1. Store individual scores per signature (not just highest)
2. Add score history/tracking
3. Add score thresholds for filtering
4. Add score visualization in Notion

## Related Files

- `amprenta_rag/ingestion/signature_matching.py` - Main implementation
- `amprenta_rag/ingestion/dataset_ingestion.py` - Calls signature matching
- `amprenta_rag/signatures/signature_scoring.py` - Score calculation logic

## Summary

✅ **Feature Status**: Fully implemented and tested  
✅ **Test Coverage**: Both test datasets verified  
✅ **Error Handling**: Robust with graceful degradation  
✅ **Logging**: Comprehensive and clear  
✅ **Idempotency**: Maintained - safe to re-run  

The Signature Match Score feature is production-ready and automatically populates during dataset ingestion.

