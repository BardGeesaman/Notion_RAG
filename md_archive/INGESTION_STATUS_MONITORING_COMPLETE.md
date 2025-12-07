# Ingestion Status & Back-Pressure Indicators v1 - Complete

**Date**: 2025-12-05  
**Status**: ✅ **COMPLETE**

## Summary

Ingestion status monitoring and pipeline health check scripts have been created and integrated into the deployment workflow. These tools provide visibility into ingestion back-pressure and system health.

---

## 1. Status Fields Audit

### Current Status Fields

**Dataset Model**:
- ❌ **No `ingestion_status` field** (inferred from feature count and timestamps)
- ❌ **No `embedding_status` field**
- ❌ **No `qc_status` field**
- ❌ **No `last_ingested_at` field**
- ✅ Has `created_at` and `updated_at` timestamps
- ✅ Has `features` relationship (can infer status from feature count)

**Program Model**:
- ❌ No status fields
- ✅ Has `created_at` and `updated_at` timestamps

**Experiment Model**:
- ❌ No status fields
- ✅ Has `created_at` and `updated_at` timestamps

**Literature Model**:
- ✅ Has `embedding_status` field
- ✅ Has `last_ingested_at` field

**Email Model**:
- ✅ Has `embedding_status` field
- ✅ Has `last_ingested_at` field

### Status Inference Logic

Since `Dataset` doesn't have explicit status fields, the monitoring script infers status from:
- **`complete`**: Has features (feature_count > 0)
- **`in_progress`**: Has data source but no features, created < 6 hours ago
- **`pending`**: Has data source but no features, created > 6 hours ago
- **`incomplete`**: No data source (missing file_paths, file_urls, external_ids)

---

## 2. Created Scripts

### `scripts/list_ingestion_status.py`

**Purpose**: Provides a comprehensive "back-pressure" snapshot of ingestion status.

**Features**:
- Lists all datasets with inferred status
- Groups by omics type and status
- Shows feature counts and last updated timestamps
- Identifies datasets stuck in "pending" for > 6 hours
- Reports on Programs, Experiments, Literature, and Emails
- Shows embedding status for Literature and Emails

**Output Example**:
```
Summary by Status:
  complete       :   6
  in_progress    :   2
  pending        :   5

⚠️  WARNING: 5 datasets stuck in 'pending' for > 6 hours
```

**Usage**:
```bash
python scripts/list_ingestion_status.py
```

### `scripts/check_pipeline_health.py`

**Purpose**: Quick smoke test to verify core pipeline functionality.

**Checks**:
1. ✅ **Database Connectivity**: Verifies Postgres connection and basic counts
2. ✅ **RAG Query Functionality**: Tests RAG query (if available)
3. ✅ **Stuck Ingestion Jobs**: Detects datasets stuck in ingestion > 6 hours

**Exit Codes**:
- `0`: All checks passed
- `1`: One or more checks failed

**Output Example**:
```
✅ PASS     Database Connectivity
✅ PASS     RAG Query
❌ FAIL     Stuck Jobs
```

**Usage**:
```bash
python scripts/check_pipeline_health.py
```

---

## 3. Documentation Updates

### `docs/DEPLOYMENT_GUIDE.md`

**Added Sections**:
1. **Ingestion Status Monitoring**: Instructions for using `list_ingestion_status.py`
2. **Pipeline Health Checks**: Instructions for using `check_pipeline_health.py`
3. **Monitoring Routine**: Recommended workflow after deployments
4. **Automated Monitoring**: Example cron job setup
5. **Regular Tasks**: Added ingestion status review to maintenance section

**Key Documentation Points**:
- Run `list_ingestion_status.py` after deployments or batch ingestion
- Run `check_pipeline_health.py` to verify core functionality
- Review for stuck datasets and backlog
- Set up automated monitoring with cron jobs

---

## 4. Current Status Snapshot

**From Test Run** (2025-12-05):

**Datasets**:
- Total: 13
- Complete: 6 (46.2%)
- In Progress: 2 (15.4%)
- Pending: 5 (38.5%)
- Stuck (> 6 hours): 5 datasets

**Breakdown by Omics Type**:
- **lipidomics**: 1 complete, 2 pending
- **metabolomics**: 2 complete, 1 pending
- **proteomics**: 1 complete, 2 pending
- **transcriptomics**: 2 complete, 1 in_progress

**Issues Identified**:
- 5 datasets stuck in "pending" for > 6 hours
  - 2 internal lipidomics datasets (27.7 hours old)
  - 1 MetaboLights metabolomics dataset (22.6 hours old)
  - 2 PRIDE proteomics datasets (21.2 hours old)

---

## 5. Future Enhancements

### Priority 1: Add Status Fields to Dataset Model

**Recommended Migration**:
```python
# Add to Dataset model:
ingestion_status = Column(String(50), nullable=True, default="pending")
embedding_status = Column(String(50), nullable=True, default="Not Embedded")
qc_status = Column(String(50), nullable=True)
last_ingested_at = Column(DateTime, nullable=True)
```

**Benefits**:
- Explicit status tracking (no inference needed)
- Better visibility into ingestion pipeline
- Ability to track QC status separately
- Clear timestamps for last ingestion

### Priority 2: Visual Dashboard

**Streamlit Dashboard Enhancement**:
- Add "Ingestion Status" page
- Show back-pressure indicators
- Visual status breakdown by omics type
- Alert badges for stuck datasets

### Priority 3: Automated Alerts

**Enhancement**:
- Email/Slack notifications for stuck jobs
- Configurable thresholds (e.g., alert if > 10 stuck datasets)
- Integration with monitoring systems (Prometheus, Grafana)

---

## 6. Usage Recommendations

### After Deployments

1. **Check ingestion status**:
   ```bash
   python scripts/list_ingestion_status.py
   ```
   - Review for backlog
   - Check for stuck datasets
   - Verify feature coverage

2. **Run health checks**:
   ```bash
   python scripts/check_pipeline_health.py
   ```
   - Verify core functionality
   - Confirm no stuck jobs
   - Check database connectivity

### Weekly Reviews

- Review `list_ingestion_status.py` output weekly
- Identify patterns in stuck datasets
- Track feature coverage trends
- Review literature/email embedding status

### Automated Monitoring

Set up cron job for hourly health checks:
```bash
# /etc/cron.hourly/amprenta-health-check
#!/bin/bash
cd /opt/amprenta-rag
source venv/bin/activate
python scripts/check_pipeline_health.py >> /var/log/amprenta/health.log 2>&1
```

---

## 7. Technical Notes

### Fixed Issues
- ✅ Fixed `PG_UUID` → `UUID` in `models.py` (feature_pathway_map table)
- ✅ Scripts work correctly with current model structure
- ✅ Status inference logic handles edge cases

### Known Limitations
- Status is inferred, not explicit (until Dataset model is updated)
- RAG query check may be skipped if Notion sync is disabled (expected)
- Stuck job detection uses 6-hour threshold (configurable in script)

---

## Conclusion

**Ingestion Status & Back-Pressure Indicators v1 is complete.** The monitoring scripts provide visibility into ingestion pipeline health and can be integrated into deployment workflows and automated monitoring systems.

**Next Steps**:
1. Add explicit status fields to Dataset model (when Implementor is ready)
2. Enhance Streamlit dashboard with visual status indicators
3. Set up automated monitoring with cron jobs
4. Consider adding email/Slack alerts for critical issues
