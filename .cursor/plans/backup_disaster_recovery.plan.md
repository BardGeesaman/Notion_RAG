# Automated Backup & Disaster Recovery

## Overview

Implement automated database backups to S3, Point-in-Time Recovery (via RDS), and user-initiated project exports. Leverages the newly implemented Celery job queue for scheduled backups.

## Current State

- PostgreSQL is sole database (via `amprenta_rag/database/base.py`)
- AWS infrastructure exists in `deploy/aws/terraform/` but no S3 bucket
- Celery job queue implemented in `amprenta_rag/jobs/` - can schedule backup tasks
- Export functionality exists in `amprenta_rag/export/export_engine.py`
- No backup-specific code exists
- RDS has no backup settings configured (P1 issue from Reviewer)

## Architecture

```
┌─────────────────────────────────────────────────────────────────┐
│                    SCHEDULED BACKUPS                             │
├─────────────────────────────────────────────────────────────────┤
│  Celery Beat ──► backup_database task ──► pg_dump               │
│                         │                                        │
│                         ▼                                        │
│                   gzip compress ──► S3 Upload (SSE-KMS)         │
└─────────────────────────────────────────────────────────────────┘

┌─────────────────────────────────────────────────────────────────┐
│                    RDS AUTOMATED BACKUPS                         │
├─────────────────────────────────────────────────────────────────┤
│  AWS RDS ──► backup_retention_period=7 ──► PITR enabled         │
│              backup_window="03:00-04:00"                         │
└─────────────────────────────────────────────────────────────────┘

┌─────────────────────────────────────────────────────────────────┐
│                    USER EXPORT                                   │
├─────────────────────────────────────────────────────────────────┤
│  POST /api/v1/backup/export ──► PackageBuilder ──► ZIP+manifest │
│                                                    │             │
│                                                    ▼             │
│                                            Streaming download    │
└─────────────────────────────────────────────────────────────────┘

┌─────────────────────────────────────────────────────────────────┐
│                    ADMIN DASHBOARD                               │
├─────────────────────────────────────────────────────────────────┤
│  backup_admin.py: Manual backup | History | Export | Restore    │
└─────────────────────────────────────────────────────────────────┘
```

## Implementation Plan

### Batch 1: S3 Infrastructure and Client (7 tests)

**Terraform - S3 Bucket** (`deploy/aws/terraform/s3.tf`):
- S3 bucket with SSE-KMS encryption
- KMS key with rotation enabled
- Lifecycle: Glacier after 30 days, expire after 365 days
- IAM policy for ECS task access
- Block public access

**Terraform - RDS Backup** (update `deploy/aws/terraform/rds.tf`):
- Add `backup_retention_period = 7`
- Add `backup_window = "03:00-04:00"`

**Terraform - Variables** (update `deploy/aws/terraform/variables.tf`):
- Add `backup_retention_period` variable (default: 7)
- Add `backup_retention_days` variable (default: 365)

**Python Config** (update `amprenta_rag/config.py`):
- Add `BackupConfig` dataclass
- Env vars: `BACKUP_S3_BUCKET`, `BACKUP_S3_ENABLED`, `BACKUP_RETENTION_DAYS`, `BACKUP_LOCAL_DIR`

**S3 Client** (`amprenta_rag/backup/s3_client.py`):
- `upload_file(key, data, metadata)` - with SSE-KMS
- `download_file(key)` - returns bytes
- `list_backups(prefix)` - returns list of backup keys
- `delete_backup(key)` - removes backup
- Local filesystem fallback when `BACKUP_S3_ENABLED=false`

**Tests** (`amprenta_rag/tests/backup/test_s3_client.py`):
- 7 tests: upload, download, list, delete, local fallback, config loading, error handling

### Batch 2: Database Backup Engine (10 tests)

**Database Model** (update `amprenta_rag/database/models.py`):
- Add `BackupRecord` model:
  - id, backup_type (full/incremental), status (pending/running/completed/failed)
  - file_path, file_size_bytes, checksum_sha256
  - started_at, completed_at, error_message
  - metadata (JSON: tables, row_counts, pg_version)

**Migration**:
- Alembic migration for BackupRecord table

**Backup Engine** (`amprenta_rag/backup/backup_engine.py`):
- `create_backup(db, backup_type)` - pg_dump to compressed file, upload to S3
- `restore_backup(backup_id)` - download from S3, pg_restore
- `verify_backup(backup_id)` - checksum validation
- `generate_manifest(backup_id)` - JSON with tables, row counts, checksums
- `get_backup_status(backup_id)` - return BackupRecord

**Tests** (`amprenta_rag/tests/backup/test_backup_engine.py`):
- 10 tests: create, restore, verify, manifest, status, error handling, compression

### Batch 3: Celery Scheduled Backup Task (8 tests)

**Celery Task** (`amprenta_rag/jobs/tasks/backup.py`):
- `run_database_backup()` - full backup task with progress tracking
- `cleanup_old_backups()` - retention policy enforcement
- `verify_latest_backup()` - periodic verification task

**Schedule Configuration** (update `amprenta_rag/jobs/schedules.py`):
- Add daily full backup at 1:00 AM UTC
- Add weekly cleanup task (Sundays at 4:00 AM)

**Tests** (`amprenta_rag/tests/jobs/test_backup_tasks.py`):
- 8 tests: scheduled backup, cleanup, retry on failure, progress tracking

### Batch 4: API Endpoints and Project Export (9 tests)

**Project Export Service** (`amprenta_rag/backup/project_export.py`):
- `export_project(program_ids, experiment_ids, compound_ids)` - selective export
- Includes related datasets, features, signatures
- ZIP with manifest (SHA256 checksums)
- Reuse patterns from `amprenta_rag/export/export_engine.py`

**API Router** (`amprenta_rag/api/routers/backup.py`):
- `POST /api/v1/backup/database` - Trigger manual backup (admin only)
- `GET /api/v1/backup/history` - List backup history with pagination
- `GET /api/v1/backup/{id}` - Get backup details
- `GET /api/v1/backup/{id}/download` - Download backup file (streaming)
- `POST /api/v1/backup/export` - Create project export ZIP
- `GET /api/v1/backup/export/{id}` - Download export ZIP

**Register Router** (update `amprenta_rag/api/main.py`):
- Add backup router to app

**Tests** (`amprenta_rag/tests/api/test_backup_api.py`):
- 9 tests: all endpoints, auth, streaming, error handling

### Batch 5: Dashboard UI (9 tests)

**Dashboard Page** (`scripts/dashboard/pages/backup_admin.py`):
- Tab 1: Manual Backup (trigger button, progress indicator, download link)
- Tab 2: Scheduled Backups (enable/disable toggle, schedule display, history table)
- Tab 3: Project Export (entity selector, generate button, download)
- Tab 4: Restore Guide (step-by-step instructions, disaster recovery checklist)

**Register Page** (update `scripts/dashboard/core/config.py`):
- Add "Backup Admin" to PAGE_REGISTRY under Admin category

**Tests** (`amprenta_rag/tests/e2e/test_backup_admin_page.py`):
- 9 tests: page load, all tabs, manual backup trigger, export workflow

## Key Files Summary

| File | Purpose |
|------|---------|
| `deploy/aws/terraform/s3.tf` | S3 bucket + KMS for backups |
| `deploy/aws/terraform/rds.tf` | RDS backup settings (P1 fix) |
| `amprenta_rag/backup/__init__.py` | Package init |
| `amprenta_rag/backup/s3_client.py` | S3 upload/download with fallback |
| `amprenta_rag/backup/backup_engine.py` | pg_dump/restore logic |
| `amprenta_rag/backup/project_export.py` | User export builder |
| `amprenta_rag/jobs/tasks/backup.py` | Celery backup tasks |
| `amprenta_rag/api/routers/backup.py` | REST endpoints |
| `scripts/dashboard/pages/backup_admin.py` | Admin UI |

## Test Summary

| Batch | Tests | Coverage |
|-------|-------|----------|
| 1 - S3 Infrastructure | 7 | S3 client, config, fallback |
| 2 - Backup Engine | 10 | pg_dump, restore, verify, manifest |
| 3 - Celery Tasks | 8 | Scheduled backup, cleanup, retry |
| 4 - API & Export | 9 | Endpoints, auth, streaming |
| 5 - Dashboard | 9 | E2E page tests |
| **Total** | **43** | |

## Dependencies

Add to `requirements.txt`:
- `boto3>=1.34.0` (AWS SDK)

Already have: `celery[redis]`, `psycopg2-binary`

## Feature Flags

- `BACKUP_S3_ENABLED=true` - Use S3 storage (false = local filesystem)
- `BACKUP_S3_BUCKET` - S3 bucket name
- `BACKUP_LOCAL_DIR` - Local backup directory (default: `./backups`)
- `BACKUP_RETENTION_DAYS` - Days to retain backups (default: 365)

## P1 Issues (from Reviewer)

1. **RDS Automated Backup** - Add `backup_retention_period` and `backup_window` to rds.tf
2. **Encryption Key Management** - Use AWS KMS (SSE-KMS) instead of custom encryption

## P2 Deferred (Phase 2)

- Backup failure notifications (SNS/email)
- Automated restore verification to staging
- Cross-region S3 replication for DR

