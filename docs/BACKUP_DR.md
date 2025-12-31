# Backup & Disaster Recovery Guide

## Overview

The Amprenta RAG platform implements a comprehensive backup and disaster recovery strategy designed to protect critical research data and ensure business continuity.

### Capabilities

- **Automated Daily Backups**: Full PostgreSQL database backups with configurable retention
- **Point-in-Time Recovery (PITR)**: WAL archiving for granular recovery options
- **Health Monitoring**: Continuous backup system monitoring with admin notifications
- **Multi-Storage Support**: S3-compatible storage and local filesystem options
- **Data Integrity**: Automated verification and checksum validation
- **Project Export**: Selective data export for migration and archival

### Recovery Objectives

- **Recovery Point Objective (RPO)**: ≤ 24 hours for standard backups, ≤ 15 minutes with WAL archiving
- **Recovery Time Objective (RTO)**: ≤ 4 hours for full system restore, ≤ 1 hour for project-level restore

## Daily Backup System

### Schedule Configuration

The backup system operates on the following automated schedule:

```python
# Configured in amprenta_rag/jobs/schedules.py
CELERYBEAT_SCHEDULE = {
    'daily-database-backup': {
        'task': 'amprenta_rag.jobs.tasks.backup.run_database_backup',
        'schedule': crontab(hour=1, minute=0),  # 1:00 AM UTC
        'args': ('full',),
        'options': {'queue': 'scheduled'}
    }
}
```

### Backup Types

**Full Backups (Default)**
- Complete PostgreSQL database dump using `pg_dump`
- Includes all schemas, data, indexes, and constraints
- Compressed and checksummed for integrity verification
- Average size: 500MB - 2GB (depending on data volume)

**Incremental Backups (Future Enhancement)**
- WAL-based incremental changes
- Reduced storage footprint and backup time
- Requires WAL archiving configuration

### Storage Options

**S3-Compatible Storage (Recommended)**
```python
# Configuration in amprenta_rag/config.py
class BackupConfig:
    storage_type: str = "s3"
    s3_bucket: str = "amprenta-backups"
    s3_prefix: str = "database/"
    s3_region: str = "us-east-1"
```

**Local Filesystem**
```python
class BackupConfig:
    storage_type: str = "local"
    local_path: str = "/var/backups/amprenta"
```

## Retention Policy

### Default Configuration

- **Retention Period**: 365 days (configurable via `BACKUP_RETENTION_DAYS`)
- **Cleanup Schedule**: Sundays at 4:00 AM UTC
- **Minimum Retention**: 7 days (safety constraint)

### Cleanup Process

```python
# Weekly cleanup task
'weekly-backup-cleanup': {
    'task': 'amprenta_rag.jobs.tasks.backup.cleanup_old_backups',
    'schedule': crontab(hour=4, minute=0, day_of_week=0),  # Sunday 4:00 AM
    'options': {'queue': 'scheduled'}
}
```

The cleanup process:
1. Identifies backups older than retention period
2. Removes backup files from storage
3. Deletes corresponding database records
4. Logs cleanup statistics

## Health Monitoring & Alerting

### Daily Health Checks

```python
# Health check schedule
'daily-backup-health-check': {
    'task': 'amprenta_rag.jobs.tasks.backup.backup_health_check',
    'schedule': crontab(hour=6, minute=0),  # 6:00 AM UTC
    'options': {'queue': 'scheduled'}
}
```

### Monitored Metrics

**System Health Indicators:**
- Recent backup availability (within 48 hours)
- Failed backup count (last 7 days)
- Total backup count and storage usage
- Backup engine availability

**Alert Conditions:**
- **Warning**: No backup in 48 hours, recent failures detected
- **Critical**: Backup creation failures, verification failures, system errors

### Admin Notifications

All backup alerts are automatically sent to admin users via the notification system:

```python
# Example notification
create_admin_notification(
    title="Backup System Warning",
    message="No recent backup found (48+ hours)",
    notification_type="warning"
)
```

## Point-in-Time Recovery (PITR) Setup

### PostgreSQL Configuration

For production deployments, enable WAL archiving to support point-in-time recovery:

#### On-Premises PostgreSQL

Add to `postgresql.conf`:
```conf
# Enable WAL archiving
wal_level = replica
archive_mode = on
archive_command = 'aws s3 cp %p s3://your-backup-bucket/wal/%f'
archive_timeout = 300  # Archive every 5 minutes

# Increase WAL retention for safety
max_wal_size = 2GB
min_wal_size = 1GB

# Enable checksums for data integrity
data_checksums = on
```

Add to `pg_hba.conf` for replication access:
```conf
# Allow replication connections
host replication backup_user 10.0.0.0/8 md5
```

#### Docker Compose Configuration

```yaml
services:
  postgres:
    image: postgres:15
    environment:
      - POSTGRES_DB=amprenta
      - POSTGRES_USER=postgres
      - POSTGRES_PASSWORD=${DB_PASSWORD}
    volumes:
      - postgres_data:/var/lib/postgresql/data
      - ./config/postgresql.conf:/etc/postgresql/postgresql.conf
    command: postgres -c config_file=/etc/postgresql/postgresql.conf
```

### AWS RDS Configuration

For AWS RDS deployments:

#### Parameter Group Settings

Create custom parameter group with:
```
log_statement = all
log_min_duration_statement = 1000
shared_preload_libraries = pg_stat_statements
max_connections = 200
```

#### Automated Backups

```python
# RDS backup configuration
{
    "BackupRetentionPeriod": 35,  # 35 days
    "BackupWindow": "03:00-04:00",  # UTC
    "PreferredMaintenanceWindow": "sun:04:00-sun:05:00",
    "DeleteAutomatedBackups": False,
    "DeletionProtection": True
}
```

#### Cross-Region Replication

```python
# Read replica configuration
{
    "SourceDBInstanceIdentifier": "amprenta-prod",
    "DBInstanceIdentifier": "amprenta-prod-replica",
    "DBInstanceClass": "db.t3.medium",
    "AvailabilityZone": "us-west-2a",
    "MultiAZ": False,
    "PubliclyAccessible": False
}
```

### WAL Archiving Script

Create a robust WAL archiving script:

```bash
#!/bin/bash
# /usr/local/bin/archive_wal.sh

WAL_FILE="$1"
WAL_PATH="$2"
S3_BUCKET="your-backup-bucket"
S3_PREFIX="wal"

# Upload to S3 with retry logic
aws s3 cp "$WAL_PATH" "s3://$S3_BUCKET/$S3_PREFIX/$WAL_FILE" \
    --storage-class STANDARD_IA \
    --metadata "archived=$(date -u +%Y-%m-%dT%H:%M:%SZ)"

# Verify upload
aws s3 ls "s3://$S3_BUCKET/$S3_PREFIX/$WAL_FILE" > /dev/null
if [ $? -eq 0 ]; then
    echo "WAL file $WAL_FILE archived successfully"
    exit 0
else
    echo "Failed to archive WAL file $WAL_FILE"
    exit 1
fi
```

## Restore Procedures

### 1. Full Database Restore

#### From Application Backup

```bash
# 1. Stop the application
docker-compose down

# 2. Download backup file
aws s3 cp s3://amprenta-backups/database/backup_20231201_010000.sql.gz ./

# 3. Restore database
gunzip backup_20231201_010000.sql.gz
psql -h localhost -U postgres -d amprenta < backup_20231201_010000.sql

# 4. Restart application
docker-compose up -d
```

#### Using Backup Engine

```python
from amprenta_rag.backup.backup_engine import BackupEngine

# Initialize engine
engine = BackupEngine()

# List available backups
backups = engine.list_backups(limit=10)

# Restore from specific backup
backup_id = backups[0].id
engine.restore_backup(backup_id, target_database="amprenta_restored")
```

### 2. Point-in-Time Recovery

```bash
# 1. Stop PostgreSQL
systemctl stop postgresql

# 2. Backup current data directory
mv /var/lib/postgresql/data /var/lib/postgresql/data.backup

# 3. Restore base backup
mkdir /var/lib/postgresql/data
tar -xzf base_backup_20231201.tar.gz -C /var/lib/postgresql/data

# 4. Create recovery configuration
cat > /var/lib/postgresql/data/recovery.conf << EOF
restore_command = 'aws s3 cp s3://your-backup-bucket/wal/%f %p'
recovery_target_time = '2023-12-01 15:30:00 UTC'
recovery_target_action = 'promote'
EOF

# 5. Start PostgreSQL (will enter recovery mode)
systemctl start postgresql

# 6. Monitor recovery progress
tail -f /var/log/postgresql/postgresql.log
```

### 3. Project-Level Export/Import

#### Export Project Data

```python
from amprenta_rag.api.routers.backup import create_project_export

# Create export
export_result = create_project_export(
    project_ids=["project-uuid"],
    include_data=True,
    include_files=True
)

# Download export
export_id = export_result["export_id"]
# Use GET /api/v1/backup/export/{export_id} to download
```

#### Import Project Data

```bash
# 1. Extract export ZIP
unzip project_export_20231201.zip

# 2. Run import script
python scripts/import_project.py \
    --export-dir ./project_export_20231201 \
    --target-database amprenta \
    --dry-run  # Remove for actual import
```

## Verification & Integrity

### Weekly Backup Verification

```python
# Verification schedule
'weekly-backup-verify': {
    'task': 'amprenta_rag.jobs.tasks.backup.verify_latest_backup',
    'schedule': crontab(hour=5, minute=0, day_of_week=0),  # Sunday 5:00 AM
    'options': {'queue': 'scheduled'}
}
```

### Verification Process

The verification task performs:
1. **File Integrity**: Checksum validation against stored checksums
2. **SQL Validation**: Basic SQL parsing of backup files
3. **Restore Test**: Periodic restore to temporary database (monthly)

### Manual Verification

```python
from amprenta_rag.backup.backup_engine import BackupEngine

# Verify specific backup
engine = BackupEngine()
is_valid = engine.verify_backup(backup_id)

if is_valid:
    print("Backup integrity verified")
else:
    print("Backup verification failed")
```

## Monitoring & Observability

### Backup Metrics

Monitor these key metrics:

```python
# Backup success rate (last 30 days)
SELECT 
    COUNT(*) FILTER (WHERE status = 'completed') * 100.0 / COUNT(*) as success_rate
FROM backup_records 
WHERE created_at >= NOW() - INTERVAL '30 days';

# Average backup size and duration
SELECT 
    AVG(file_size_bytes) as avg_size,
    AVG(EXTRACT(EPOCH FROM (completed_at - started_at))) as avg_duration_seconds
FROM backup_records 
WHERE status = 'completed' AND created_at >= NOW() - INTERVAL '7 days';

# Storage usage trend
SELECT 
    DATE(created_at) as backup_date,
    SUM(file_size_bytes) as total_size
FROM backup_records 
WHERE status = 'completed'
GROUP BY DATE(created_at)
ORDER BY backup_date DESC
LIMIT 30;
```

### Health Dashboard

Create monitoring dashboard with:
- Backup success/failure rates
- Storage usage trends
- Recovery time metrics
- Alert frequency

## Troubleshooting

### Common Issues

#### 1. Backup Failures

**Symptoms**: `BackupEngineError` in logs, failed backup status

**Diagnosis**:
```bash
# Check disk space
df -h /var/backups

# Check S3 credentials
aws s3 ls s3://your-backup-bucket/

# Check PostgreSQL connectivity
psql -h localhost -U postgres -c "SELECT version();"
```

**Solutions**:
- Ensure sufficient disk space (2x database size)
- Verify S3 credentials and bucket permissions
- Check PostgreSQL connection settings
- Review backup engine configuration

#### 2. WAL Archiving Issues

**Symptoms**: WAL files accumulating in `pg_wal/` directory

**Diagnosis**:
```sql
-- Check archiver status
SELECT * FROM pg_stat_archiver;

-- Check current WAL file
SELECT pg_current_wal_lsn();
```

**Solutions**:
- Verify `archive_command` script permissions
- Check S3 connectivity and credentials
- Ensure sufficient storage space
- Review PostgreSQL logs for archiver errors

#### 3. Restore Failures

**Symptoms**: Restore process hangs or fails with errors

**Diagnosis**:
```bash
# Check backup file integrity
gzip -t backup_file.sql.gz

# Verify PostgreSQL version compatibility
psql --version
```

**Solutions**:
- Ensure backup file is not corrupted
- Verify PostgreSQL version compatibility
- Check target database permissions
- Review restore logs for specific errors

### Alert Response Procedures

#### Backup System Warning

1. **Immediate Actions**:
   - Check backup system health dashboard
   - Review recent backup logs
   - Verify storage availability

2. **Investigation Steps**:
   - Check disk space on backup storage
   - Verify S3 credentials and permissions
   - Review PostgreSQL connectivity
   - Check for recent configuration changes

3. **Resolution**:
   - Address identified issues
   - Manually trigger backup if needed
   - Update monitoring if false positive

#### Backup System Critical

1. **Immediate Actions**:
   - Escalate to on-call engineer
   - Check system status dashboard
   - Verify database accessibility

2. **Emergency Procedures**:
   - Manual backup creation if automated system fails
   - Activate disaster recovery procedures if needed
   - Document incident for post-mortem

### Recovery Testing

#### Monthly Recovery Drills

Perform monthly recovery testing:

```bash
#!/bin/bash
# Monthly recovery test script

# 1. Create test environment
docker run --name postgres-test -e POSTGRES_PASSWORD=test -d postgres:15

# 2. Restore latest backup
LATEST_BACKUP=$(aws s3 ls s3://amprenta-backups/database/ | sort | tail -n 1 | awk '{print $4}')
aws s3 cp "s3://amprenta-backups/database/$LATEST_BACKUP" ./test_restore.sql.gz
gunzip test_restore.sql.gz

# 3. Perform restore
docker exec -i postgres-test psql -U postgres < test_restore.sql

# 4. Verify data integrity
docker exec postgres-test psql -U postgres -c "SELECT COUNT(*) FROM users;"

# 5. Cleanup
docker rm -f postgres-test
rm test_restore.sql
```

## Security Considerations

### Encryption

- **In-Transit**: All backup transfers use TLS/SSL encryption
- **At-Rest**: S3 server-side encryption (SSE-S3 or SSE-KMS)
- **Database**: Enable PostgreSQL TLS for backup connections

### Access Control

```python
# S3 bucket policy for backup access
{
    "Version": "2012-10-17",
    "Statement": [
        {
            "Sid": "BackupAccess",
            "Effect": "Allow",
            "Principal": {
                "AWS": "arn:aws:iam::ACCOUNT:role/AmprenterBackupRole"
            },
            "Action": [
                "s3:GetObject",
                "s3:PutObject",
                "s3:DeleteObject"
            ],
            "Resource": "arn:aws:s3:::amprenta-backups/*"
        }
    ]
}
```

### Compliance

- **Data Retention**: Configure retention policies per regulatory requirements
- **Audit Logging**: Enable CloudTrail for S3 access logging
- **Encryption Keys**: Use customer-managed KMS keys for sensitive data

## Configuration Reference

### Environment Variables

```bash
# Backup Configuration
BACKUP_STORAGE_TYPE=s3
BACKUP_S3_BUCKET=amprenta-backups
BACKUP_S3_PREFIX=database/
BACKUP_RETENTION_DAYS=365

# Database Configuration
DB_HOST=localhost
DB_PORT=5432
DB_NAME=amprenta
DB_USER=postgres
DB_PASSWORD=secure_password

# Monitoring
BACKUP_HEALTH_CHECK_ENABLED=true
BACKUP_NOTIFICATION_ENABLED=true
```

### Backup Engine Configuration

```python
from amprenta_rag.config import BackupConfig

config = BackupConfig(
    storage_type="s3",
    s3_bucket="amprenta-backups",
    s3_prefix="database/",
    retention_days=365,
    compression_enabled=True,
    checksum_enabled=True,
    notification_enabled=True
)
```

---

**Document Version**: 1.0  
**Last Updated**: December 2024  
**Next Review**: March 2025
