"""Celery tasks for database backup operations."""

import logging
from datetime import datetime, timedelta, timezone
from typing import Dict, Optional

from amprenta_rag.backup.backup_engine import BackupEngine, BackupEngineError
from amprenta_rag.config import get_config
from amprenta_rag.database.models import BackupRecord, ProjectExport
from amprenta_rag.database.session import db_session
from amprenta_rag.jobs.celery_app import celery_app
from amprenta_rag.notifications.service import create_admin_notification

logger = logging.getLogger(__name__)


@celery_app.task(bind=True, max_retries=2, default_retry_delay=300, queue='scheduled')
def run_database_backup(self, backup_type: str = "full") -> Dict[str, str]:
    """Execute database backup with progress tracking.
    
    Args:
        backup_type: Type of backup to create ('full' or 'incremental')
        
    Returns:
        Dict with status, backup_id, and any error information
    """
    try:
        logger.info(f"Starting {backup_type} database backup")
        
        # Initialize backup engine
        engine = BackupEngine()
        
        # Create backup
        backup_id = engine.create_backup(backup_type)
        
        # Verify backup was created successfully
        backup_record = engine.get_backup_status(backup_id)
        if not backup_record:
            raise BackupEngineError(f"Backup record {backup_id} not found after creation")
        
        if backup_record.status == "completed":
            logger.info(f"Database backup {backup_id} completed successfully")
            return {
                "status": "completed",
                "backup_id": str(backup_id),
                "backup_type": backup_type,
                "file_size_bytes": str(backup_record.file_size_bytes or 0),
            }
        else:
            error_msg = backup_record.error_message or "Unknown error"
            raise BackupEngineError(f"Backup failed: {error_msg}")
            
    except Exception as exc:
        logger.exception(f"Database backup failed: {exc}")
        
        # Retry on transient failures
        if self.request.retries < self.max_retries:
            logger.info(f"Retrying backup (attempt {self.request.retries + 1}/{self.max_retries})")
            raise self.retry(exc=exc, countdown=300 * (2 ** self.request.retries))
        
        # Send notification to admin users after all retries exhausted
        create_admin_notification(
            title="Database Backup Failed",
            message=f"Database backup failed after {self.max_retries} retries: {str(exc)}",
            notification_type="critical",
        )
        
        return {
            "status": "failed",
            "backup_type": backup_type,
            "error": str(exc),
        }


@celery_app.task(queue='scheduled')
def cleanup_old_backups(retention_days: Optional[int] = None) -> Dict[str, str]:
    """Clean up old backup records and files based on retention policy.
    
    Args:
        retention_days: Number of days to retain backups (uses config default if None)
        
    Returns:
        Dict with cleanup results
    """
    try:
        if retention_days is None:
            retention_days = get_config().backup.retention_days
        
        cutoff_date = datetime.now(timezone.utc) - timedelta(days=retention_days)
        
        logger.info(f"Cleaning up backups older than {retention_days} days (before {cutoff_date})")
        
        deleted_count = 0
        error_count = 0
        
        with db_session() as db:
            # Find old backup records
            old_backups = db.query(BackupRecord).filter(
                BackupRecord.created_at < cutoff_date
            ).all()
            
            logger.info(f"Found {len(old_backups)} old backup records to clean up")
            
            # Initialize backup engine for file deletion
            engine = BackupEngine()
            
            for backup in old_backups:
                try:
                    # Delete backup file from storage if it exists
                    if backup.file_path:
                        try:
                            engine.backup_client.delete_backup(backup.file_path)
                            logger.debug(f"Deleted backup file: {backup.file_path}")
                        except Exception as e:
                            logger.warning(f"Failed to delete backup file {backup.file_path}: {e}")
                            error_count += 1
                    
                    # Delete database record
                    db.delete(backup)
                    deleted_count += 1
                    
                except Exception as e:
                    logger.error(f"Failed to delete backup record {backup.id}: {e}")
                    error_count += 1
            
            db.commit()
        
        logger.info(f"Cleanup completed: {deleted_count} backups deleted, {error_count} errors")
        
        return {
            "status": "completed",
            "deleted_count": str(deleted_count),
            "error_count": str(error_count),
            "retention_days": str(retention_days),
        }
        
    except Exception as exc:
        logger.exception(f"Backup cleanup failed: {exc}")
        return {
            "status": "failed",
            "error": str(exc),
        }


@celery_app.task(queue='scheduled')
def verify_latest_backup() -> Dict[str, str]:
    """Verify the integrity of the most recent backup.
    
    Returns:
        Dict with verification results
    """
    try:
        logger.info("Starting backup verification")
        
        # Find the most recent completed backup
        with db_session() as db:
            latest_backup = db.query(BackupRecord).filter(
                BackupRecord.status == "completed"
            ).order_by(BackupRecord.created_at.desc()).first()
            
            if not latest_backup:
                logger.warning("No completed backups found for verification")
                return {
                    "status": "skipped",
                    "reason": "No completed backups found",
                }
            
            backup_id = latest_backup.id
            db.expunge(latest_backup)  # Detach from session
        
        # Verify the backup
        engine = BackupEngine()
        is_valid = engine.verify_backup(backup_id)
        
        if is_valid:
            logger.info(f"Backup {backup_id} verification successful")
            return {
                "status": "verified",
                "backup_id": str(backup_id),
                "backup_date": latest_backup.created_at.isoformat(),
            }
        else:
            logger.error(f"Backup {backup_id} verification failed")
            
            # Send notification to admin users
            create_admin_notification(
                title="Backup Verification Failed",
                message=f"Backup verification failed for backup {backup_id} created on {latest_backup.created_at.isoformat()}",
                notification_type="critical",
            )
            
            return {
                "status": "failed",
                "backup_id": str(backup_id),
                "error": "Backup verification failed",
            }
            
    except Exception as exc:
        logger.exception(f"Backup verification failed: {exc}")
        
        # Send notification to admin users
        create_admin_notification(
            title="Backup Verification Error",
            message=f"Backup verification encountered an error: {str(exc)}",
            notification_type="critical",
        )
        
        return {
            "status": "failed",
            "error": str(exc),
        }


@celery_app.task(queue='scheduled')
def backup_health_check() -> Dict[str, str]:
    """Check backup system health and recent backup status.
    
    Returns:
        Dict with health check results
    """
    try:
        logger.info("Running backup system health check")
        
        health_info = {
            "status": "healthy",
            "checks": [],
        }
        
        with db_session() as db:
            # Check for recent backups (within last 48 hours)
            recent_cutoff = datetime.now(timezone.utc) - timedelta(hours=48)
            recent_backup = db.query(BackupRecord).filter(
                BackupRecord.created_at >= recent_cutoff,
                BackupRecord.status == "completed"
            ).first()
            
            if recent_backup:
                health_info["checks"].append("Recent backup found")
            else:
                health_info["status"] = "warning"
                health_info["checks"].append("No recent backup found (48h)")
            
            # Check for failed backups in last week
            week_cutoff = datetime.now(timezone.utc) - timedelta(days=7)
            failed_count = db.query(BackupRecord).filter(
                BackupRecord.created_at >= week_cutoff,
                BackupRecord.status == "failed"
            ).count()
            
            if failed_count == 0:
                health_info["checks"].append("No recent failures")
            else:
                health_info["status"] = "warning"
                health_info["checks"].append(f"{failed_count} failed backups in last week")
            
            # Check total backup count
            total_backups = db.query(BackupRecord).count()
            health_info["total_backups"] = str(total_backups)
            health_info["checks"].append(f"Total backups: {total_backups}")
        
        # Check backup engine availability
        try:
            BackupEngine()
            health_info["checks"].append("Backup engine available")
        except Exception as e:
            health_info["status"] = "error"
            health_info["checks"].append(f"Backup engine error: {e}")
        
        logger.info(f"Backup health check completed: {health_info['status']}")
        
        # Send notifications for warnings and errors
        if health_info["status"] == "warning":
            create_admin_notification(
                title="Backup System Warning",
                message=f"Backup system warning detected: {', '.join(health_info['checks'])}",
                notification_type="warning",
            )
        
        if health_info["status"] == "error":
            create_admin_notification(
                title="Backup System Error",
                message=f"Backup system error: {health_info.get('error', 'Unknown error')}",
                notification_type="critical",
            )
        
        return health_info
        
    except Exception as exc:
        logger.exception(f"Backup health check failed: {exc}")
        return {
            "status": "error",
            "error": str(exc),
        }


@celery_app.task(queue='scheduled')
def cleanup_expired_exports() -> Dict[str, str]:
    """Clean up expired project exports."""
    try:
        logger.info("Starting cleanup of expired project exports")
        
        current_time = datetime.now(timezone.utc)
        deleted_count = 0
        error_count = 0
        
        with db_session() as db:
            # Find expired exports
            expired_exports = db.query(ProjectExport).filter(
                ProjectExport.expires_at < current_time
            ).all()
            
            logger.info(f"Found {len(expired_exports)} expired exports to clean up")
            
            # Initialize backup engine for file deletion
            engine = BackupEngine()
            
            for export in expired_exports:
                try:
                    # Delete export file from storage if it exists
                    if export.file_path:
                        try:
                            engine.backup_client.delete_backup(export.file_path)
                            logger.debug(f"Deleted export file: {export.file_path}")
                        except Exception as e:
                            logger.warning(f"Failed to delete export file {export.file_path}: {e}")
                            error_count += 1
                    
                    # Delete database record
                    db.delete(export)
                    deleted_count += 1
                    
                except Exception as e:
                    logger.error(f"Failed to delete export record {export.id}: {e}")
                    error_count += 1
            
            db.commit()
        
        logger.info(f"Export cleanup completed: {deleted_count} exports deleted, {error_count} errors")
        
        return {
            "status": "completed",
            "deleted_count": str(deleted_count),
            "error_count": str(error_count),
        }
        
    except Exception as exc:
        logger.exception(f"Export cleanup failed: {exc}")
        return {
            "status": "failed",
            "error": str(exc),
        }
