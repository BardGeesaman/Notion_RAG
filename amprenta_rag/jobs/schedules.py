"""Celery Beat schedule configuration."""

from celery.schedules import crontab

# Static schedules for external source synchronization
CELERYBEAT_SCHEDULE = {
    'sync-chembl-daily': {
        'task': 'amprenta_rag.jobs.tasks.scheduled.run_external_sync',
        'schedule': crontab(hour=2, minute=0),  # 2:00 AM daily
        'args': ['chembl', 'incremental'],
        'options': {'queue': 'scheduled'},
    },
    'sync-pubchem-daily': {
        'task': 'amprenta_rag.jobs.tasks.scheduled.run_external_sync',
        'schedule': crontab(hour=3, minute=0),  # 3:00 AM daily
        'args': ['pubchem', 'incremental'],
        'options': {'queue': 'scheduled'},
    },
    # Weekly digest generation (Mondays at 8:00 AM)
    'weekly-digest-check': {
        'task': 'amprenta_rag.jobs.tasks.scheduled.check_digest_schedules',
        'schedule': crontab(hour=8, minute=0, day_of_week=1),  # Monday 8:00 AM
        'options': {'queue': 'scheduled'},
    },
    # Daily database backup (1:00 AM UTC)
    'daily-database-backup': {
        'task': 'amprenta_rag.jobs.tasks.backup.run_database_backup',
        'schedule': crontab(hour=1, minute=0),  # 1:00 AM daily
        'args': ['full'],
        'options': {'queue': 'scheduled'},
    },
    # Weekly backup cleanup (Sundays at 4:00 AM UTC)
    'weekly-backup-cleanup': {
        'task': 'amprenta_rag.jobs.tasks.backup.cleanup_old_backups',
        'schedule': crontab(hour=4, minute=0, day_of_week=0),  # Sunday 4:00 AM
        'options': {'queue': 'scheduled'},
    },
    # Weekly backup verification (Sundays at 5:00 AM UTC)
    'weekly-backup-verify': {
        'task': 'amprenta_rag.jobs.tasks.backup.verify_latest_backup',
        'schedule': crontab(hour=5, minute=0, day_of_week=0),  # Sunday 5:00 AM
        'options': {'queue': 'scheduled'},
    },
    # Daily backup health check (6:00 AM UTC)
    'daily-backup-health-check': {
        'task': 'amprenta_rag.jobs.tasks.backup.backup_health_check',
        'schedule': crontab(hour=6, minute=0),  # 6:00 AM daily
        'options': {'queue': 'scheduled'},
    },
    # Hourly cleanup of expired project exports
    'hourly-export-cleanup': {
        'task': 'amprenta_rag.jobs.tasks.backup.cleanup_expired_exports',
        'schedule': crontab(minute=0),  # Every hour at the top of the hour
        'options': {'queue': 'scheduled'},
    },
}

# Note: Dynamic harvest schedules are managed by database-driven periodic checks
# rather than static Celery Beat schedules, as they can be created/modified at runtime
