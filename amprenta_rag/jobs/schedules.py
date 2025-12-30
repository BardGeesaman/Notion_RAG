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
}

# Note: Dynamic harvest schedules are managed by database-driven periodic checks
# rather than static Celery Beat schedules, as they can be created/modified at runtime
