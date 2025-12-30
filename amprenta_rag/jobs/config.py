"""Celery configuration."""

import os
from kombu import Queue

# Broker and result backend
broker_url = os.getenv('CELERY_BROKER_URL', 'redis://localhost:6379/0')
result_backend = os.getenv('CELERY_RESULT_BACKEND', 'redis://localhost:6379/0')

# Serialization (JSON only for security)
task_serializer = 'json'
result_serializer = 'json'
accept_content = ['json']

# Result settings
result_expires = int(os.getenv('CELERY_RESULT_EXPIRES', '86400'))  # 24 hours

# Worker settings
worker_prefetch_multiplier = 1  # For long-running tasks
worker_concurrency = int(os.getenv('CELERY_WORKER_CONCURRENCY', '4'))

# Task routing
task_default_queue = 'default'
task_routes = {
    'amprenta_rag.jobs.tasks.genomics.*': {'queue': 'default'},
    'amprenta_rag.jobs.tasks.docking.*': {'queue': 'high'},
    'amprenta_rag.jobs.tasks.extraction.*': {'queue': 'high'},
    'amprenta_rag.jobs.tasks.sync.*': {'queue': 'low'},
    'amprenta_rag.jobs.tasks.single_cell.*': {'queue': 'default'},
}

# Queue definitions
task_queues = (
    Queue('default', routing_key='default'),
    Queue('high', routing_key='high'),
    Queue('low', routing_key='low'),
    Queue('scheduled', routing_key='scheduled'),
)

# Testing mode
task_always_eager = os.getenv('CELERY_TASK_ALWAYS_EAGER', 'false').lower() == 'true'

# Celery Beat schedule configuration
from amprenta_rag.jobs.schedules import CELERYBEAT_SCHEDULE

beat_schedule = CELERYBEAT_SCHEDULE
beat_scheduler = 'celery.beat:PersistentScheduler'
