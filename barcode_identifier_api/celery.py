import os

os.environ.setdefault('DJANGO_SETTINGS_MODULE', 'barcode_identifier_api.settings')

from celery import Celery
from celery.schedules import crontab

# Set the default Django settings module for the 'celery' program.

app = Celery('barcode_identifier_api')

# Using a string here means the worker doesn't have to serialize
# the configuration object to child processes.
# - namespace='CELERY' means all celery-related configuration keys
#   should have a `CELERY_` prefix.
app.config_from_object('django.conf:settings', namespace='CELERY')

# Load task modules from all registered Django apps.
app.autodiscover_tasks()

app.conf.beat_schedule = {
    'database-scheduler': {
        'task': 'barcode_blastn.tasks.update_database',
        # 'schedule': crontab(minute=0, hour=0, day_of_week='sunday'), # update all databases Saturday midnight
        'schedule': crontab(minute='0', hour='0'), # update daily at midnight 
    }
}
