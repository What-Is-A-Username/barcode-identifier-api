#!/bin/sh 
# Script run to finish building the app container 

set -e 

# run commands needed to start server
python manage.py wait_for_db
python manage.py migrate 

celery -A barcode_identifier_api worker --loglevel=INFO -B -c 1 -Q BarcodeQueue.fifo -s /var/data/celerybeat-schedule
chown nobody:nogroup "/var/data/celerybeat-schedule"
