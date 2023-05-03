#!/bin/sh 
# Script run to finish building the app container 

set -e 

# run commands needed to start server
python manage.py wait_for_db
# collect static files from all Django apps into STATIC_ROOT (--noinput prevents prompt)
python manage.py collectstatic --noinput

echo Creating superuser with $DJANGO_SUPERUSER_NAME $DJANGO_SUPERUSER_EMAIL $DJANGO_SUPERUSER_PASSWORD

python manage.py createsuperuser --noinput

# run on socket on port 9000 (nginx serves 9000)
# workers: # of workers working concurrently in the container accepting requests
# run as master daemon (foreground)
# enable multithreading
# run as wsgi service
uwsgi --socket :9000 --workers 1 --master --enable-threads --module barcode_identifier_api.wsgi
