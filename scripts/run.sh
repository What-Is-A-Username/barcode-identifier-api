#!/bin/sh 
# Script run to finish building the app container 

set -e 

# run commands needed to start server
python manage.py wait_for_db
# collect static files from all Django apps into STATIC_ROOT (--noinput prevents prompt)
python manage.py collectstatic --noinput

echo Deleting existing migrations ...
# Delete remaining migrations
rm -rf ./barcode_blastn/migrations/**
mkdir -p ./barcode_blastn/migrations
# Add __init__ file 
echo "" > ./barcode_blastn/migrations/__init__.py

echo Setting up database and creating necessary migrations ... 

python manage.py makemigrations

echo Running any migrations ...

python manage.py migrate

echo Creating superuser for $DJANGO_SUPERUSER_USERNAME with email $DJANGO_SUPERUSER_EMAIL and password $DJANGO_SUPERUSER_PASSWORD

python manage.py init_superuser --username $DJANGO_SUPERUSER_USERNAME --email $DJANGO_SUPERUSER_EMAIL --password $DJANGO_SUPERUSER_PASSWORD

echo "Setting up /var/data/library folder to store run and library data."
mkdir -p /var/data/library
chmod -R 764 /var/data/library
chown -R appuser:appgroup /var/data/library
echo "/var/data/library complete."

echo "Setting up /var/data/runs folder to store run and library data"
mkdir -p /var/data/runs
chmod -R 764 /var/data/runs
chown -R appuser:appgroup /var/data/runs
echo "/var/data/runs complete."

echo "Setting up /vol/static/runs folder to store run and library data"
mkdir -p /vol/static/runs
chmod -R 764 /vol/static/runs
chown -R appuser:appgroup /vol/static/runs
echo "/vol/static/runs complete."

# run on socket on port 9000 (nginx serves 9000)
# workers: # of workers working concurrently in the container accepting requests
# run as master daemon (foreground)
# enable multithreading
# run as wsgi service
uwsgi --socket :9000 --workers 1 --master --enable-threads --module barcode_identifier_api.wsgi
