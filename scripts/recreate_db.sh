#!bin/sh
# recreate_db.sh: Drop the existing database and remake it.
# This method should be used as a last resort, since data is wiped, or a convenience
# script if you don't mind losing data currently stored.

sudo su postgres
psql
drop database barcode_identifier_db;
create database barcode_identifier_db with owner <user_name>;
\q
exit
python manage.py makemigrations
python manage.py migrate
psql -h 127.0.0.1 -d barcode_identifier_db -U admin -f db.sql