#!/bin/sh
# docker_recreate_db.sh: Drop the existing database in the docker postgres container
# and remake it. Requires that the container is already running. This method should
# be used as a last resort, since data is wiped, or a convenience script if you don't mind
# losing data currently stored.

docker exec -it barcode-identifier-api-db-1 psql -U admin -d postgres -c "DROP DATABASE barcode_identifier_db;"
docker exec -it barcode-identifier-api-db-1 psql -U admin -d postgres -c "CREATE DATABASE barcode_identifier_db WITH OWNER admin;"